// The memory headroom under the cgroup limits of this process (Slurm --mem, a container limit), for the window admission
// of diaspextract.cpp. Header-only and OpenMS-free so tests/test_cgroup_mem.cpp can run it on a fake /proc and /sys tree.
//
// headroom = limit - (usage - reclaimable page cache), the least over the leaf cgroup and its ancestors (a job's limit
// sits above its step's cgroup). The page cache a job's reads leave behind is charged to its cgroup (memory.current,
// memory.usage_in_bytes) but the kernel reclaims it before it kills anything, as MemAvailable assumes; without
// subtracting it, a job that has streamed its .d would see no headroom and run windows one at a time.
// The hierarchies are found in /proc/self/mountinfo, so a v1 memory controller mounted together with others
// ("cpu,memory"), a hybrid layout (v2 at /sys/fs/cgroup/unified) and a container's mount root are all handled; where
// mountinfo cannot be read, the standard mounts (/sys/fs/cgroup, /sys/fs/cgroup/memory) are assumed.
#pragma once

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace diaspextract
{
  namespace cgroup
  {
    /// One mounted cgroup hierarchy that carries the memory controller.
    struct MemHierarchy
    {
      bool v2 = false;
      std::string mount;   ///< where it is mounted (/sys/fs/cgroup/memory)
      std::string root;    ///< the cgroup path at that mount point ("/" outside a container)
    };

    inline bool hasToken_(const std::string& list, const std::string& tok)
    {
      std::stringstream ss(list);
      for (std::string t; std::getline(ss, t, ',');) if (t == tok) return true;
      return false;
    }

    /// The memory hierarchies listed in a /proc/self/mountinfo text: every cgroup2 mount and every cgroup (v1) mount
    /// whose super options name "memory".
    inline std::vector<MemHierarchy> parseMountinfo(std::istream& in)
    {
      std::vector<MemHierarchy> out;
      for (std::string line; std::getline(in, line);)
      {
        const size_t dash = line.find(" - ");
        if (dash == std::string::npos) continue;
        std::istringstream pre(line.substr(0, dash)), post(line.substr(dash + 3));
        std::string id, parent, dev, root, mount, fstype, source, super;
        if (!(pre >> id >> parent >> dev >> root >> mount) || !(post >> fstype >> source >> super)) continue;
        if (fstype == "cgroup2") out.push_back({true, mount, root});
        else if (fstype == "cgroup" && hasToken_(super, "memory")) out.push_back({false, mount, root});
      }
      return out;
    }

    /// The reclaimable page cache in a memory.stat text: active and inactive file pages (v1: the hierarchical
    /// total_* counters where present). 0 where absent.
    inline long long parseReclaimable(std::istream& in, bool v2)
    {
      long long local = 0, total = 0; bool have_total = false;
      std::string key; long long v;
      while (in >> key >> v)
      {
        if (key == "inactive_file" || key == "active_file") local += v;
        else if (!v2 && (key == "total_inactive_file" || key == "total_active_file")) { total += v; have_total = true; }
      }
      return have_total ? total : local;
    }

    /// A number file (memory.max, memory.current): its value, or -1 for "max", empty or unreadable.
    inline long long readNum_(const std::string& p)
    {
      std::ifstream f(p); std::string s;
      return (f >> s) && !s.empty() && std::isdigit((unsigned char)s[0]) ? std::atoll(s.c_str()) : -1;
    }

    /// The headroom in bytes, SIZE_MAX where no limit applies or nothing is readable. @p sysroot prefixes every path
    /// (empty on a real system; a fake tree in the tests).
    inline size_t headroom(const std::string& sysroot = "")
    {
      std::vector<MemHierarchy> hs;
      { std::ifstream mi(sysroot + "/proc/self/mountinfo"); if (mi) hs = parseMountinfo(mi); }
      if (hs.empty()) hs = {{true, "/sys/fs/cgroup", "/"}, {false, "/sys/fs/cgroup/memory", "/"}};
      size_t room = SIZE_MAX;
      std::ifstream cg(sysroot + "/proc/self/cgroup");
      for (std::string line; std::getline(cg, line);)
      {
        const size_t a = line.find(':'), b = line.find(':', a == std::string::npos ? a : a + 1);
        if (a == std::string::npos || b == std::string::npos) continue;
        const std::string ctl = line.substr(a + 1, b - a - 1);
        const bool v2 = line.compare(0, 2, "0:") == 0 && ctl.empty();
        if (!v2 && !hasToken_(ctl, "memory")) continue;
        const std::string path = line.substr(b + 1);
        for (const MemHierarchy& h : hs)
        {
          if (h.v2 != v2) continue;
          // the cgroup path as seen from this mount: strip the mount's root (a container's /docker/<id>)
          std::string rel = path;
          if (h.root != "/" && rel.compare(0, h.root.size(), h.root) == 0) rel = rel.substr(h.root.size());
          else if (h.root != "/") continue;   // this mount does not show our cgroup
          const char* lim = v2 ? "/memory.max" : "/memory.limit_in_bytes";
          const char* use = v2 ? "/memory.current" : "/memory.usage_in_bytes";
          for (std::string p = rel;; p = p.substr(0, p.rfind('/')))
          {
            const std::string dir = sysroot + h.mount + (p == "/" ? "" : p);
            const long long l = readNum_(dir + lim), u = readNum_(dir + use);
            if (l > 0 && u >= 0 && l < (1LL << 60))
            {
              long long cache = 0;
              { std::ifstream st(dir + "/memory.stat"); if (st) cache = parseReclaimable(st, v2); }
              const long long used = std::max(0LL, u - std::min(u, cache));
              room = std::min(room, (size_t)std::max(0LL, l - used));
            }
            if (p.empty() || p == "/" || p.rfind('/') == std::string::npos) break;
          }
        }
      }
      return room;
    }
  }
}
