// Tests of src/CgroupMem.h on fake /proc and /sys trees in a temporary directory: page cache is not counted as used,
// v1 controllers mounted together ("cpu,memory"), a hybrid layout (v2 at /sys/fs/cgroup/unified), a container mount
// root, the ancestor walk and the fallback mounts. Needs neither OpenMS nor Linux.
#include "CgroupMem.h"
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;
namespace cg = diaspextract::cgroup;

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { std::fprintf(stderr, "FAIL line %d: ", __LINE__); std::fprintf(stderr, __VA_ARGS__); std::fputc('\n', stderr); ++g_fail; } } while (0)

static void put(const fs::path& p, const std::string& text)
{
  fs::create_directories(p.parent_path());
  std::ofstream(p) << text;
}

static const long long G = 1LL << 30;

int main()
{
  const fs::path base = fs::temp_directory_path() / "spx_cgroup_mem_test";
  fs::remove_all(base);

  { // mountinfo parsing: v2, combined v1 cpu,memory, a v1 hierarchy without memory, a non-cgroup mount
    std::istringstream mi(
      "30 23 0:26 / /sys/fs/cgroup/unified rw,nosuid shared:4 - cgroup2 cgroup2 rw,nsdelegate\n"
      "31 23 0:27 / /sys/fs/cgroup/cpu,memory rw,nosuid shared:5 - cgroup cgroup rw,cpu,memory\n"
      "32 23 0:28 / /sys/fs/cgroup/pids rw,nosuid shared:6 - cgroup cgroup rw,pids\n"
      "33 1 8:1 / / rw,relatime shared:1 - ext4 /dev/sda1 rw\n");
    const auto hs = cg::parseMountinfo(mi);
    CHECK(hs.size() == 2, "%zu hierarchies", hs.size());
    CHECK(hs.size() == 2 && hs[0].v2 && hs[0].mount == "/sys/fs/cgroup/unified", "v2 mount");
    CHECK(hs.size() == 2 && !hs[1].v2 && hs[1].mount == "/sys/fs/cgroup/cpu,memory", "v1 mount");
  }
  { // memory.stat: v2 counts active+inactive file; v1 prefers the hierarchical totals
    std::istringstream v2("anon 100\nfile 900\nactive_file 300\ninactive_file 500\n");
    CHECK(cg::parseReclaimable(v2, true) == 800, "v2 reclaimable");
    std::istringstream v1("cache 900\nactive_file 1\ninactive_file 2\ntotal_active_file 300\ntotal_inactive_file 500\n");
    CHECK(cg::parseReclaimable(v1, false) == 800, "v1 reclaimable");
  }
  { // v2 with page cache: limit 64 G, usage 63 G of which 40 G file cache -> 41 G headroom; the parent's limit is lower
    const fs::path r = base / "v2";
    put(r / "proc/self/mountinfo", "30 23 0:26 / /sys/fs/cgroup rw - cgroup2 cgroup2 rw\n");
    put(r / "proc/self/cgroup", "0::/slurm/job1/step0\n");
    const fs::path j = r / "sys/fs/cgroup/slurm/job1";
    put(j / "memory.max", std::to_string(64 * G) + "\n");
    put(j / "memory.current", std::to_string(63 * G) + "\n");
    put(j / "memory.stat", "anon " + std::to_string(23 * G) + "\nactive_file " + std::to_string(10 * G) + "\ninactive_file " + std::to_string(30 * G) + "\n");
    put(j / "step0/memory.max", "max\n");
    put(j / "step0/memory.current", std::to_string(63 * G) + "\n");
    CHECK(cg::headroom(r.string()) == (size_t)(41 * G), "v2 headroom %zu", cg::headroom(r.string()));
    put(r / "sys/fs/cgroup/slurm/memory.max", std::to_string(30 * G) + "\n");   // an ancestor's tighter limit
    put(r / "sys/fs/cgroup/slurm/memory.current", std::to_string(10 * G) + "\n");
    CHECK(cg::headroom(r.string()) == (size_t)(20 * G), "ancestor headroom %zu", cg::headroom(r.string()));
  }
  { // v1 with memory mounted together with cpu, under a container mount root
    const fs::path r = base / "v1";
    put(r / "proc/self/mountinfo", "31 23 0:27 /docker/abc /sys/fs/cgroup/cpu,memory rw - cgroup cgroup rw,cpu,memory\n");
    put(r / "proc/self/cgroup", "4:cpu,memory:/docker/abc\n3:pids:/docker/abc\n");
    const fs::path m = r / "sys/fs/cgroup/cpu,memory";
    put(m / "memory.limit_in_bytes", std::to_string(8 * G) + "\n");
    put(m / "memory.usage_in_bytes", std::to_string(7 * G) + "\n");
    put(m / "memory.stat", "total_inactive_file " + std::to_string(2 * G) + "\ntotal_active_file " + std::to_string(1 * G) + "\n");
    CHECK(cg::headroom(r.string()) == (size_t)(4 * G), "v1 combined headroom %zu", cg::headroom(r.string()));
  }
  { // hybrid: v2 unified without the memory controller, memory on v1; usage below the cache counts as 0
    const fs::path r = base / "hybrid";
    put(r / "proc/self/mountinfo",
        "30 23 0:26 / /sys/fs/cgroup/unified rw - cgroup2 cgroup2 rw\n"
        "31 23 0:27 / /sys/fs/cgroup/memory rw - cgroup cgroup rw,memory\n");
    put(r / "proc/self/cgroup", "0::/user.slice\n5:memory:/user.slice\n");
    put(r / "sys/fs/cgroup/unified/user.slice/cgroup.procs", "1\n");
    const fs::path m = r / "sys/fs/cgroup/memory/user.slice";
    put(m / "memory.limit_in_bytes", std::to_string(16 * G) + "\n");
    put(m / "memory.usage_in_bytes", std::to_string(1 * G) + "\n");
    put(m / "memory.stat", "total_inactive_file " + std::to_string(2 * G) + "\n");
    CHECK(cg::headroom(r.string()) == (size_t)(16 * G), "hybrid headroom %zu", cg::headroom(r.string()));
  }
  { // no mountinfo: the standard v2 mount; unlimited v1 (a huge limit) and "max" are no limit
    const fs::path r = base / "fallback";
    put(r / "proc/self/cgroup", "0::/a\n");
    put(r / "sys/fs/cgroup/a/memory.max", "max\n");
    put(r / "sys/fs/cgroup/a/memory.current", "5\n");
    CHECK(cg::headroom(r.string()) == SIZE_MAX, "unlimited");
    put(r / "sys/fs/cgroup/a/memory.max", "1000\n");
    CHECK(cg::headroom(r.string()) == 995, "fallback mount %zu", cg::headroom(r.string()));
    CHECK(cg::headroom((base / "nothing").string()) == SIZE_MAX, "no files");
  }
  fs::remove_all(base);
  if (g_fail) { std::fprintf(stderr, "%d check(s) failed\n", g_fail); return 1; }
  std::printf("cgroup_mem: all checks passed\n");
  return 0;
}
