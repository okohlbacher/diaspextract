// Copyright (c) 2026, Oliver Kohlbacher and contributors
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

// PortableMath.h: within 1 ulp of the platform libm, the special cases, and one hash over 10^6 fixed inputs that must be
// the SAME on every platform (that is the header's whole point; the platform libm's own hash differs on each).
#include "../src/PortableMath.h"
#if defined(__GNUC__) && ! defined(__clang__)
  #pragma GCC optimize("fp-contract=off") // the inputs must not be fused either (clang: the pragma in main)
#endif
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>

static int fails = 0;
#define CHECK(c, ...)                    \
  do                                     \
  {                                      \
    if (! (c))                           \
    {                                    \
      std::printf("FAIL  " __VA_ARGS__); \
      std::printf("\n");                 \
      ++fails;                           \
    }                                    \
  } while (0)

static int64_t ulps(double a, double b)
{
  if (a == b) { return 0; }
  int64_t ia, ib;
  std::memcpy(&ia, &a, 8);
  std::memcpy(&ib, &b, 8);
  if (ia < 0) { ia = INT64_MIN - ia; }
  if (ib < 0) { ib = INT64_MIN - ib; }
  return ia > ib ? ia - ib : ib - ia;
}

int main()
{
#if defined(__clang__)
  #pragma clang fp contract(off) // the inputs below must not be fused either: the header guards only its own arithmetic
#endif
  const double inf = INFINITY;
  CHECK(pmath::exp(0.0) == 1.0, "exp(0)");
  CHECK(ulps(pmath::exp(1.0), 2.718281828459045) <= 1, "exp(1)");
  CHECK(pmath::exp(-inf) == 0.0 && pmath::exp(inf) == inf && std::isnan(pmath::exp(NAN)), "exp specials");
  CHECK(pmath::exp(710.0) == inf && pmath::exp(-746.0) == 0.0 && pmath::exp(-740.0) > 0.0, "exp range ends");
  CHECK(pmath::log(1.0) == 0.0 && pmath::log(0.0) == -inf && std::isnan(pmath::log(-1.0)) && pmath::log(inf) == inf, "log specials");
  CHECK(std::fabs(pmath::log(5e-324) - (-744.4400719213812)) < 1e-12, "log of the smallest subnormal");
  CHECK(pmath::log1p(0.0) == 0.0 && pmath::log1p(1e-300) == 1e-300 && pmath::log1p(-1.0) == -inf, "log1p specials");
  CHECK(pmath::tanh(0.0) == 0.0 && pmath::tanh(30.0) == 1.0 && pmath::tanh(-30.0) == -1.0 && pmath::atanh(0.0) == 0.0, "tanh/atanh specials");
  CHECK(pmath::expm1(0.0) == 0.0 && pmath::expm1(1e-300) == 1e-300 && pmath::expm1(-800.0) == -1.0, "expm1 specials");
  CHECK(pmath::pow(0.0, 2.0) == 0.0 && pmath::pow(4.0, 0.5) == 2.0 && pmath::pow(2.0, 10.0) == 1024.0, "pow");

  int64_t worst[6] = {0, 0, 0, 0, 0, 0};
  uint64_t s = 88172645463325252ULL, h = 1469598103934665603ULL;
  auto eat = [&](double v) {
    uint64_t b;
    std::memcpy(&b, &v, 8);
    h = (h ^ b) * 1099511628211ULL;
  };
  for (int i = 0; i < 1000000; ++i)
  {
    s ^= s << 13;
    s ^= s >> 7;
    s ^= s << 17;
    const double u = (double)(s >> 11) / 9007199254740992.0; // [0, 1)
    const double xe = 1400.0 * u - 720.0, xl = std::ldexp(0.5 + u, (int)(s % 2000) - 1000), x1 = 1e-12 + 1e4 * u * u * u;
    const double e = pmath::exp(xe), l = pmath::log(xl), l1 = pmath::log1p(x1), p = pmath::pow(0.01 + 9.0 * u, 2.5);
    const double xr = 1.998 * u - 0.999, th = pmath::tanh(8.0 * xr), at = pmath::atanh(xr);
    eat(e);
    eat(l);
    eat(l1);
    eat(p);
    eat(th);
    eat(at);
    {
      const int64_t d = ulps(th, std::tanh(8.0 * xr));
      if (d > worst[4]) { worst[4] = d; }
    }
    {
      const int64_t d = ulps(at, std::atanh(xr));
      if (d > worst[5]) { worst[5] = d; }
    }
    if (std::isnormal(e))
    {
      const int64_t d = ulps(e, std::exp(xe));
      if (d > worst[0]) { worst[0] = d; }
    }
    {
      const int64_t d = ulps(l, std::log(xl));
      if (d > worst[1]) { worst[1] = d; }
    }
    {
      const int64_t d = ulps(l1, std::log1p(x1));
      if (d > worst[2]) { worst[2] = d; }
    }
    {
      const int64_t d = ulps(p, std::pow(0.01 + 9.0 * u, 2.5));
      if (d > worst[3]) { worst[3] = d; }
    }
  }
  std::printf("max ulps vs this libm: exp %lld log %lld log1p %lld pow %lld tanh %lld atanh %lld; hash %016llx\n", (long long)worst[0],
              (long long)worst[1], (long long)worst[2], (long long)worst[3], (long long)worst[4], (long long)worst[5], (unsigned long long)h);
  CHECK(worst[0] <= 1 && worst[1] <= 1, "exp/log within 1 ulp of the libm");
  CHECK(worst[2] <= 4 && worst[3] <= 64 && worst[4] <= 8 && worst[5] <= 8, "log1p within 4, tanh/atanh 8, pow 64 ulps");
  // the cross-platform pin: macOS arm64 (Apple clang), Linux x86-64 (GCC), Windows x64 (clang-cl) must all print it
  const uint64_t pin = 0x3cba1a03b8ba67e3ULL;
  CHECK(h == pin, "hash %016llx is not the pinned %016llx: this platform computes different bits", (unsigned long long)h, (unsigned long long)pin);
  std::printf(fails ? "[portable_math] %d FAILED\n" : "[portable_math] selftest passed\n", fails);
  return fails ? 1 : 0;
}
