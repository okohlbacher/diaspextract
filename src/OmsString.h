// OpenMS 3.7 replaced OpenMS::String with std::string. The tool builds its log lines and the mzML header stamps with
// String(x); this keeps those conversions, formatted by StringUtils::toStr, which carries over String's number
// formatting -- so the stamps in the output do not change.
#ifndef DIASPEXTRACT_OMS_STRING_H
#define DIASPEXTRACT_OMS_STRING_H

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <string>
#include <type_traits>
#include <utility>

namespace diaspextract
{
struct String : std::string
{
  String() = default;
  String(const std::string& s) : std::string(s) {}
  String(std::string&& s) : std::string(std::move(s)) {}
  String(const char* s) : std::string(s) {}
  String(char c) : std::string(1, c) {}
  template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, char> && !std::is_same_v<T, bool>, int> = 0>
  String(T v) : std::string(OpenMS::StringUtils::toStr(v)) {}
  String(bool b) : std::string(OpenMS::StringUtils::toStr(static_cast<int>(b))) {}
  String(float v, bool full_precision = true) : std::string(OpenMS::StringUtils::toStr(v, full_precision)) {}
  String(double v, bool full_precision = true) : std::string(OpenMS::StringUtils::toStr(v, full_precision)) {}
  String(long double v, bool full_precision = true) : std::string(OpenMS::StringUtils::toStr(v, full_precision)) {}
};
} // namespace diaspextract

#endif
