// Copyright (c) 2026, Oliver Kohlbacher and contributors
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

// OpenMS 3.7 replaced OpenMS::String with std::string. The tool builds its log lines and the mzML header stamps with
// String(x); this keeps those conversions, formatted by StringUtils::toStr, which carries over String's number
// formatting -- so the stamps in the output do not change.
#pragma once

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <string>
#include <type_traits>
#include <utility>

/// @brief Namespace of the DIAspeXtract tool.
namespace DIAspeXtract
{
/**
  @brief A std::string with the converting constructors of the former OpenMS::String.

  Numbers and characters are converted with OpenMS::StringUtils::toStr, so log lines and mzML header stamps keep
  the formatting OpenMS::String gave them.
*/
struct String : std::string
{
  /// @brief Constructs an empty string.
  String() = default;
  /**
    @brief Copies a std::string.
    @param[in] s The string to copy.
  */
  String(const std::string& s): std::string(s) // NOLINT(runtime/explicit): implicit by design, as OpenMS::String was
  {
  }
  /**
    @brief Moves from a std::string.
    @param[in,out] s The string to move from; left in a valid but unspecified state.
  */
  String(std::string&& s): std::string(std::move(s)) // NOLINT(runtime/explicit): implicit by design, as OpenMS::String was
  {
  }
  /**
    @brief Copies a null-terminated C string.
    @param[in] s The C string to copy; must not be null.
  */
  String(const char* s): std::string(s) // NOLINT(runtime/explicit): implicit by design, as OpenMS::String was
  {
  }
  /**
    @brief Constructs a string of one character.
    @param[in] c The character.
  */
  String(char c): std::string(1, c) // NOLINT(runtime/explicit): implicit by design, as OpenMS::String was
  {
  }
  /**
    @brief Converts an integer (any integral type except char and bool) to its decimal text.
    @param[in] v The integer value.
  */
  template<class T, std::enable_if_t<std::is_integral_v<T> && ! std::is_same_v<T, char> && ! std::is_same_v<T, bool>, int> = 0>
  String(T v): std::string(OpenMS::StringUtils::toStr(v)) // NOLINT(runtime/explicit): implicit by design, as OpenMS::String was
  {
  }
  /**
    @brief Converts a bool to "1" or "0".
    @param[in] b The value.
  */
  String(bool b): std::string(OpenMS::StringUtils::toStr(static_cast<int>(b))) // NOLINT(runtime/explicit): implicit by design, as OpenMS::String was
  {
  }
  /**
    @brief Converts a float to text with OpenMS::StringUtils::toStr.
    @param[in] v The value.
    @param[in] full_precision Whether all significant digits are written (default true).
  */
  String(float v, bool full_precision = true): // NOLINT(runtime/explicit)
      std::string(OpenMS::StringUtils::toStr(v, full_precision))
  {
  }
  /**
    @brief Converts a double to text with OpenMS::StringUtils::toStr.
    @param[in] v The value.
    @param[in] full_precision Whether all significant digits are written (default true).
  */
  String(double v, bool full_precision = true): // NOLINT(runtime/explicit)
      std::string(OpenMS::StringUtils::toStr(v, full_precision))
  {
  }
  /**
    @brief Converts a long double to text with OpenMS::StringUtils::toStr.
    @param[in] v The value.
    @param[in] full_precision Whether all significant digits are written (default true).
  */
  String(long double v, bool full_precision = true): // NOLINT(runtime/explicit)
      std::string(OpenMS::StringUtils::toStr(v, full_precision))
  {
  }
};
} // namespace DIAspeXtract
