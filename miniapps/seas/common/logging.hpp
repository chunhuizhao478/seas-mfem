// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_LOGGING_HPP
#define MFEM_SEAS_LOGGING_HPP

#include <iostream>
#include <string>

namespace mfem
{
namespace seas
{

/// Verbosity levels for SEAS logging.
enum class LogLevel { Silent = 0, Error = 1, Warning = 2, Info = 3, Debug = 4 };

/// @brief Simple logging wrapper for rank-0 output.
///
/// Centralizes the `if (rank == 0) cout << ...` pattern.
/// Only the root rank produces output; other ranks are silent.
class Logger
{
public:
   Logger(int rank = 0, LogLevel level = LogLevel::Info)
      : rank_(rank), level_(level) {}

   void SetLevel(LogLevel level) { level_ = level; }
   LogLevel GetLevel() const { return level_; }
   bool IsRoot() const { return rank_ == 0; }

   /// Log at Info level (default simulation output)
   template <typename... Args>
   void Info(Args &&...args) const
   {
      static_assert(sizeof...(Args) > 0, "Logger::Info requires at least one argument");
      if (rank_ == 0 && level_ >= LogLevel::Info)
      {
         (std::cout << ... << std::forward<Args>(args));
         std::cout << "\n";
      }
   }

   /// Log at Warning level
   template <typename... Args>
   void Warn(Args &&...args) const
   {
      static_assert(sizeof...(Args) > 0, "Logger::Warn requires at least one argument");
      if (rank_ == 0 && level_ >= LogLevel::Warning)
      {
         std::cout << "[WARNING] ";
         (std::cout << ... << std::forward<Args>(args));
         std::cout << "\n";
      }
   }

   /// Log at Error level
   template <typename... Args>
   void Error(Args &&...args) const
   {
      static_assert(sizeof...(Args) > 0, "Logger::Error requires at least one argument");
      if (rank_ == 0 && level_ >= LogLevel::Error)
      {
         std::cerr << "[ERROR] ";
         (std::cerr << ... << std::forward<Args>(args));
         std::cerr << "\n";
      }
   }

   /// Log at Debug level
   template <typename... Args>
   void Debug(Args &&...args) const
   {
      static_assert(sizeof...(Args) > 0, "Logger::Debug requires at least one argument");
      if (rank_ == 0 && level_ >= LogLevel::Debug)
      {
         std::cout << "[DEBUG] ";
         (std::cout << ... << std::forward<Args>(args));
         std::cout << "\n";
      }
   }

private:
   int rank_;
   LogLevel level_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LOGGING_HPP
