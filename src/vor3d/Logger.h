#pragma once

#include <iostream>
#include <ctime>

// ANSI color codes
#define COLOR_RESET   "\033[0m"
#define COLOR_GREEN   "\033[32m"
#define COLOR_RED     "\033[31m"
#define COLOR_YELLOW  "\033[33m"
#define COLOR_BLUE    "\033[34m"

#define LOG_PASS(msg) std::cout << COLOR_GREEN << "\u2714 Pass - " << msg << COLOR_RESET << std::endl
#define LOG_FAIL(msg) std::cerr << COLOR_RED   << "\u2718 Fail - " << msg << COLOR_RESET << std::endl
#define LOG_INFO(msg) std::cout << COLOR_BLUE  << "\u2139 Info - " << msg << COLOR_RESET << std::endl
#define LOG_WARN(msg) std::cout << COLOR_YELLOW<< "\u26A0 Warn - " << msg << COLOR_RESET << std::endl

#define LOG_TIME() \
    do { \
        std::time_t now = std::time(nullptr); \
        std::string timeStr = std::ctime(&now); \
        timeStr.pop_back(); /* remove trailing newline */ \
        std::cout << COLOR_BLUE << "\u2022 Time - " << timeStr << COLOR_RESET << std::endl; \
    } while(0)