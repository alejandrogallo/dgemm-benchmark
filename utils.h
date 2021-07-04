#pragma once
#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <map>
#include <string>
#define Q(...) #__VA_ARGS__
#define QUOTE(...) Q(__VA_ARGS__)
#define SHOW_MACRO(x) #x ": " QUOTE(x)
#define SHOW_VAR(x) #x ": " << x

// simple logger
#define LOG if (rank == 0) std::cout

struct Timer {
  using Clock = std::chrono::high_resolution_clock;
  using Event = std::chrono::time_point<Clock>;
  std::chrono::duration<double> duration;
  Event _start;
  Event start() { return _start = Clock::now(); }
  Event stop() {
    Event const _end = Clock::now();
    duration += _end - _start;
    return _end;
  }
  double count() const { return duration.count(); }
};

using Timings = std::map<std::string, Timer>;
