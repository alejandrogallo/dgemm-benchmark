#pragma once
#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <map>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "hauta.h"
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
  void start() { _start = Clock::now(); }
  void stop() { duration += Clock::now() - _start; }
  double count() const { return duration.count(); }
};

using Timings = std::map<std::string, Timer>;

struct Averager {
  std::vector<double> values;
  void push(double flopValue) { values.push_back(flopValue); }
  double count() const {
    return
      std::accumulate(values.begin(), values.end(), 0.0) / values.size();
  }
  size_t size() const { return values.size(); }
  double sigma() const {
    const double mu = count();
    double sigma = 0.0;
    for (auto x: values) {
      sigma += (x - mu) * (x - mu);
    }
    sigma /= values.size();
    return std::sqrt(sigma);
  }
};

using Averages = std::map<std::string, Averager>;
