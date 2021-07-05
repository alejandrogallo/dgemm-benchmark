#pragma once

#include <string>
#include <algorithm>
#include <vector>
#include <iostream>

#define INSTANTIATE_VALUE(__type, __reader)           \
  template <> __type value(Args &as, Flag &f) {       \
    auto const v(value<std::string>(as, f));          \
    return __reader(v.c_str());                       \
  }

namespace hauta {
  typedef const std::string Flag;
  typedef std::string Arg;
  typedef const std::vector<Arg> Args;

  bool isFlagPresent(Args &a, Flag &f) {
    return std::find(a.begin(), a.end(), f) != a.end();
  }

  // value
  template<typename F> F value(Args &args, Flag &f);

  template<> std::string value(Args &a, Flag &f) {
    const auto it(std::find(a.begin(), a.end(), f));
    if (!isFlagPresent(a, f)) {
      std::cerr << "Expecting flag " << f << "\n";
      throw "";
    }
    return std::string(*(it+1));
  }

  INSTANTIATE_VALUE(size_t, std::atoi)
  INSTANTIATE_VALUE(int,    std::atoi)
  INSTANTIATE_VALUE(double, std::atof)
  INSTANTIATE_VALUE(float,  std::atof)

  template<typename F> F value(Args &args, Flag &f, F const def) {
    return isFlagPresent(args, f)
         ? value<F>(args, f)
         : def
         ;
  }

  template <typename F>
  F value(int argc, char **argv, Flag& f, const F def) {
    return value<F>(Args {argv, argv + argc}, f, def);
  }
  template <typename F>
  F value(int argc, char **argv, Flag& f) {
    return value<F>(Args {argv, argv + argc}, f);
  }

}

#undef INSTANTIATE_VALUE
