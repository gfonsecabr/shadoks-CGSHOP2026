#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <queue>
#include <cassert>
#include <array>
#include <random>
#include <stdexcept>
#include <chrono>
#include <unistd.h>
#include <sstream>
#include <cmath>
#include <boost/unordered_set.hpp> // hash_combine
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <deque>
#include <map>

int verbose = 2;


const auto beginTime = std::chrono::high_resolution_clock::now();

inline double elapsed() {
  auto end = std::chrono::high_resolution_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end - beginTime);

  return dur.count() / 1000.0;
}

using i64 = long long int;

inline std::string timeString(std::chrono::system_clock::time_point tp = std::chrono::system_clock::now()) {
  std::time_t tt = std::chrono::system_clock::to_time_t(tp);
  struct std::tm * ptm = std::localtime(&tt);
  std::ostringstream oss;
  oss << std::put_time(ptm, "%Y%m%d-%H%M%S");
  std::string s = oss.str();
  return s;
}

static std::random_device dev;
static std::mt19937 rgen(dev()); // Replace dev() with an int for deterministic behavior
// static std::mt19937 rgen(1);

inline std::string exec(std::string cmd) {
  // std::cout << std::endl << cmd << std::endl;
  std::array<char, 128> buffer;
  std::string result;
  auto voidpclose = [](FILE *f){pclose(f);};
  std::unique_ptr<FILE, decltype(voidpclose)> pipe(popen(cmd.c_str(), "r"), voidpclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
    // std::cout << result << std::endl;
  }
  return result;
}

inline std::string concatenate(int argc, char **argv) {
  std::string s;
  for(int i = 0; i < argc; i++) {
    s += argv[i];
    if(i != argc-1)
      s += " ";
  }
  return s;
}

inline bool replace(std::string& str, const std::string& from, const std::string& to) {
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}


inline bool fileExists(const std::string &name) {
  std::ifstream f(name);
  return f.good();
}
