#pragma once
#include "tools.hpp"
#include "triangulation.hpp"

struct VecHash {
  size_t operator()(const std::vector<int>& v) const noexcept {
    // Standard hashing technique: combine seed with each element
    // using a 64-bit rotate and a constant (similar to Boost)
    size_t h = v.size();
    for (int x : v) {
      // Mix: shift-left 1, shift-right 7, XOR with x hashed.
      h ^= std::hash<int>()(x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    }
    return h;
  }
};

class CrossLB {
  using IEdge = std::array<short,2>;
  inline static std::unordered_map<std::vector<int>, int, VecHash> presolved;

  static void cutborders(std::vector<int>& oc) {
    if (oc.size() == 1) {
      if (oc[0] == 1 || oc[0] == 2) {
        oc.erase(oc.begin());
      } else {
        oc[0] = oc[0] - 2;
      }
    } else {
      if (oc.size() == 0) return;
      if (oc[0] == 1) {
        oc.erase(oc.begin());
      } else {
        oc[0] = oc[0] - 1;
      }
      if (oc.empty()) return;
      if (oc.back() == 1) {
        oc.pop_back();
      } else {
        oc.back() = oc.back() - 1;
      }
    }
  }

  static void reducepairs(std::vector<int>& oc, const std::vector<int>& intouchables) {
    if (oc.size() >= 1) {
      size_t i = 0;
      size_t n = oc.size();
      while (i < n) {
        if (oc[i] > 1) {
          double diff = static_cast<double>(oc[i] - intouchables[i]) / 2.0;
          oc[i] = static_cast<int>(ceil(diff)) + intouchables[i];
        }
        ++i;
      }
    }
  }

  static void switch_op(std::vector<int>& oc, std::vector<std::pair<int,int>>& ones, std::vector<int>& intouchables, int start) {
    int n = (int)oc.size();

    // if oc[start]==1 and oc[start+1]==1:
    if (oc[start] == 1 && oc[start+1] == 1) {

      if (start == 0) {
        // oc.pop(0)
        oc.erase(oc.begin());

        // oc[1] = oc[1] + 1
        oc[1] = oc[1] + 1;

        // intouchables[0] = 1
        intouchables[0] = 1;

        // intouchables[1] = intouchables[1] + 1
        intouchables[1] = intouchables[1] + 1;

        // for j in range(len(ones)): ones[j][0]-- ; ones[j][1]--
        for (int j = 0; j < (int)ones.size(); ++j) {
          ones[j].first  = ones[j].first  - 1;
          ones[j].second = ones[j].second - 1;
        }
      }

      else if (start == n-2) {
        // oc.pop(start+1)
        oc.erase(oc.begin() + (start + 1));

        // oc[start-1] += 1
        oc[start-1] = oc[start-1] + 1;

        // intouchables[start-1] += 1
        intouchables[start-1] = intouchables[start-1] + 1;

        // intouchables[start] = 1
        intouchables[start] = 1;
      }

      else {
        // oc[start-1] += 1
        oc[start-1] = oc[start-1] + 1;

        // oc[start+2] += 1
        oc[start+2] = oc[start+2] + 1;

        // intouchables[start-1] += 1
        intouchables[start-1] = intouchables[start-1] + 1;

        // intouchables[start+1] += 1
        intouchables[start+1] = intouchables[start+1] + 1;

        // oc.pop(start)
        oc.erase(oc.begin() + start);

        // oc.pop(start)
        oc.erase(oc.begin() + start);

        // intouchables.pop(start)
        intouchables.erase(intouchables.begin() + start);

        // intouchables.pop(start)
        intouchables.erase(intouchables.begin() + start);

        // for j in range(len(ones)): ones[j][0] -= 2 ; ones[j][1] -= 2
        for (int j = 0; j < (int)ones.size(); ++j) {
          ones[j].first  = ones[j].first  - 2;
          ones[j].second = ones[j].second - 2;
        }
      }
    }

    // else:  # oc[start+1] > 1
    else {

      if (start == 0 && oc[0] == 1) {
        // used in 121 combined flips
        if (intouchables[0] == 0 && intouchables[1] < 2) {

          // oc.insert(0,1)
          oc.insert(oc.begin(), 1);

          // oc[2] = oc[2] - 1
          oc[2] = oc[2] - 1;

          // intouchables.insert(0,1)
          intouchables.insert(intouchables.begin(), 1);

          // intouchables[1] = 1
          intouchables[1] = 1;

          // update ones (all +1,+1)
          for (int j = 0; j < (int)ones.size(); ++j) {
            ones[j].first  = ones[j].first  + 1;
            ones[j].second = ones[j].second + 1;
          }
        }
      }

      else if (start == 0) {
        if (oc[1] == 1) {
          // oc[0] -= 1
          oc[0] = oc[0] - 1;

          // intouchables[1] = 1
          intouchables[1] = 1;

          if ((int)oc.size() > 2) {
            // oc[2] += 1
            oc[2] = oc[2] + 1;

            // intouchables[2] = 1
            intouchables[2] = 1;
          }
          else if ((int)oc.size() == 1) {   // likely useless
            oc.push_back(1);
            intouchables.push_back(1);
          }
        }
      }

      else {
        // start > 0

        if (start == n-2 && oc[n-1] == 1) {
          // used in 121 combined flips
          if (intouchables[n-2] > 0 && intouchables[n-1] > 1) {

            // oc.append(1)
            oc.push_back(1);

            // intouchables[n-1] = 1
            intouchables[n-1] = 1;

            // intouchables[n] = 1
            intouchables.push_back(1);

            // oc[n-2] -= 1
            oc[n-2] = oc[n-2] - 1;
          }
        }

        else if (start == n-2) {
          if (oc[n-2] == 1) {
            // oc[n-1] -= 1
            oc[n-1] = oc[n-1] - 1;

            // oc[n-3] += 1
            oc[n-3] = oc[n-3] + 1;

            // intouchables[n-3] += 1
            intouchables[n-3] = intouchables[n-3] + 1;

            // intouchables[n-2] = 1
            intouchables[n-2] = 1;
          }
        }

        else {
          if (oc[start] != 1 && oc[start+1] != 1) {
            // nothing to do
          }
          else if (oc[start] == 1) {
            if (intouchables[start] == 0) {

              // oc[start-1] += 1
              oc[start-1] = oc[start-1] + 1;

              // oc[start+1] -= 1
              oc[start+1] = oc[start+1] - 1;

              // intouchables[start] = 1
              intouchables[start] = 1;

              // intouchables[start-1] += 1
              intouchables[start-1] = intouchables[start-1] + 1;
            }
          }

          else if (oc[start+1] == 1) {
            if (intouchables[start+1] == 0) {

              // oc[start] += 1
              oc[start] = oc[start] + 1;

              // oc[start+2] -= 1
              oc[start+2] = oc[start+2] - 1;

              // intouchables[start+1] = 1
              intouchables[start+1] = 1;

              // intouchables[start] += 1
              intouchables[start] = intouchables[start] + 1;
            }
          }
        }
      }
    }
  }


  static std::vector<int> permute3(std::vector<int>& oc, int verbose = 0) {
    std::vector<std::pair<int,int>> ones;
    int start = -1;
    int end = -1;
    int n = static_cast<int>(oc.size());
    for (int i = 0; i < n; ++i) {
      int g = oc[i];
      if (g == 1) {
        if (start == -1) {
          start = i;
          end = -1;
        }
      } else {
        if (start != -1) {
          end = i - 1;
          ones.emplace_back(start, end);
          start = -1;
        }
      }
    }
    if (end == -1 && start != -1) {
      end = static_cast<int>(oc.size()) - 1;
      ones.emplace_back(start, end);
    }

    std::vector<int> intouchables(oc.size(), 0);
    if (ones.size() > 0) {
      if (oc.size() > 2) {
        // i=0  (unused in Python except as loop variable)
        for (auto pr : ones) {
          start = pr.first;
          end = pr.second;
          if (end - start == 0) { // isolated one
            if (start == 0 || start == (int)oc.size() - 1) {
              // a = 0 (do nothing)
            } else if (start == 1) {
              double u = 0.0;
              if (oc[0] - 1 > 0) u = log2(static_cast<double>(oc[0] - 1));
              if (oc.size() > 3) {
                if (floor(u) == u && u > 0) {
                  switch_op(oc, ones, intouchables, 0);
                }
              } else if (oc.size() == 3) { // 8 case X1Y
                int parity0 = oc[0] % 2;
                int parity2 = oc[2] % 2;
                if (parity0 == 1 && parity2 == 1) {
                  int uu0 = oc[0] / 2;
                  int uu1 = oc[2] / 2;
                  if (uu0 % 2 == 1 && uu1 % 2 == 1) {
                    switch_op(oc, ones, intouchables, 0);
                  }
                } else if (parity0 == 1 && parity2 == 0) {
                  int uu0 = oc[0] / 2;
                  int uu1 = oc[2] / 2;
                  if (uu0 % 2 == 1 && uu1 % 2 == 0) {
                    switch_op(oc, ones, intouchables, 0);
                  }
                } else if (parity0 == 0 && parity2 == 1) {
                  int uu0 = oc[0] / 2;
                  int uu1 = oc[2] / 2;
                  if (uu0 % 2 == 0 && uu1 % 2 == 1) {
                    switch_op(oc, ones, intouchables, 1);
                  }
                }
              }
            } else if (start == (int)oc.size() - 2) {
              double u = 0.0;
              if (oc.back() - 1 > 0) u = log2(static_cast<double>(oc.back() - 1));
              if (floor(u) == u && u > 0) {
                switch_op(oc, ones, intouchables, start);
              }
            }
          } else if (end - start == 1) {
            switch_op(oc, ones, intouchables, start);
          } else if (end - start == 2) {
            if (start == 0) {
              switch_op(oc, ones, intouchables, start);
            } else if (end == (int)oc.size() - 1) {
              switch_op(oc, ones, intouchables, start + 1);
            } else if (start == 1 && end == (int)oc.size() - 2) {
            //   int u = static_cast<int>(ceil(log2(static_cast<double>(oc[0]))));
            //   int v = static_cast<int>(ceil(log2(static_cast<double>(oc[0] + 1))));
            //   int uu = static_cast<int>(ceil(log2(static_cast<double>(oc[end+1]))));
            //   int vv = static_cast<int>(ceil(log2(static_cast<double>(oc[end+1] + 1))));
              if (oc[0] >= oc[end+1]) {
                switch_op(oc, ones, intouchables, 1);
              } else {
                switch_op(oc, ones, intouchables, end-1);
              }
            } else if (start == 1) {
              int u = static_cast<int>(ceil(log2(static_cast<double>(oc[0]))));
              int v = static_cast<int>(ceil(log2(static_cast<double>(oc[0] + 1))));
              if (u != v && u > 0) {
                switch_op(oc, ones, intouchables, 2);
              } else {
                switch_op(oc, ones, intouchables, 1);
              }
            } else if (end == (int)oc.size() - 2) {
              int u = static_cast<int>(ceil(log2(static_cast<double>(oc[end+1]))));
              int v = static_cast<int>(ceil(log2(static_cast<double>(oc[end+1] + 1))));
              if (u != v && u > 0) {
                switch_op(oc, ones, intouchables, start);
              } else {
                switch_op(oc, ones, intouchables, start+1);
              }
            } else {
              if (oc[start-1] >= oc[end+1]) {
                switch_op(oc, ones, intouchables, start+0);
              } else {
                switch_op(oc, ones, intouchables, start+1);
              }
            }
          } else if (end - start == 3) {
            if (start == 0) {
              switch_op(oc, ones, intouchables, start+2);
            } else if (end == (int)oc.size() - 1) {
              switch_op(oc, ones, intouchables, start);
            } else {
              switch_op(oc, ones, intouchables, start+1);
            }
          } else if (end - start == 4) {
            if (start == 0) {
              switch_op(oc, ones, intouchables, start+2);
            } else if (end == (int)oc.size() - 1) {
              switch_op(oc, ones, intouchables, start+1);
            } else {
              switch_op(oc, ones, intouchables, start);
              start = start - 2;
              switch_op(oc, ones, intouchables, start+3);
            }
          } else if ((end - start - 1) % 4 == 0) {
            for (int k = 0; k < 1 + (end - start - 1) / 4; ++k) {
              switch_op(oc, ones, intouchables, start + 2*k);
            }
          } else if ((end - start - 1) % 4 == 2) {
            if (start == 0) {
              for (int k = 0; k < (end - start + 1) / 4; ++k) {
                switch_op(oc, ones, intouchables, start + 2 + 2*k);
              }
            } else if (end == (int)oc.size() - 1) {
              for (int k = 0; k < (end - start + 1) / 4; ++k) {
                switch_op(oc, ones, intouchables, start + 2*k);
              }
            } else {
              for (int k = 0; k < 1 + (end - start - 1) / 4; ++k) {
                switch_op(oc, ones, intouchables, start + 1 + 2*k);
              }
            }
          } else if ((end - start - 1) % 4 == 1) {
            if (start == 0) {
              for (int k = 0; k < 1 + (end - start + 1) / 4; ++k) {
                switch_op(oc, ones, intouchables, start + 2*k);
              }
            } else if (end == (int)oc.size() - 1) {
              for (int k = 0; k < (end - start + 1) / 4; ++k) {
                switch_op(oc, ones, intouchables, start + 2*k);
              }
            } else {
              if (oc[start-1] >= oc[end+1]) {
                for (int k = 0; k < (end - start + 1) / 4; ++k) {
                  switch_op(oc, ones, intouchables, start + 2*k);
                }
              } else {
                for (int k = 0; k < (end - start + 1) / 4; ++k) {
                  switch_op(oc, ones, intouchables, start + 2*k + 1);
                }
              }
            }
          } else if ((end - start - 1) % 4 == 3) {
            if (start == 0) {
              for (int k = 0; k < (end - start + 1) / 4; ++k) {
                switch_op(oc, ones, intouchables, start + 1 + 2*k);
              }
            } else if (end == (int)oc.size() - 1) {
              for (int k = 0; k < (end - start + 1) / 4; ++k) {
                switch_op(oc, ones, intouchables, start + 2*k);
              }
            } else {
              if (oc[start-1] >= oc[end+1]) {
                for (int k = 0; k < (end - start + 1) / 4; ++k) {
                  switch_op(oc, ones, intouchables, start + 2*k);
                }
              } else {
                for (int k = 0; k < (end - start + 1) / 4; ++k) {
                  switch_op(oc, ones, intouchables, start + 2*k + 1);
                }
              }
            }
          }
        } // end for ones

        // onetwoone handling
        std::vector<int> onetwoone;
        if (oc.size() > 4) {
          int firstonetwoone = 0;
          int s = 0;
          while (s < (int)oc.size() - 2) {
            if (oc[s] == 1 && oc[s+1] == 2 && oc[s+2] == 1) {
              if (intouchables[s] == 0 && intouchables[s+2] == 0 && intouchables[s+1] == 0) {
                if (firstonetwoone == 0) firstonetwoone = 1;
                switch_op(oc, ones, intouchables, s);
                if (s == 0) {
                  switch_op(oc, ones, intouchables, 2);
                } else {
                  switch_op(oc, ones, intouchables, s+1);
                }
              }
            }
            ++s;
          }
        }

      } // end if oc.size()>2
    } // end if ones.size()>0

    return intouchables;
  }

public:
  static int lb(std::vector<int> oc1) {
    int countflips = 1;
    if (oc1.empty()) {
      return 1;
    }

    while (!oc1.empty()) {
      auto it = presolved.find(oc1);
      if (it != presolved.end()) {
        return it->second + countflips - 1;
      }
      cutborders(oc1);
      it = presolved.find(oc1);
      if (it != presolved.end()) {
        return it->second + countflips - 1;
      }
      std::vector<int> intouchables = permute3(oc1);
      reducepairs(oc1, intouchables);
      countflips = countflips + 1;
    }
    return countflips;
  }

  static void load_lbfile(const std::string& filename) {
    int lines = 0;
    if(!presolved.empty())
      return;
    std::ifstream file(filename);
    if (!file.is_open()) {
      // If file doesn't exist we silently continue (like Python would throw)
      // But we keep behavior robust.
      std::cerr << "Warning: could not open " << filename << " (presolved remains empty)\n";
      return;
    }

    std::string line;
    while (getline(file, line)) {
      // trim
      auto l = line;
      // remove leading/trailing whitespace
      auto is_space = [](char c){ return c==' ' || c=='\t' || c=='\r' || c=='\n'; };
      size_t start = 0;
      while (start < l.size() && is_space(l[start])) ++start;
      size_t end = l.size();
      while (end > start && is_space(l[end-1])) --end;
      if (start >= end) continue;
      l = l.substr(start, end-start);

      // split by " = "
      std::string delim = " = ";
      size_t pos = l.find(delim);
      if (pos == std::string::npos) continue;
      std::string left = l.substr(0, pos);
      std::string right = l.substr(pos + delim.size());

      // left is like "[1, 2, 3]" in Python: strip [] and spaces
      if (!left.empty() && left.front() == '[' && left.back() == ']') {
        left = left.substr(1, left.size()-2);
      }
      // if left is empty => empty tuple - handle
      std::vector<int> numbers;
      if (!left.empty()) {
        std::stringstream ss(left);
        std::string token;
        while (getline(ss, token, ',')) {
          // trim token
          size_t a = 0; while (a < token.size() && is_space(token[a])) ++a;
          size_t b = token.size(); while (b > a && is_space(token[b-1])) --b;
          if (b > a) {
            std::string t = token.substr(a, b-a);
            try {
              int x = stoi(t);
              numbers.push_back(x);
            } catch (...) {
              // ignore parse errors
            }
          }
        }
      }

      // parse right as int
      try {
        int val = stoi(right);
        if(!presolved.contains(numbers) || presolved.at(numbers) > val)
          presolved[numbers] = val;
        std::vector<int> numbers2(numbers.rbegin(),numbers.rend());
        if(!presolved.contains(numbers2) || presolved.at(numbers2) > val)
          presolved[numbers2] = val;
        lines++;
      } catch (...) {
        // ignore parse errors
      }
    }
    std::cout << "Read crossing lower bound file with " << lines << " lines, getting " << presolved.size() << " entries\n";
  }

  static int simpleCrossLB(const std::vector<IEdge> &crossings) {
    if(crossings.size() <= 2)
      return crossings.size();

    int lb = 1+std::floor(std::log2(crossings.size()));
    return lb;
  }

  static int fancyCrossLB(const std::vector<IEdge> &crossings) {
    if(crossings.size() <= 2)
      return crossings.size();

    std::vector<int> sides = buildSides(crossings);
    std::vector<int> sidesRev(sides.rbegin(),sides.rend());
    return std::min(CrossLB::lb(sides), CrossLB::lb(sidesRev));
  }

  static short commonElement(IEdge a, IEdge b) {
    std::vector<short> v_common;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(v_common));
    assert(v_common.size() == 1);
    return v_common.at(0);
  }

  static std::vector<int> buildSides(const std::vector<IEdge> &crossings) {
    assert(crossings.size() >= 2);

    short common = commonElement(crossings[0], crossings[1]);

    std::vector<int> ret;
    int x = 1;
    for(size_t i = 1; i + 1 < crossings.size(); i++) {
      short newCommon = commonElement(crossings[i], crossings[i+1]);
      if(common != newCommon) {
        common = newCommon;
        ret.push_back(x);
        x = 0;
      }
      x++;
    }
    ret.push_back(x);
    // for(int i : ret)
    //   std:: cout << i << " + ";
    // std::cout << " (1) = " << crossings.size() << std::endl;

    return ret;
  }
};
