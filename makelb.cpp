// g++ --static -Wall -Ofast -o makelb -std=c++20 makelb.cpp
#include <vector>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <unistd.h>
#include <chrono>

const auto beginTime = std::chrono::high_resolution_clock::now();

inline double elapsed() {
  auto end = std::chrono::high_resolution_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end - beginTime);

  return dur.count() / 1000.0;
}

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

  std::unordered_map<std::vector<int>, int, VecHash> presolved;

void cutborders(std::vector<int>& oc) {
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

void reducepairs(std::vector<int>& oc, const std::vector<int>& intouchables) {
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

void switch_op(std::vector<int>& oc, std::vector<std::pair<int,int>>& ones, std::vector<int>& intouchables, int start) {
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


std::vector<int> permute3(std::vector<int>& oc, int verbose = 0) {
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

int lbOneWay(std::vector<int> oc1) {
  int countflips = 1;
  if (oc1.empty()) {
    return 1;
  }

  while (!oc1.empty()) {
    cutborders(oc1);
    auto it = presolved.find(oc1);
    if (it != presolved.end()) {
      return it->second + countflips - 1;
    }
    std::vector<int> intouchables = permute3(oc1);
    reducepairs(oc1, intouchables);
    countflips = countflips + 1;
  }
  return countflips;
}

int lb(std::vector<int> oc1) {
  std::vector<int> ocr(oc1.rbegin(),oc1.rend()); // Reverse
  return std::min(lbOneWay(oc1), lbOneWay(ocr));
}

void load_lbfile(const std::string& filename, bool force = false) {
  int lines = 0;
  if(!force && !presolved.empty())
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
  std::cout << "Read lower bound file with " << lines << " lines, getting " << presolved.size() << " entries\n";
}

// -----------------------------------------------------------------------------
// Data structure matching Python: (oc0, newscore+1, octemp, binaire)
struct SpecialEntry {
    std::vector<int> oc0;       // s[0]
    int flips;             // s[1]
    std::vector<int> octemp;    // s[2]
    std::vector<int> binaire;   // s[3] (0/1 std::vector)
};

// binary_search_in_special(special, target)
static int binary_search_in_special(const std::vector<SpecialEntry>& special,
                                    const std::vector<int>& target)
{
    int left = 0, right = (int)special.size() - 1;
    while (left <= right) {
        int mid = (left + right) / 2;
        if (special[mid].oc0 == target) return mid;
        if (special[mid].oc0 < target) left = mid + 1;   // lexicographic compare
        else right = mid - 1;
    }
    // not found: return insertion point encoded like Python code expects
    return -left - 1;
}

// translate_word_to_multiplicities(word)
static std::vector<int> translate_word_to_multiplicities(const std::vector<int>& word) {
    std::vector<int> oc = {0};
    int letter;
    if (!word.empty() && word[0] == 0) letter = 0;
    else letter = 1;

    int i = 0;
    for (int w : word) {
        if (w == letter) {
            oc[i] = oc[i] + 1;
        } else {
            letter = 1 - letter;
            i = i + 1;
            oc.push_back(1);
        }
    }
    return oc;
}

// translate_multiplicities_to_word(oc)
static std::vector<int> translate_multiplicities_to_word(const std::vector<int>& oc) {
    std::vector<int> word;
    int letter = 1;
    int j = 0; // kept line-by-line (even if not used for logic)
    for (int k : oc) {
        for (int i = 0; i < k; ++i) {
            word.push_back(letter);
            j = j + 1;
        }
        letter = 1 - letter;
    }
    return word;
}

// convert_integer(n): run-length encoding of the binary representation of n
static std::vector<int> convert_integer(int n) {
    if (n == 0) return {};

    // binaire = bin(n)[2:]
    std::string binaire;
    {
        std::string tmp;
        while (n > 0) {
            tmp.push_back((n & 1) ? '1' : '0');
            n >>= 1;
        }
        std::reverse(tmp.begin(), tmp.end());
        binaire = tmp;
    }

    std::vector<int> result;
    int i = 0;
    int longueur = (int)binaire.size();
    while (i < longueur) {
        char current_char = binaire[i];
        int c = 1;
        while (i + 1 < longueur && binaire[i + 1] == current_char) {
            i = i + 1;
            c = c + 1;
        }
        result.push_back(c);
        i = i + 1;
    }
    return result;
}

// bfs(b0,b1)
static std::vector<SpecialEntry> bfs(int b0, int b1) {
    std::vector<SpecialEntry> special;
    std::cout << "computes " << b1 << "-" << b0 << " lower bounds\n";

    std::vector<int> used(30, 0);

    for (int m = b0; m < b1; ++m) {
        std::vector<int> oc0 = convert_integer(m);

        // check whether log2(s+2) fits with heuristic
        int score = lbOneWay(oc0);

        int s = 0;
        for (int x : oc0) s += x;
        int h = (int)ceil(log2((double)s + 2.0));

        if (score > h) {
            cutborders(oc0);
            int n = (int)oc0.size();

            if ((int)oc0.size() > 2) { // otherwise heuristic is fine
                // bbb = list(product([0,1], repeat=n-1))
                int bits = n - 1;
                long long total = 1LL << bits;

                for (long long mask = 0; mask < total; ++mask) {
                    int useless = 0;
                    std::vector<int> word = translate_multiplicities_to_word(oc0);

                    for (int i = 0; i < (int)word.size(); ++i) {
                        used[i] = 0; // reinitialize used
                    }

                    int indexinbinaire = 0;

                    // first step: do the switch
                    for (int i = 0; i < (int)word.size() - 1; ++i) {
                        if (word[i] != word[i + 1] && used[i] == 0) {
                            int bit = (int)((mask >> (bits - 1 - indexinbinaire)) & 1LL);

                            if (bit == 1 && used[i] == 0) {
                                // we switch the two letters (if possible)
                                used[i] = 1;
                                used[i + 1] = 1;
                                int letter = word[i + 1];
                                word[i + 1] = word[i];
                                word[i] = letter;
                            } else if (bit == 1 && used[i] == 1) {
                                useless = 1;
                                break;
                            }

                            indexinbinaire = indexinbinaire + 1;
                        }
                    }

                    if (useless == 0) {
                        // second step: reduce lengths LL=>L, RR=>R
                        if ((int)word.size() > 1) {
                            int i = (int)word.size() - 2;
                            while (i >= 0) {
                                if (word[i] == word[i + 1] && used[i] == 0 && used[i + 1] == 0) {
                                    word.erase(word.begin() + (i + 1)); // pop(i+1)
                                    used[i] = 1; // in this second phase we start from the end...
                                }
                                i = i - 1;
                            }
                        }

                        // our new word is computed
                        std::vector<int> octemp = translate_word_to_multiplicities(word);
                        int newscore = lbOneWay(octemp);

                        if (newscore + 1 < score) {
                            if (newscore < score - 2) {
                                std::cout << "diff2\n";
                            }

                            std::vector<int> oc1 = oc0;
                            reverse(oc1.begin(), oc1.end());
                            if (oc1 < oc0) {
                                oc0 = oc1;
                            }

                            int pos = binary_search_in_special(special, oc0);
                            if (pos < 0) {
                                SpecialEntry e;
                                e.oc0 = oc0;
                                e.flips = newscore + 1;
                                e.octemp = octemp;

                                // store binaire as explicit std::vector<int> (0/1)
                                e.binaire.assign(bits, 0);
                                for (int t = 0; t < bits; ++t) {
                                    e.binaire[t] = (int)((mask >> (bits - 1 - t)) & 1LL);
                                }

                                special.insert(special.begin() + (-pos - 1), std::move(e));
                            } else {
                                int previousnumberofflips = special[pos].flips;
                                if (previousnumberofflips > newscore + 1) {
                                    special[pos].oc0 = oc0;
                                    special[pos].flips = newscore + 1;
                                    special[pos].octemp = octemp;
                                    special[pos].binaire.assign(bits, 0);
                                    for (int t = 0; t < bits; ++t) {
                                        special[pos].binaire[t] = (int)((mask >> (bits - 1 - t)) & 1LL);
                                    }
                                }
                            }
                        }

                        if (newscore + 1 == h) {
                            break;
                        }
                    }
                }
            }
        }
    }

    // special.sort(key=lambda v: sum(v[0]))
    sort(special.begin(), special.end(),
         [](const SpecialEntry& a, const SpecialEntry& b) {
             int sa = std::accumulate(a.oc0.begin(), a.oc0.end(), 0);
             int sb = std::accumulate(b.oc0.begin(), b.oc0.end(), 0);
             return sa < sb;
         });

    return special;
}

// -----------------------------------------------------------------------------
int main(int argc, char **argv) {
    std::string filename = "lb.txt";
    std::cout << "Removing " << filename << std::endl;
    unlink(filename.c_str());
    int c = 27; // we deal with segments having at most c crossings

    if(argc >= 2)
      c = atoi(argv[1]);

    std::cout << "c = " << c << std::endl;

    for (int k = 3; k < c; ++k) {
        double t0 = elapsed();
        std::cout << std::endl << elapsed() << "s ";
        std::cout << "Building segments with " << (k + 1)
             << " crossings (" << k << " initial letters => "
             << (k - 2) << " letters in " << filename << ")\n";

        load_lbfile(filename, true);

        int b0 = 1 << (k - 1);
        int b1 = 1 << k;

        std::vector<SpecialEntry> special = bfs(b0, b1);

        // spe = [(s[0], s[1]) for s in special]
        std::vector<std::pair<std::vector<int>, int>> spe;
        spe.reserve(special.size());
        for (const auto& s : special) {
            spe.push_back({s.oc0, s.flips});
        }

        // append to lb.txt
        {
            std::ofstream f(filename, std::ios::app);
            for (const auto& name : spe) {
                f << "[";
                for (size_t i = 0; i < name.first.size(); ++i) {
                    f << name.first[i];
                    if (i + 1 < name.first.size()) f << ", ";
                }
                f << "] = " << name.second << " \n";
            }
        }

        std::cout << spe.size() << " special cases\n";
        load_lbfile(filename, true);
        std::cout << "Iteration time: " <<elapsed() - t0 << std::endl;
    }

    return 0;
}
