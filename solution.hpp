#pragma once
#include "tools.hpp"
#include "point.hpp"
#include "triangulation.hpp"
#include "instance.hpp"

template <class Number>
class Solution {
  using Path = std::vector<Triangulation<Number>>;
  Instance<Number> *inst;
  std::vector<Path> paths;
  Triangulation<Number> center;
  int length = 0;
  inline static std::map<float,int> log;

public:
  Solution(Instance<Number> &_inst)
    : inst(&_inst), paths(_inst.triangulations.size()) {
  }

  Solution(Instance<Number> &_inst, const std::vector<std::vector<std::vector<IEdge>>> &flips)
    : inst(&_inst), paths(_inst.triangulations.size()) {
    int t = paths.size();
    for(int i = 0; i < t; i++) {
      paths[i].push_back(inst->triangulations[i]);
      for(const auto & pflip : flips.at(i)) {
        auto triangulation = paths[i].back();
        for(IEdge e : pflip) {
          triangulation.flip(e);
        }
        paths[i].push_back(triangulation);
      }
      std::reverse(paths[i].begin(),paths[i].end());
    }

    center = paths[0].front();

    calculateLength();
  }

  Instance<Number> &getInstance() const {
    return *inst;
  }

  const Triangulation<Number> &getCenter() const {
    return center;
  }

  int getLength() const {
    return length;
  }

  std::vector<int> getLengths() {
    std::vector<int> lengths(paths.size());
    for(size_t i = 0; i < lengths.size(); i++) {
      lengths[i] = getPath(i).size() -1;
    }
    return lengths;
  }

  const Path &getPath(int i) const {
    return paths[i];
  }

  std::string makeFilename(std::string label = "") {
    std::string fn = inst->name + "." + timeString(beginTime) + label + ".sol.json";
    return fn;
  }

  std::vector<int> sortedLastFlips() const {
    int t = paths.size();
    std::vector<int> myLengths(t);
    for(int i = 0; i < t; i++)
      myLengths[i] = getPath(i).size() > 1 ? paths[i][0].countFlips(paths[i][1]) : length;
    std::sort(myLengths.begin(), myLengths.end());

    return myLengths;
  }

  std::vector<std::pair<int,int>> decoratedSortedLastFlips() const {
    int t = paths.size();
    std::vector<std::pair<int,int>> myLengths(t);
    for(int i = 0; i < t; i++)
      myLengths[i] = {getPath(i).size() > 1 ? paths[i][0].countFlips(paths[i][1]) : length, i};
    std::shuffle(myLengths.begin(), myLengths.end(), rgen);
    std::stable_sort(myLengths.begin(), myLengths.end());

    return myLengths;
  }

  auto operator<(const Solution<Number> &other) const {
    if(getLength() < other.getLength())
      return true;
    if(getLength() > other.getLength())
      return false;

    int t = paths.size();
    std::vector<int> myLengths = sortedLastFlips();
    std::vector<int> otherLengths = other.sortedLastFlips();

    return myLengths < otherLengths;
  }

  bool stepTo(int r) {
    if(getPath(r).size() <= 1)
      return false;
    center = getPath(r)[1];
    auto edges = center.getEdges();

    for(int i = 0; i < int(paths.size()); i++) {
      if(i == r) {
        paths[i].erase(paths[i].begin()); // Remove front
      }
      else {
        paths[i].insert(paths[i].begin(), center); // Insert on front
      }
    }
    length += paths.size() - 2;

    return true;
  }

  void setCenter(const Triangulation<Number> &newCenter) {
    center = newCenter;
    paths = std::vector<Path>(paths.size());
    length = 0;
  }

  void calculateCenter() {
    center = paths.at(0).at(0);
  }

  void setPath(int c, const Path &newPath) {
    length += (int)newPath.size() - (int)paths[c].size();
    paths[c] = newPath;
  }


  int improvePath(int c, const Path &newPath) {
    if(paths[c].empty()) {
      length += newPath.size() - 1;
      paths[c] = newPath;
      return 2;
    }

    if(paths[c].size() > newPath.size()) {
      length += newPath.size() - paths[c].size();
      paths[c] = newPath;
      return 2;
    }

    if(paths[c].size() < newPath.size())
      return 0;

    if(paths[c].size() >= 2) {
      int fn = (paths[c][0] - paths[c][1]).size();
      int newfn = (newPath[0] - newPath[1]).size();
      if(fn > newfn) {
        paths[c] = newPath;
        return 1;
      }
    }

    return false;
  }

  bool improve(Solution <Number> &newSol, bool autosave=true) {
    if(getLength() == 0 || newSol < *this) {
      *this = newSol;
      if(autosave) save();
      return true;
    }

    return false;
  }

  std::vector<std::vector<IEdge>> getFlipVector(int c) {
    std::vector<std::vector<IEdge>> ret;
    for(int i = paths[c].size()-1; i >= 1; i--) {
      ret.push_back(paths[c][i] - paths[c][i-1]);
    }
    return ret;
  }

  std::vector<std::vector<std::vector<IEdge>>> getFlipVectors() {
    const int t = inst->triangulations.size();
    std::vector<std::vector<std::vector<IEdge>>> ret;

    for(int c = 0; c < t; c++) {
      ret.push_back(getFlipVector(c));
    }
    return ret;
  }

  void save(std::string fn = "", std::string label = "") {
    if(fn.empty()) fn = makeFilename(label);
    std::cout << "Saving solution of length " << length << " to " << fn << std::endl;

    for(Path &path : paths)
      assert(path.size() >= 1);

    std::ofstream f(fn);
    f << "{\n";

    f << "  \"content_type\": \"CGSHOP2026_Solution\",\n";
    f << "  \"instance_uid\": \"" << inst->name << "\",\n";

    f << "  \"objective_value\": " << length << ",\n";

    f << "  \"meta\": {\n";

    for(const auto &[k,v] : inst->meta)
      f << "    \"" << k << "\": \"" << v << "\",\n";

    f << "    \"saved\": [\n";
    for(auto [t,l] : log)
      f << "      [" << t << ", " << l << "],\n";
    f << "      [" << elapsed() << ", " << length << "]\n";
    f << "    ]\n";
    f << "  },\n";

    log[elapsed()] = length;

    const int t = inst->triangulations.size();

    f << "  \"flips\": [\n";
    for(int c = 0; c < t; c++) {
      f << "    [\n";
      for(int i = paths[c].size()-1; i >= 1; i--) {
        f << "      [ ";
        auto flips = paths[c][i] - paths[c][i-1];
        for(int j = 0; j < int(flips.size()); j++) {
          int u = flips[j][0];
          int v = flips[j][1];
          f << "[" << u << "," << v << "]";
          if(j != int(flips.size()) - 1) f << ",";
        }
        f << " ]" << (i != 1 ? "," : "") << "\n";
      }
      f << "    ]" << (c != t-1 ? "," : "") << "\n";
    }
    f << "  ]\n";

    f << "}\n";
  }

  void saveOld(std::string fn = "") {
    if(fn.empty()) fn = makeFilename();
    std::cout << "Saving solution of length " << length << " to " << fn << std::endl;

    for(Path &path : paths)
      assert(path.size() >= 1);

    const int n = inst->points.size();
    std::ofstream f(fn);
    f << "{\n";

    f << "  \"instance_uid\": \"" << inst->name << "\",\n";

    f << "  \"length\": " << length << ",\n";

    f << "  \"log\": [\n";
    for(auto [t,l] : log)
      f << "    [" << t << ", " << l << "],\n";
    f << "    [" << elapsed() << ", " << length << "]\n";
    f << "],\n";

    log[elapsed()] = length;

    f << "  \"points_x\": [";
    for(int i = 0; i < n; i++) {
      f << inst->points[i].x();
      if(i != n-1) f << ",";
    }
    f << "],\n";

    f << "  \"points_y\": [";
    for(int i = 0; i < n; i++) {
      f << (*center.points)[i].y();
      if(i != n-1) f << ",";
    }
    f << "],\n";

    const int t = inst->triangulations.size();
    f << "  \"triangulations\": [\n";
    for(int i = 0; i < t; i++) {
      f << "    [ ";
      auto edges = paths[i].back().getEdges();
      for(int j = 0; j < int(edges.size()); j++) {
        int u = edges[j][0];
        int v = edges[j][1];
        f << "[" << u << "," << v << "]";
        if(j != int(edges.size()) - 1) f << ",";
      }
      f << " ]" << (i != t-1 ? "," : "") << "\n";
    }
    f << "],\n";

    f << "  \"target\": [ ";
    auto centerEdges = center.getEdges();
    for(int j = 0; j < int(centerEdges.size()); j++) {
      int u = centerEdges[j][0];
      int v = centerEdges[j][1];
      f << "[" << u << "," << v << "]";
      if(j != int(centerEdges.size()) - 1) f << ",";
    }
    f << " ],\n";

    f << "  \"flips\": [\n";
    for(int c = 0; c < t; c++) {
      f << "    [\n";
      for(int i = paths[c].size()-1; i >= 1; i--) {
        f << "      [ ";
        auto flips = paths[c][i] - paths[c][i-1];
        for(int j = 0; j < int(flips.size()); j++) {
          int u = flips[j][0];
          int v = flips[j][1];
          f << "[" << u << "," << v << "]";
          if(j != int(flips.size()) - 1) f << ",";
        }
        f << " ]" << (i != 1 ? "," : "") << "\n";
      }
      f << "    ]" << (c != t-1 ? "," : "") << "\n";
    }
    f << "  ]\n";

    f << "}\n";
  }

  Instance<Number> subInstance(int dist) {
    std::vector<int> lengths(paths.size());
    for(int i = 0; i < int(paths.size()); i++) {
      lengths[i] = std::min(dist, (int)paths[i].size()-1);
    }
    return subInstance(lengths);
  }

  Instance<Number> subInstance(std::vector<int> lengths) {
    Instance<Number> sub = *inst;
    for(int i = 0; i < int(paths.size()); i++) {
      sub.triangulations[i] = paths[i].at(lengths[i]);
    }
    sub.name += "_sub";
    return sub;
  }

  void improveFromSub(const Solution <Number> &subSol) {
    for(int i = 0; i < int(paths.size()); i++) {
      auto oldPath = paths[i];
      paths[i] = subSol.paths[i];
      for(auto &t :paths[i])
        t.points = &inst->points;
      auto edges = paths[i].back().getEdges();
      bool found = false;
      // Copy the rest of the path
      for(auto &t :oldPath) {
        if(found)
          paths[i].push_back(t);
        else if(t.getEdges() == edges)
          found = true;
      }
    }
    calculateCenter();
    calculateLength();
  }

  void calculateLength() {
    length = 0;
    for(const auto &path : paths)
      if(!path.empty())
        length += path.size() - 1;
  }
};

template <class Number>
inline std::vector<std::vector<std::vector<IEdge>>> readSolutionFlips(const std::string &fn, std::string &uid) {
  rapidjson::Document doc = Instance<Number>::read_json(fn);
  uid = doc["instance_uid"].GetString();
  std::vector<std::vector<std::vector<IEdge>>> v_flips;
  int length = 0;

  for (auto &path : doc["flips"].GetArray()) {
    std::vector<std::vector<IEdge>> v_path;
    for (auto &pflip: path.GetArray()) {
      std::vector<IEdge> v_pflip;
      for (auto &flip: pflip.GetArray()) {
        v_pflip.push_back({(short)flip[0].GetInt(),(short)flip[1].GetInt()});
      }
      length++;
      v_path.push_back(v_pflip);
    }
    v_flips.push_back(v_path);
  }


  std::cout << "Solution " << fn << " with length " << length << std::endl;

  return v_flips;
}
