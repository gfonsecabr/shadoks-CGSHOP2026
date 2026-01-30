#pragma once
#include "tools.hpp"
#include "point.hpp"
#include "segment.hpp"
#include "triangle.hpp"
#include "graph.hpp"


template <class T>
void printVectors(const T &v) {
  for(const auto &t : v) {
    std::cout << "{";
    for(const auto &x : t)
      std::cout << x << " ";
    std::cout << "} ";
  }
  std::cout << std::endl;
}

template<class T>
inline void hash_combine(std::size_t& seed, const T& val) {
  std::hash<T> hasher;
  seed ^= hasher(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


using IEdge = std::array<short int,2>;
using ITriangle = std::array<short int,3>;

namespace std {
  template<>
  struct hash<IEdge> {
    inline size_t operator()(const IEdge& val) const {
      return std::hash<std::uint32_t>()(*reinterpret_cast<const std::uint32_t *>(val.data()));
    }
  };
}

namespace std {
  template<>
  struct hash<ITriangle> {
    inline size_t operator()(const ITriangle& val) const {
      size_t seed = 0;
      hash_combine(seed, val[0]);
      hash_combine(seed, val[1]);
      hash_combine(seed, val[2]);
      return seed;
    }
  };
}

template <class Number>
class Triangulation {
public:
  std::vector<Point<Number>> *points;

protected:
  std::unordered_map<IEdge, std::array<ITriangle,2>> adj;
  mutable std::vector<std::vector<IEdge>> trianglesContaining;

public:
  Triangulation() {}

  Triangulation(std::vector<Point<Number>> *_points, const std::vector<IEdge> &edges)
  : points(_points) {
    Graph<short int> g(edges);
    std::unordered_map<IEdge, std::vector<ITriangle>> adj2;
    std::vector<ITriangle> cycles = g.triangles();
    for(const ITriangle &t : cycles) {
      const std::array<IEdge,3> tedges{{{t[0],t[1]}, {t[1],t[2]}, {t[0],t[2]}}};
      Triangle<Number> tri = {(*_points)[t[0]],(*_points)[t[1]],(*_points)[t[2]]};

      for(const IEdge &edge : tedges) {
        if(!adj2.contains(edge))
          adj2[edge] = std::vector<ITriangle>{t};
        else {
          int p = missing(t, edge);
          auto &aedge = adj2.at(edge);
          bool empty = true;
          for(auto &t2: aedge) {
            int p2 = missing(t2, edge);
            if(tri.contains((*_points)[p2])) {
              empty = false;
              break;
            }
          }
          if(empty) {
            std::erase_if(aedge,[_points,p](ITriangle t2){
              Triangle<Number> tri2 = {(*_points)[t2[0]],(*_points)[t2[1]],(*_points)[t2[2]]};
              return tri2.contains((*_points)[p]);
            });
            aedge.push_back(t);
          }
        }
      }
    }

    for(auto &[edge,tv] : adj2) {
      // std::cout << edge[0] << "," << edge[1] << ": ";
      // printVectors(tv);
      assert(tv.size() <= 2);
      assert(tv.size() > 0);
      for(ITriangle t : tv) {
        if(adj.contains(edge)) {
          adj[edge][1] = t;
        }
        else {
          adj[edge] = {t, {-1,-1,-1}};
        }
      }
    }
  }

  bool empty() const {
    return adj.empty();
  }

  bool operator!=(const Triangulation<Number> &other) const {
    return getEdges() != other.getEdges();
  }

  bool operator==(const Triangulation<Number> &other) const {
    return getEdges() == other.getEdges();
  }

  bool contains(IEdge e) const {
    return adj.contains(e);
  }

  void addTriangle(ITriangle t) {
    assert(t[0] < t[1] && t[1] < t[2]);
    const std::array<IEdge,3> edges{{{t[0],t[1]}, {t[1],t[2]}, {t[0],t[2]}}};
    for(const IEdge &edge : edges) {
      if(adj.contains(edge)) {
        adj[edge][1] = t;
      }
      else {
        adj[edge] = {t, {-1,-1,-1}};
      }
    }
  }

  std::optional<IEdge> canFlip(IEdge oldEdge) const {
    std::array<ITriangle,2> oldTriangles = adj.at(oldEdge);
    if(oldTriangles[1][0] == -1)
      return {};

    IEdge newEdge{missing(oldTriangles[0],oldEdge), missing(oldTriangles[1],oldEdge)};
    if(newEdge[0] > newEdge[1])
      std::swap(newEdge[0],newEdge[1]);
    ITriangle ti1 = {newEdge[0], newEdge[1], oldEdge[0]};
    ITriangle ti2 = {newEdge[0], newEdge[1], oldEdge[1]};
    Triangle<Number> t1 = {(*points)[ti1[0]], (*points)[ti1[1]], (*points)[ti1[2]]};
    Triangle<Number> t2 = {(*points)[ti2[0]], (*points)[ti2[1]], (*points)[ti2[2]]};

    if(t1.canFlip(t2))
      return newEdge;

    return {};
  }

  void flip(IEdge oldEdge) {
    assert(trianglesContaining.empty());
    std::array<ITriangle,2> oldTriangles = adj.at(oldEdge);
    IEdge newEdge{missing(oldTriangles[0],oldEdge), missing(oldTriangles[1],oldEdge)};
    if(newEdge[0] > newEdge[1])
      std::swap(newEdge[0],newEdge[1]);

    adj.erase(oldEdge);
    ITriangle ti1 = {newEdge[0], newEdge[1], oldEdge[0]};
    ITriangle ti2 = {newEdge[0], newEdge[1], oldEdge[1]};
    std::sort(ti1.begin(), ti1.end());
    std::sort(ti2.begin(), ti2.end());
    adj[newEdge] = {ti1,ti2};

    for(const ITriangle &t : adj.at(newEdge)) {
      for(IEdge e : std::array<IEdge,3>{{ {t[0],t[1]}, {t[1],t[2]}, {t[0],t[2]} }}) {
        if(e != newEdge) {
          replaceTriangle(e, oldEdge , t);
        }
      }
    }
  }

  std::unordered_set<ITriangle> getTriangles() const {
    std::unordered_set<ITriangle> ret;
    for(const auto &[_,t0t1] : adj) {
      ret.insert(t0t1[0]);
      if(t0t1[1][0] >= 0)
        ret.insert(t0t1[1]);
    }
    return ret;
  }

  std::vector<IEdge> getEdges(bool sorted = true) const {
    std::vector<IEdge> ret;
    for(const auto &[e,_] : adj) {
      ret.push_back(e);
    }
    if(sorted)
      std::sort(ret.begin(),ret.end());
    return ret;
  }

  void printTriangles() const {
    auto triangles = getTriangles();
    for(const ITriangle &t :triangles)
      std::cout << "(" << t[0] << "," << t[1] << "," << t[2] << ") " ;
  }

  std::vector<IEdge> flipsTo(const Triangulation &triangulation) const {
    std::vector<IEdge> flips;
    for(const auto &[e,_] : adj) {
      if(!triangulation.adj.contains(e))
        flips.push_back(e);
    }
    return flips;
  }

  std::vector<IEdge> operator-(const Triangulation &triangulation) const {
    return flipsTo(triangulation);
  }

  int countFlips(const Triangulation &triangulation) const {
    int ret = 0;
    for(const auto &[e,_] : adj) {
      if(!triangulation.adj.contains(e))
        ret++;
    }
    return ret;
  }

  void clearTrianglesContaining() {
    trianglesContaining.clear();
  }

  void buildTrianglesContaining() const {
    trianglesContaining.clear();
    trianglesContaining.resize(Triangulation<Number>::points->size());
    auto triangles = Triangulation<Number>::getTriangles();
    for(const ITriangle &t : triangles) {
      trianglesContaining[t[0]].push_back({t[1],t[2]});
      trianglesContaining[t[1]].push_back({t[0],t[2]});
      trianglesContaining[t[2]].push_back({t[0],t[1]});
    }
  }

  int countCrossings(IEdge e0) {
    return getCrossings().size();
  }

  bool cross(IEdge e0, const std::unordered_set<IEdge> &forbidden) {
    if(trianglesContaining.empty())
      buildTrianglesContaining();
    Segment<Number> s0((*points)[e0[0]],(*points)[e0[1]]);

    const short startVertex = e0[0], endVertex = e0[1];
    ITriangle startTriangle = {-1,-1,-1};
    IEdge crossingEdge;
    for(IEdge e :trianglesContaining.at(startVertex)) {
      Segment<Number> s((*points)[e[0]],(*points)[e[1]]);
      if(s0.cross(s)) {
        startTriangle = {startVertex,e[0],e[1]};
        std::sort(startTriangle.begin(), startTriangle.end());
        crossingEdge = e;
        break;
      }
    }

    if(startTriangle[0] == -1) // No triangle is crossed
      return false;

    if(forbidden.contains(crossingEdge))
      return true;
    ITriangle currentTriangle = startTriangle;
    while(currentTriangle[0] != endVertex && currentTriangle[1] != endVertex && currentTriangle[2] != endVertex) {
      auto &trianglePair = adj.at(crossingEdge);
      currentTriangle = (trianglePair[0] != currentTriangle ? trianglePair[0] : trianglePair[1]);
      for(IEdge e : std::array<IEdge,3>{{{currentTriangle[0],currentTriangle[1]},
                                         {currentTriangle[1],currentTriangle[2]},
                                         {currentTriangle[0],currentTriangle[2]}}}) {
        if(e != crossingEdge) {
          Segment<Number> s((*points)[e[0]],(*points)[e[1]]);
          if(s0.cross(s)) {
            crossingEdge = e;
            if(forbidden.contains(crossingEdge))
              return true;
            break;
          }
        }
      }
    }

    return false;
  }

  std::vector<IEdge> getCrossings(IEdge e0, int limit = INT_MAX) const {
    // if(contains(e0))
    //   return {};
    if(trianglesContaining.empty())
      buildTrianglesContaining();
    std::vector<IEdge> ret;
    Segment<Number> s0((*points)[e0[0]],(*points)[e0[1]]);

    const short startVertex = e0[0], endVertex = e0[1];
    ITriangle startTriangle = {-1,-1,-1};
    IEdge crossingEdge;
    for(IEdge e :trianglesContaining.at(startVertex)) {
      Segment<Number> s((*points)[e[0]],(*points)[e[1]]);
      if(s0.cross(s)) {
        startTriangle = {startVertex,e[0],e[1]};
        std::sort(startTriangle.begin(), startTriangle.end());
        crossingEdge = e;
        break;
      }
    }

    if(startTriangle[0] == -1) // No triangle is crossed
      return {};

    ret.push_back(crossingEdge);
    ITriangle currentTriangle = startTriangle;
    while(currentTriangle[0] != endVertex && currentTriangle[1] != endVertex && currentTriangle[2] != endVertex) {
      auto &trianglePair = adj.at(crossingEdge);
      currentTriangle = (trianglePair[0] != currentTriangle ? trianglePair[0] : trianglePair[1]);
      for(IEdge e : std::array<IEdge,3>{{{currentTriangle[0],currentTriangle[1]},
                                         {currentTriangle[1],currentTriangle[2]},
                                         {currentTriangle[0],currentTriangle[2]}}}) {
        if(e != crossingEdge) {
          Segment<Number> s((*points)[e[0]],(*points)[e[1]]);
          if(s0.cross(s)) {
            crossingEdge = e;
            ret.push_back(crossingEdge);
            if(int(ret.size()) > limit)
              return ret;
            break;
          }
        }
      }
    }

    return ret;
  }

  std::vector<IEdge> slowGetCrossings(IEdge e0) const {
    std::vector<IEdge> ret;

    Segment<Number> s0((*Triangulation<Number>::points)[e0[0]], (*Triangulation<Number>::points)[e0[1]]);
    for(const auto &[e1,_] : Triangulation<Number>::adj) {
      Segment<Number> s1((*Triangulation<Number>::points)[e1[0]], (*Triangulation<Number>::points)[e1[1]]);
      if(s0.cross(s1))
        ret.push_back(e1);
    }

    return ret;
  }

  void writeSVG(const std::string &fn, const std::unordered_map<IEdge,IEdge> &flips = {}, std::unordered_set<IEdge> red = {}) const {
    Number sz = 0;
    for(auto p : *points) {
      sz = std::max(p.x(),sz);
      sz = std::max(p.y(),sz);
    }

    std::ofstream f(fn);
    std::string prefix = R"(<?xml version="1.0" encoding="utf-8"?>
  <svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="1000.0" height="1000.0">)";
    std::string suffix = R"(
  </svg>
  )";

    f << prefix;

    for(auto [ie0,ie1] : flips) {
      std::array<Point<Number>,4> s = {(*points)[ie0[0]], (*points)[ie1[0]], (*points)[ie0[1]], (*points)[ie1[1]]};

      f << R"(<polygon stroke="none" fill="yellow" points=")";
      for(size_t j = 0; j < s.size(); j++) {
        f << 1000*s[j].x()/sz << " " << 1000*s[j].y()/sz;
        if(j != s.size() -1)
          f << ' ';
      }
      f << "\"></polygon>\n";
    }

    for(IEdge ie : getEdges()) {
      std::array<Point<Number>,2> s = {(*points)[ie[0]], (*points)[ie[1]]};

      f << R"(<polyline stroke="black" points=")";
      for(size_t j = 0; j < s.size(); j++) {
        f << 1000*s[j].x()/sz << " " << 1000*s[j].y()/sz;
        if(j != s.size() -1)
          f << ' ';
      }
      f << "\"></polyline>\n";
    }

    for(IEdge ie : red) {
      std::array<Point<Number>,2> s = {(*points)[ie[0]], (*points)[ie[1]]};

      f << R"(<polyline stroke="red" stroke-opacity=".5" stroke-width="3" points=")";
      for(size_t j = 0; j < s.size(); j++) {
        f << 1000*s[j].x()/sz << " " << 1000*s[j].y()/sz;
        if(j != s.size() -1)
          f << ' ';
      }
      f << "\"></polyline>\n";
    }

    for(auto [ie0,ie1] : flips) {
      std::array<Point<Number>,2> s = {(*points)[ie0[0]], (*points)[ie0[1]]};

      f << R"(<polyline stroke="red" stroke-dasharray="4" points=")";
      for(size_t j = 0; j < s.size(); j++) {
        f << 1000*s[j].x()/sz << " " << 1000*s[j].y()/sz;
        if(j != s.size() -1)
          f << ' ';
      }
      f << "\"></polyline>\n";

      s = {(*points)[ie1[0]], (*points)[ie1[1]]};
      f << R"(<polyline stroke="green" stroke-dasharray="4" points=")";
      for(size_t j = 0; j < s.size(); j++) {
        f << 1000*s[j].x()/sz << " " << 1000*s[j].y()/sz;
        if(j != s.size() -1)
          f << ' ';
      }
      f << "\"></polyline>\n";
    }

    // f << "<text y=\"20\">" << val << "</text>\n";

    f << suffix;
  }

  size_t hash() const {
    std::vector<IEdge> vec = getEdges();
    std::size_t seed = vec.size();
    for(auto e : vec) {
      std::size_t x = std::hash<IEdge>()(e);
      seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }


protected:
  void replaceTriangle(IEdge e, IEdge oldEdge , ITriangle t) {
    auto &trianglePair = adj.at(e);
    if(oldEdge == IEdge{trianglePair[0][0], trianglePair[0][1]} ||
       oldEdge == IEdge{trianglePair[0][1], trianglePair[0][2]} ||
       oldEdge == IEdge{trianglePair[0][0], trianglePair[0][2]})
      trianglePair[0] = t;
    else
      trianglePair[1] = t;
  }
};









