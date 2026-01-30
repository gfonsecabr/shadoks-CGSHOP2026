#pragma once
#include "tools.hpp"
#include "graph.hpp"

template <class Vertex, class W>
class IndependentSolver {
  Graph<Vertex> &g;
  std::unordered_map<Vertex, W> &weight;
  std::vector<Vertex> sol;
  std::unordered_set<Vertex> forbidden;

public:
  IndependentSolver(Graph<Vertex> &g, std::unordered_map<Vertex, W> &weight)
    : g(g), weight(weight) {}

  void solve() {
    std::vector<Vertex> vertices = g.vertices();

    std::vector<std::tuple<int,int,int,Vertex>> decorated(vertices.size());
    static std::uniform_int_distribution<int> distribution(65536);
    for(size_t i = 0; i < vertices.size(); i++) {
      decorated[i] = {-weight.at(vertices[i]), g.degree(vertices[i]), distribution(rgen), vertices[i]};
    }
    std::sort(decorated.begin(), decorated.end());

    // std:: cout << "Vertices: " << vertices.size() << std::flush;
    for(auto &t : decorated) {
      Vertex v = std::get<3>(t);
      if(weight.at(v) <= 0) //  && sol.size() > 0
        break;
      if(!forbidden.contains(v)) {
        sol.push_back(v);
        for(Vertex u : g.neighbors(v))
          forbidden.insert(u);
      }
    }
  }

  const std::vector<Vertex> &getSolution() {
    return sol;
  }
};
