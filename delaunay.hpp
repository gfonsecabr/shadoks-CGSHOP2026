#pragma once
#include "tools.hpp"
#include "triangulation.hpp"
#include "graph.hpp"
#include "independent.hpp"
#include "solution.hpp"

template <class Number>
class DelaunaySolver {
  using Path = std::vector<Triangulation<Number>>;
  Instance<Number> *inst;
  Solution<Number> solution;

public:
  DelaunaySolver(Instance<Number> &_inst)
    : inst(&_inst), solution(_inst){
  }

  const Solution<Number> &getSolution() {
    return solution;
  }

  void run() {
    for(size_t t = 0; t < inst->triangulations.size(); t++)
      solveTarget(t);
  }

  void solveTarget(int t) {
    Path path = {inst->triangulations.at(t)};
    path.back().clearTrianglesContaining();
    while(flipStep(path))
      ;

    if(solution.getCenter().empty())
      solution.setCenter(path.back());

    std::reverse(path.begin(),path.end());
    solution.improvePath(t, path);
  }

  bool flipStep(Path &path) {
    auto triangulation = path.back();
    std::unordered_map<IEdge,IEdge> flips = getFlips(triangulation);

    if(flips.empty()) {
      if(verbose >= 5) std::cout << " = " << path.size() - 1 << std::endl;
      return false;
    }

    std::unordered_map<IEdge, double> gain;
    for(auto &[f,_] : flips) gain[f] = 1;
    Graph<IEdge> graph = incompatibleFlips(flips);
    IndependentSolver<IEdge,double> solver(graph, gain);
    solver.solve();
    std::vector<IEdge> independent = solver.getSolution();

    if(verbose >= 5) std::cout << " " << independent.size() << std::flush;
    for(IEdge e :independent)
      triangulation.flip(e);

    path.push_back(triangulation);
    return true;
  }

protected:
  Graph<IEdge> incompatibleFlips(const std::unordered_map<IEdge,IEdge> &flips) const {
    std::unordered_map<ITriangle,std::vector<IEdge>> flipsUsing;
    for(auto &[e0,e1] : flips) {
      ITriangle t0 = {e0[0],e0[1],e1[0]};
      ITriangle t1 = {e0[0],e0[1],e1[1]};
      std::sort(t0.begin(), t0.end());
      std::sort(t1.begin(), t1.end());
      flipsUsing[t0].push_back(e0);
      flipsUsing[t1].push_back(e0);
    }

    Graph<IEdge> g;
    for(auto &[_,clique] : flipsUsing) {
      for(auto v : clique)
        g.addVertex(v);

      for(int i = 0; i < int(clique.size()) - 1; i++) {
        for(int j = i+1; j < int(clique.size()); j++) {
          g.addEdge(clique[i],clique[j]);
        }
      }
    }

    return g;
  }


  std::unordered_map<IEdge,IEdge> getFlips(const Triangulation<Number> &triangulation) {
    std::unordered_map<IEdge,IEdge> flips;
    std::vector<IEdge> allEdges = triangulation.getEdges();
    for(auto & e : allEdges) {
      auto f = triangulation.canFlip(e);
      if(f) {
        Triangle<Number> tri(inst->points[e[0]], inst->points[e[1]], inst->points[(*f)[0]]);
        int r = tri.incircle(inst->points[(*f)[1]]);
        if(r > 0 || (r == 0 && *f < e))
          flips[e] = *f;
      }
    }
    return flips;
  }
};
