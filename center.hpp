#pragma once
#include "tools.hpp"
#include "triangulation.hpp"
#include "graph.hpp"
#include "independent.hpp"

template <class Number>
class CenterSolver {
  Triangulation<Number> solution;
  std::vector<Triangulation<Number>> targets;

  std::unordered_map<IEdge, std::vector<std::vector<IEdge>>> crossings;
  std::unordered_map<IEdge, double> cost, weight;
  double exponent = 1.0;

public:
  CenterSolver(const Triangulation<Number> &_source, const std::vector<Triangulation<Number>> &_targets, const std::unordered_map<IEdge, double> &_weight = {})
    : solution(_source), targets(_targets), weight(_weight) {
    solution.clearTrianglesContaining();
  }

  const Triangulation<Number> &getSolution() {
    return solution;
  }

  void setExponent(double _exponent) {
    exponent = _exponent;
  }

  void greedy() {
    while(greedyStep())
      ;
  }

  bool greedyStep() {
    auto triangulation = solution;
    std::unordered_map<IEdge,IEdge> flips = getFlips(triangulation);
    std::unordered_map<IEdge, double> gain = calculateGains(flips, triangulation);
    Graph<IEdge> graph = incompatibleUsefulFlips(flips, gain);
    IndependentSolver<IEdge,double> solver(graph, gain);
    solver.solve();
    std::vector<IEdge> independent = solver.getSolution();

    if(independent.empty()) {
      int c = 0;
      for(IEdge e :solution.getEdges())
        c += getCrossings(e).size();
      if(verbose >= 4) std::cout << " = " << c << std::endl;
      return false;
    }

    if(verbose >= 4) std::cout << " " << independent.size() << "/" << graph.countVertices() << std::flush;
    for(IEdge e :independent)
      triangulation.flip(e);
    solution = triangulation;
    return true;
  }

protected:
  std::unordered_map<IEdge, double> calculateGains(const std::unordered_map<IEdge,IEdge> &flips, const Triangulation<Number> &triangulation) {
    std::unordered_map<IEdge, double> gains;
    for(auto &[e0,e1] : flips) {
      double gain = getCost(e0) - getCost(e1);
      if(gain > 0.0)
        gains[e0] = gain;
    }
    return gains;
  }

  Graph<IEdge> incompatibleUsefulFlips(const std::unordered_map<IEdge,IEdge> &flips, std::unordered_map<IEdge, double> &gain) const {
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
      for(auto v : clique) {
        if(gain.contains(v))
          g.addVertex(v);
      }

      for(int i = 0; i < int(clique.size()) - 1; i++) {
        for(int j = i+1; j < int(clique.size()); j++) {
          if(gain.contains(clique[i]) && gain.contains(clique[j]))
            g.addEdge(clique[i],clique[j]);
        }
      }
    }

    return g;
  }

  double getWeight(IEdge e) {
    if(weight.contains(e))
      return weight.at(e);
    return 1.0;
  }

  const std::vector<std::vector<IEdge>> &getCrossings(IEdge e) {
    if(!crossings.contains(e)) {
      for(auto &target : targets) {
        crossings[e].push_back(target.getCrossings(e));
      }
    }

    return crossings.at(e);
  }

  double getCost(IEdge e) {
    if(!cost.contains(e)) {
      double c = 0.0;
      for(const auto &cr : getCrossings(e)) {
        double ctarget = 0.0;
        for(auto ce : cr) {
          ctarget += getWeight(ce);
        }
        c += pow(ctarget, exponent);
      }
      cost[e] = c;
    }

    return cost.at(e);
  }

  std::unordered_map<IEdge,IEdge> getFlips(const Triangulation<Number> &triangulation) {
    std::unordered_map<IEdge,IEdge> flips;
    std::vector<IEdge> allEdges = triangulation.getEdges();
    for(auto & e : allEdges) {
      auto f = triangulation.canFlip(e);
      if(f) {
        flips[e] = *f;
      }

    }
    return flips;
  }
};
