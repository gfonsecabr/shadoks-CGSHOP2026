#pragma once
#include "tools.hpp"
#include "triangulation.hpp"
#include "graph.hpp"
#include "independent.hpp"

template <class Number>
class DistanceSolver {
  using Path = std::vector<Triangulation<Number>>;
  Triangulation<Number> source;
  Triangulation<Number> target;
  Path solution;
  bool reversed = false;
  std::unordered_map<IEdge, std::vector<IEdge>> crossings;
  std::unordered_map<IEdge, double> cost, weight;
  double exponent = 1.0;

public:
  DistanceSolver(const Triangulation<Number> &_source, const Triangulation<Number> &_target, const std::unordered_map<IEdge, double> &_weight = {})
    : source(_source), target(_target), weight(_weight) {
    source.clearTrianglesContaining();
    solution.push_back(source);
  }

  void setReverse() {
    reversed = true;
  }

  void setExponent(double _exponent) {
    exponent = _exponent;
  }

  Path getSolution() {
    if(!reversed)
      return solution;
    return Path(solution.rbegin(),solution.rend());
  }

  void greedy() {
    while(greedyStep())
      ;

    assert(solution.back().getEdges() == target.getEdges());
  }

  bool isBetter(Path &p0, Path &p1) {
    if(p1.size() == 0 || p0.size() < p1.size())
      return true;

    if(p0.size() > p1.size() || p0.size() == 0)
      return false;

    if(p0.size() == 1)
      return false;

    if(reversed)
      return p0[0].countFlips(p0[1]) < p1[0].countFlips(p1[1]);

    return p0[p0.size()-1].countFlips(p0[p0.size()-2]) < p1[p1.size()-1].countFlips(p1[p1.size()-2]);
  }

  void squeaky(int maxSteps = 16, int goal = 2) {
    Path bestSolution;

    goal = std::max(goal, 2);

    int i = 0;
    do {
      i++;
      solution = {source};
      cost.clear();
      greedy();

      if(verbose >= 4) std::cout << " " << solution.size()-1 << std::flush;

      if(solution.size() >= 2) {
        for(IEdge e : solution[solution.size()-2] - solution[solution.size()-1]) {
          for(IEdge ce : getCrossings(e)) {
            weight[ce] = 1.2*getWeight(ce);
          }
        }
      }

      if(isBetter(solution, bestSolution))
        bestSolution = solution;

      if(int(solution.size()) - 1 <= goal)
        break;
    } while(solution.size() <= bestSolution.size() && i < maxSteps);

    solution = bestSolution;
    if(verbose >= 2 && verbose < 4) std::cout << solution.size()-1 << std::flush;
  }

  bool greedyStep() {
    auto triangulation = solution.back();
    std::unordered_map<IEdge,IEdge> flips = getFlips(triangulation);
    std::unordered_map<IEdge, double> gain = calculateGains(flips, triangulation);
    Graph<IEdge> graph = incompatibleUsefulFlips(flips,gain);
    IndependentSolver<IEdge,double> solver(graph, gain);
    solver.solve();
    std::vector<IEdge> independent = solver.getSolution();

    // std::unordered_map<IEdge,IEdge> flipsToApply;
    // for(IEdge e0 :independent)
    //   flipsToApply[e0] = flips[e0];
    // std::string fn = "flips" + std::to_string(solution.size()) + ".svg";
    // triangulation.writeSVG(fn,flipsToApply);

    if(independent.empty())
      return false;

    for(IEdge e :independent)
      triangulation.flip(e);
    solution.push_back(triangulation);
    return true;
  }


protected:
  std::unordered_map<IEdge,IEdge> getFlipsAndGain(const Triangulation<Number> &triangulation, std::deque<std::pair<int,IEdge>> & levelFlips, std::unordered_map<IEdge, double> & gain) {
    std::unordered_map<IEdge,IEdge> flips;
    for(auto & [l,e] : levelFlips) {
      if(!target.contains(e) && triangulation.contains(e)) {
        auto f = triangulation.canFlip(e);
        if(f) {
          flips[e] = *f;
          gain[e] = levelFlips.size() - l + 1;
        }
      }
    }
    return flips;
  }


  std::unordered_map<IEdge,IEdge> getFlips(const Triangulation<Number> &triangulation) {
    std::unordered_map<IEdge,IEdge> flips;
    std::vector<IEdge> allEdges = triangulation.getEdges(false);
    for(auto & e : allEdges) {
      if(!target.contains(e)) {
        auto f = triangulation.canFlip(e);
        if(f) {
          flips[e] = *f;
        }
      }
    }
    return flips;
  }

  std::unordered_map<IEdge, double> calculateGains(const std::unordered_map<IEdge,IEdge> &flips, const Triangulation<Number> &triangulation) {
    std::unordered_map<IEdge, double> gains;
    for(auto &[e0,e1] : flips) {
      double gain = pow(getCost(e0), exponent) - pow(getCost(e1), exponent);
      // if(getCrossings(e1) == 0)
        // gain *= 8;
      gains[e0] = gain;
    }
    return gains;
  }

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
        if(gain.at(v) > 0)
          g.addVertex(v);
      }

      for(int i = 0; i < int(clique.size()) - 1; i++) {
        for(int j = i+1; j < int(clique.size()); j++) {
          if(gain.at(clique[i]) > 0 && gain.at(clique[j]) > 0)
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

  std::vector<IEdge> getCrossings(IEdge e) {
    if(!crossings.contains(e))
      crossings[e] = target.getCrossings(e);
    return crossings.at(e);
  }

  double getCost(IEdge e) {
    if(!cost.contains(e)) {
      double c = 0.0;
      for(auto ce : getCrossings(e)) {
        c += getWeight(ce);
      }
      cost[e] = c;
    }

    return cost.at(e);
  }

};
