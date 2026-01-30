#pragma once
#include "segment.hpp"
#include "tools.hpp"
#include "triangulation.hpp"
#include "cadicalinterface.h"
#include "EvalMaxSAT.h"
#include "instance.hpp"
#include "crosslb.hpp"

template <class Number>
class SATDistanceSolver {
  using Clause = std::vector<int>;
  using Path = std::vector<Triangulation<Number>>;
  Triangulation<Number> source;
  Triangulation<Number> target;
  Instance<Number> *inst;
  Path solution;
  std::map<std::tuple<IEdge,int>, int> edgeVars;
  std::map<std::tuple<IEdge,IEdge,int>, int> flipVars;
  std::vector<Clause> clauses;
  std::unordered_set<IEdge> happyEdges;
  bool inexact, happyConjecture ;
  int nvars;
  std::vector<Point<Number>> sortedPoints;

public:
  SATDistanceSolver(Instance<Number> &_inst, const Triangulation<Number> &_source, const Triangulation<Number> &_target, bool _inexact, bool _happyConjecture)
    : source(_source), target(_target), inst(&_inst), inexact(_inexact), happyConjecture(_happyConjecture) {
    if(inexact)
      CrossLB::load_lbfile("lb.txt");
  }

  int getLength() {
    return int(solution.size()) - 1;
  }

  Path getReverseSolution() {
    return Path(solution.rbegin(),solution.rend());
  }

  Path getSolution() {
    return solution;
  }

  bool solveDecreasingDistanceMaxSAT(int k) {
    bool solved = false;
    int ret = 0;
    for(int l = k; l >= 0; l--) {
      if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Looking for a path of length " << l << std::endl;
      ret = solveExactDistanceMaxSAT(l);
      if(ret == 1) {
        if(verbose >= 3) {
          std::cout << int(elapsed()) << "s: FOUND a path of length " <<  solution.size() - 1;
          if(solution.size() >= 2) {
            int fn = (solution.back() - solution[solution.size()-2]).size();
            std::cout << " with " << fn << " last step flips" << std::endl;
          }
          else
            std::cout << std::endl;
        }
        solved = true;
      }
      else {
        if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "IMPOSSIBLE!\n";
        break;
      }

      if(l == solution.size()-1) // MaxSAT will find a shorter solution if it exists
        break; // If no shorter solution exists, stop
    }

    return solved;
  }

  bool solveDecreasingDistanceSAT(int k) {
    bool solved = false;
    int ret = 0;
    for(int l = k; l >= 0; l--) {
      if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Looking for a path of length " << l << std::endl;
      ret = solveExactDistanceSAT(l);
      if(ret == 1) {
        if(verbose >= 3) {
          std::cout << int(elapsed()) << "s: FOUND a path of length " <<  solution.size() - 1;
          if(solution.size() >= 2) {
            int fn = (solution.back() - solution[solution.size()-2]).size();
            std::cout << " with " << fn << " last step flips" << std::endl;
          }
          else
            std::cout << std::endl;
        }
        solved = true;
      }
      else {
        if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "IMPOSSIBLE!\n";
        break;
      }
    }

    return solved;
  }

  int solveExactDistanceMaxSAT(int k, int delta = 1) {
    try {
      buildClauses(k);
    }
    catch(std::out_of_range &e) {
      return false;
    }

    int satisfiable = solveMaxSAT(k, delta);

    return satisfiable;
  }

  bool solveExactDistanceIterativeMaxSAT(int k, int delta) {
    Path partial = {target};
    for(int i = 0; i < delta; i++) {
      if(!solveExactDistanceMaxSAT(k-i)) {
        std::cout << "No solution found!\n";
        if(inexact) {
          std::cout << "Disabling  inexact formulation and retrying\n";
          i--;
          inexact = false;
          continue;
        }
        if(i == 0)
          return false;
        else
          break;
      }
      // else
        // std::cout << "solved " << k-i << " " << solution.size() << std::endl;
      if(solution.size() >= 2) {
        target = solution[solution.size()-2];
        partial.push_back(target);
        if(verbose >= 3) std::cout << int(elapsed()) << "s: " << (partial.back()-partial[partial.size()-2]).size() << " flips\n";
        if(solution.size() == 2)
          break;
      }
      else {
        std::cout << "Weird!\n";
        exit(1);
      }
    }

    auto link = partial.back();
    partial.pop_back();
    auto linkEdges = link.getEdges();

    while(!solution.empty() && solution.back().getEdges() != linkEdges)
      solution.pop_back();

    while(!partial.empty()) {
      solution.push_back(partial.back());
      partial.pop_back();
    }

    target = solution.back();
    // std::cout << "SOLVED " << k << " " << solution.size() << std::endl;
    return true;
  }

  int solveExactDistanceSAT(int k) {
    try {
      buildClauses(k);
    }
    catch(std::out_of_range &e) {
      return false;
    }

    int satisfiable = solveSAT(k);

    return satisfiable;
  }

protected:
  int solveSAT(int k) {
    Solver_cadical solver;

    for(Clause &c :clauses) {
      solver.addClause(c);
    }
    clauses = {}; // To free memory

    int ret = solver.solve();
    if(ret == 1) {
      std::vector<bool> satSol = solver.getSolution();
      std::vector<std::vector<IEdge>> edgesPath(k+1);

      for(auto [et,var] : edgeVars) {
        if(satSol.at(var)) {
          edgesPath.at(std::get<int>(et)).push_back(std::get<IEdge>(et));
        }
      }

      solution.clear();
      for(auto &edges :edgesPath) {
        solution.push_back(Triangulation<Number>(source.points, edges));
      }
    }
    return ret;
  }

  int solveMaxSAT(int k, int minimizeDelta = 1) {
    EvalMaxSAT solver;

    std::vector<int> solverVars(nvars);
    int maxVar = 0;
    for(int &var : solverVars) {
      var = solver.newVar();
      maxVar = std::max(var, maxVar);
    }

    for(Clause &c :clauses) {
      std::vector<int> cSolver;

      for(int lit : c) {
        if(lit < 0)
          cSolver.push_back(-solverVars.at(abs(lit)-1));
        else
          cSolver.push_back(solverVars.at(abs(lit)-1));
      }
      if(!cSolver.empty())
        solver.addClause(cSolver);
    }

    clauses = {}; // To free memory

    // Soft clauses
    int soft = 0;
    for(auto [f,var] : flipVars) {
      int t = std::get<int>(f);
      if(t >= k - minimizeDelta) {
        solver.addClause({-solverVars.at(var-1)}, 1);
        soft++;
        // std::cout << '.';
      }
    }

    if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Soft clauses: " << soft << std::endl;

    if(minimizeDelta > 1) {
      solver.setTargetComputationTime(25); // Target time
      solver.setBoundRefTime(.3, 50.0);
    }
    else {
      solver.setTargetComputationTime(10); // Target time
      solver.setBoundRefTime(0.0, 0.0);
    }
    solver.disableOptimize();

// TODO Set previous solution
//     std::vector<bool> bsol;
//     if(!sol.empty()) {
//       bsol.resize(maxVar+1);
//       for(int v : sol) {
//         bsol[vars[v]] = true;
//       }
//
//       solver.setSolution(bsol, sol.size());
//     }

    // solver.printInfo();
    int ret = solver.solve();

    if(ret == 1 && solver.getCost() < std::numeric_limits<t_weight>::max()) {
      std::vector<std::vector<IEdge>> edgesPath(k+1);

      for(auto [et,var] : edgeVars) {
        if(solver.getValue(solverVars.at(var-1)))
          edgesPath.at(std::get<int>(et)).push_back(std::get<IEdge>(et));
      }

      solution.clear();
      for(auto &edges :edgesPath) {
        solution.push_back(Triangulation<Number>(source.points, edges));
      }
      // If we have 0 flips in the last step
      if(solution.size() >= 2 && solution[solution.size()-1].getEdges() == solution[solution.size()-2].getEdges()) {
        solution.pop_back();
      }

    }
    return ret;
  }

  void buildVariables(int k, bool flipVariables = true) {
    nvars = 0;
    edgeVars.clear();
    flipVars.clear();
    const int n = inst->points.size();

    if(inexact || happyConjecture) {
      happyEdges.clear();
      for(IEdge e : source.getEdges())
        if(target.contains(e))
          happyEdges.insert(e);
      if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Happy edges: " << happyEdges.size() << std::endl;
    }

    std::vector<std::unordered_set<int>> flipableEdges(n);

    for(short u = 0; u < n-1; u++) {
      for(short v = u+1; v < n; v++) {
        if(!emptySegment(u,v))
          continue;

        if((inexact || happyConjecture) && crossHappy({u,v}))
          continue;

        auto xs = source.getCrossings({u,v}, 1<<k);
        auto xt = target.getCrossings({u,v}, 1<<k);

        for(int t = 0; t <= k; t++) {
          // if((t > 10 || xs < (1<<t)) &&
          //   (k-t > 10 || xt < (1<<(k-t)))
          if(CrossLB::simpleCrossLB(xs) <= t && CrossLB::simpleCrossLB(xt) <= k-t
             && (!inexact || (CrossLB::fancyCrossLB(xs) <= t && CrossLB::fancyCrossLB(xt) <= k-t))) {
            if(!edgeVars.contains({{u,v},t})) {
              edgeVars[{{u,v},t}] = ++nvars;
              flipableEdges[u].insert(v);
            }
          }
        }
      }
    }

    for(IEdge e :happyEdges) {
      flipableEdges[e[0]].erase(e[1]);
    }

    if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Edge variables: " << edgeVars.size() << std::endl;

    if(!flipVariables)
      return;

    // Variables for flips
    for(short u = 0; u < n; u++) {
      for(short v :flipableEdges[u]) {
        if((inexact || happyConjecture) && target.contains({u,v}))
          continue;
        for(short u2 = 0; u2 < n; u2++) {
          if(u2 == u || u2 == v || flipableEdges[u2].empty() || !inst->emptyTriangle(u,v,u2))
            continue;
          for(short v2 : flipableEdges[u2]) {
            if((inexact || happyConjecture) && source.contains({u2,v2}))
              continue;

            if(v2 != u && v2 != v && convexQuadrilateral(u,u2,v,v2) && inst->emptyTriangle(u,v,v2)) {
              IEdge uu2 = u<u2 ? IEdge{u,u2} : IEdge{u2,u};
              IEdge vu2 = v<u2 ? IEdge{v,u2} : IEdge{u2,v};
              IEdge uv2 = u<v2 ? IEdge{u,v2} : IEdge{v2,u};
              IEdge vv2 = v<v2 ? IEdge{v,v2} : IEdge{v2,v};
              for(int t = 0; t < k; t++) {
                if(edgeVars.contains({{u,v},t}) &&
                    edgeVars.contains({uu2,t}) &&
                    edgeVars.contains({vu2,t}) &&
                    edgeVars.contains({uv2,t}) &&
                    edgeVars.contains({vv2,t}) &&
                    edgeVars.contains({{u2,v2},t+1}) &&
                    edgeVars.contains({uu2,t+1}) &&
                    edgeVars.contains({vu2,t+1}) &&
                    edgeVars.contains({uv2,t+1}) &&
                    edgeVars.contains({vv2,t+1})) {
                  if(!flipVars.contains({{u,v},{u2,v2},t}))
                    flipVars[{{u,v},{u2,v2},t}] = ++nvars;
                }
              }
            }
          }
        }
      }
    }

    if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Flip variables: " << flipVars.size() << std::endl;
  }

  void buildClauses(int k) {
    buildVariables(k);
    clauses.clear();

    if(inexact || happyConjecture) {
      for(IEdge e : happyEdges)
        for(int t = 1; t < k; t++)
          clauses.push_back({edgeVars.at({e,t})});
    }

    // Clauses for start and target edges
    for(IEdge e :source.getEdges())
      clauses.push_back({edgeVars.at({e,0})});
    for(IEdge e :target.getEdges())
      clauses.push_back({edgeVars.at({e,k})});

    // Data structures
    std::vector<std::tuple<int, IEdge, IEdge, int>> indexedFlipVars;
    std::vector<std::tuple<int, IEdge, IEdge, int>> indexedBackFlipVars;
    for(auto [f,var] : flipVars) {
      indexedFlipVars.push_back({std::get<2>(f), std::get<0>(f), std::get<1>(f),var});
      indexedBackFlipVars.push_back({std::get<2>(f), std::get<1>(f), std::get<0>(f),var});
    }
    std::sort(indexedFlipVars.begin(), indexedFlipVars.end());
    std::sort(indexedBackFlipVars.begin(), indexedBackFlipVars.end());

    // Clauses for flips
    for(auto [f,var] : flipVars) {
      short u = std::get<0>(f)[0];
      short v = std::get<0>(f)[1];
      IEdge uv = std::get<0>(f);
      short u2 = std::get<1>(f)[0];
      short v2 = std::get<1>(f)[1];
      IEdge u2v2 = std::get<1>(f);
      int t = std::get<2>(f);
      IEdge uu2 = u<u2 ? IEdge{u,u2} : IEdge{u2,u};
      IEdge vu2 = v<u2 ? IEdge{v,u2} : IEdge{u2,v};
      IEdge uv2 = u<v2 ? IEdge{u,v2} : IEdge{v2,u};
      IEdge vv2 = v<v2 ? IEdge{v,v2} : IEdge{v2,v};

      // All edges except u2v2 must be present at time t
      clauses.push_back({ -var, edgeVars.at({uv,t}) });
      clauses.push_back({ -var, edgeVars.at({uu2,t}) });
      clauses.push_back({ -var, edgeVars.at({vu2,t}) });
      clauses.push_back({ -var, edgeVars.at({uv2,t}) });
      clauses.push_back({ -var, edgeVars.at({vv2,t}) });

      // All edges except uv must be present at time t+1
      // This replaces that the quadrilaterals of two flips cannot intersect
      clauses.push_back({ -var, edgeVars.at({u2v2,t+1}) });
      clauses.push_back({ -var, edgeVars.at({uu2,t+1}) });
      clauses.push_back({ -var, edgeVars.at({vu2,t+1}) });
      clauses.push_back({ -var, edgeVars.at({uv2,t+1}) });
      clauses.push_back({ -var, edgeVars.at({vv2,t+1}) });

      // Edges uv is flipped into u2v2 at time t+1
      if(edgeVars.contains({uv,t+1}))
        clauses.push_back({ -var, -edgeVars.at({uv,t+1}) });
      clauses.push_back({ -var, edgeVars.at({u2v2,t+1}) });
    }

    // Edges cannot disappear without a flip
    for(auto [e,var] : edgeVars) {
      IEdge uv = std::get<0>(e);
      int t = std::get<1>(e);
      if(t == k)
        continue;
      Clause c = {-var};

      if(edgeVars.contains({uv,t+1}))
        c.push_back(edgeVars.at({uv,t+1}));

      auto it = std::upper_bound(indexedFlipVars.begin(), indexedFlipVars.end(), std::make_tuple(t,uv,IEdge{-1,-1},-1));
      for(; it != indexedFlipVars.end() && std::get<0>(*it) == t && std::get<1>(*it) == uv; ++it)
        c.push_back(std::get<3>(*it));

      clauses.push_back(c);
    }


    // Edges cannot appear without a flip
    for(auto [e,var] : edgeVars) {
      IEdge u2v2 = std::get<0>(e);
      int tplus1 = std::get<1>(e);
      if(tplus1 == 0)
        continue;
      Clause c = {-var};

      if(edgeVars.contains({u2v2,tplus1-1}))
        c.push_back(edgeVars.at({u2v2,tplus1-1}));

      auto it = std::upper_bound(indexedBackFlipVars.begin(), indexedBackFlipVars.end(), std::make_tuple(tplus1-1,u2v2,IEdge{-1,-1},-1));
      for(; it != indexedBackFlipVars.end() && std::get<0>(*it) == tplus1-1 && std::get<1>(*it) == u2v2; ++it)
        c.push_back(std::get<3>(*it));

      clauses.push_back(c);
    }

    if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Clauses: " << clauses.size() << std::endl;
  }

  bool crossHappy(IEdge iuv) {
    if(happyEdges.size() > 10*sqrt(inst->points.size())) {
      return target.cross(iuv, happyEdges);
    }

    Segment<Number> uv = {inst->points[iuv[0]], inst->points[iuv[1]]};
    for(IEdge ie :happyEdges) {
      Segment<Number> e = {inst->points[ie[0]], inst->points[ie[1]]};
      if(uv.cross(e))
        return true;
    }

    return false;
  }

  bool emptySegment(int ia, int ib) const {
    Segment<Number> t = {inst->points[ia], inst->points[ib]};
    for(Point<Number> p : inst->points) {
      if(t.containsInside(p))
        return false;
    }

    return true;
  }

  bool convexQuadrilateral(int ia, int ib, int ic, int id) const {
    Segment<Number> ac = {inst->points[ia], inst->points[ic]};
    Segment<Number> bd = {inst->points[ib], inst->points[id]};

    if(ac.cross(bd))
      return true;

    return false;
  }

};
