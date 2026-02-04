#pragma once
#include "segment.hpp"
#include "tools.hpp"
#include "triangulation.hpp"
#include "graph.hpp"
#include "independent.hpp"
#include "instance.hpp"
#include "solution.hpp"
#include <sstream>
#include <string>
#include <unistd.h>
#include "cadicalinterface.h"
#include "EvalMaxSAT.h"
#include "crosslb.hpp"

template <class Number>
class SATCenterSolver {
  using Clause = std::vector<int>;
  using Path = std::vector<Triangulation<Number>>;
  Instance<Number> *inst;
  std::vector<Triangulation<Number>> sources;
  std::vector<int> lengths;
  Solution<Number> solution;
  std::map<std::tuple<IEdge,int,int>, int> edgeVars;
  std::map<std::tuple<IEdge,IEdge,int,int>, int> flipVars;
  std::vector<Clause> clauses;
  std::unordered_set<IEdge> happyEdges;
  bool inexact = true, happyConjecture = true;
  int nvars;
  double centerDist = 0, pathDist = 0;

public:
  std::string instName;
  SATCenterSolver(Instance<Number> &_inst, const std::vector<int> &_lengths, double _centerDist, const Solution<Number> &_solution, double _pathDist, bool _inexact, bool _happyConjecture)
    : inst(&_inst), sources(_inst.triangulations), lengths(_lengths), solution(_solution), centerDist(_centerDist), pathDist(_pathDist), inexact(_inexact), happyConjecture(_happyConjecture) {
    if(inexact)
      CrossLB::load_lbfile("lb.txt");
  }

  SATCenterSolver(Instance<Number> &_inst, const std::vector<int> &_lengths, bool _inexact, bool _happyConjecture)
    : inst(&_inst), sources(_inst.triangulations), lengths(_lengths), solution(_inst), inexact(_inexact), happyConjecture(_happyConjecture) {
    if(inexact)
      CrossLB::load_lbfile("lb.txt");
  }

  Solution<Number> &getSolution() {
    return solution;
  }

  int solve(double timeout = 0) {
    std::cout << "Solving for lengths: ";
    for(int l : lengths)
      std::cout << l << " ";
    std::cout << std::endl;
    try {
      buildClauses();
    }
    catch(std::out_of_range &e) {
      return false;
    }

    int satisfiable = solveSAT(timeout);

    return satisfiable;
  }

  int solveSAT(double timeout = 0) {
    Solver_cadical solver;

    flipVars = {}; // To free memory
    for(Clause &c :clauses) {
      solver.addClause(c);
    }
    clauses = {}; // To free memory
    int ret = timeout > 0 ? solver.solveWithTimeout(timeout) : solver.solve();
    if(ret == 1) {
      std::vector<bool> satSol = solver.getSolution();
      std::vector<std::vector<std::vector<IEdge>>> edgesPath(lengths.size());
      for(int i = 0; i < int(lengths.size()); i++) {
        edgesPath[i] = std::vector<std::vector<IEdge>>(lengths[i]+1);
      }

      for(auto [eti,var] : edgeVars) {
        if(satSol.at(var))
          edgesPath.at(std::get<2>(eti)).at(std::get<1>(eti)).push_back(std::get<IEdge>(eti));
      }

      for(int i = 0; i < int(lengths.size()); i++) {
        std::vector<Triangulation<Number>> path;
        for(int j = edgesPath[i].size() - 1; j >= 0; j--) {
          auto &edges = edgesPath[i][j];
          path.push_back(Triangulation<Number>(&inst->points, edges));
        }
        if(i == 0)
          solution.setCenter(path[0]);
        solution.improvePath(i, path);
      }
    }
    return ret;
  }

  const std::map<std::tuple<IEdge,int,int>, int> &getEdgeVars() const {
    return edgeVars;
  }

  void updateEdgeVars(const std::map<std::tuple<IEdge,int,int>, int> &evars, int ired) {
    edgeVars.clear();

    // Copy edges from non-reduced path
    for(auto [eti,v] : evars) {
      IEdge e = std::get<0>(eti);
      int t = std::get<1>(eti);
      int i = std::get<2>(eti);
      if(i != ired) {
        edgeVars[{e,t,i}] = v;
      }
    }

    // Trim unused variables
    nvars = 0;
    std::vector<std::pair<int,std::tuple<IEdge,int,int>>> varVector;
    for(auto [eti,v] : edgeVars)
      varVector.push_back({v,eti});
    std::stable_sort(varVector.begin(),varVector.end());
    for(int previous = -1; auto [oldVar,eti] : varVector) {
      if(previous != oldVar) {
        nvars++;
        previous = oldVar;
      }
      edgeVars.at(eti) = nvars;
    }

    // Add edges for trimmed path
    for(auto [eti,v] : evars) {
      IEdge e = std::get<0>(eti);
      int t = std::get<1>(eti);
      int i = std::get<2>(eti);
      if(i == ired) {
        int ic = ired != 0 ? 0 : 1;
        if(t >= lengths.at(ired)) {
          t = std::min(t,lengths[i]);
          if(edgeVars.contains({e,lengths[ic],ic}))
            edgeVars[{e,t,i}] = edgeVars.at({e,lengths[ic],ic});
        }
        else {
          if(!edgeVars.contains({e,t,i}))
            edgeVars[{e,t,i}] = ++nvars;
          if(t >= 2) {
            t--;
            if(!edgeVars.contains({e,t,i}))
              edgeVars[{e,t,i}] = ++nvars;
          }
        }

      }
    }


    if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Updated edge variables: " << nvars << ", " << edgeVars.size() << std::endl;
  }

  void buildEdgeVariables() {
    if(!edgeVars.empty())
      return;

    nvars = 0;
    edgeVars.clear();
    const int n = inst->points.size();

    if(inexact || happyConjecture) {
      happyEdges.clear();
      for(IEdge e : sources[0].getEdges()) {
        bool isHappy = true;
        for(int i = 1; i < int(sources.size()); i++) {
          if(!sources[i].contains(e)) {
            isHappy = false;
            break;
          }
        }
        if(isHappy)
          happyEdges.insert(e);
      }

      if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Happy edges: " << happyEdges.size() << std::endl;
    }

    std::unordered_map<IEdge,int> targetEdges;

    for(short u = 0; u < n-1; u++) {
      for(short v = u+1; v < n; v++) {
        if(!emptySegment(u,v))
          continue;
        if((inexact  || happyConjecture) && crossHappy({u,v}))
          continue;

        std::vector<IEdge> xoc;
        if(centerDist > 0)
          xoc = solution.getCenter().getCrossings({u,v});

        for(int i = 0; i < int(lengths.size()); i++) {
          auto xs = sources.at(i).getCrossings({u,v}, 1 << lengths[i]);
          for(int t = 0; t <= lengths[i]; t++) {
            if(CrossLB::simpleCrossLB(xs) <= t
               && (!inexact || CrossLB::fancyCrossLB(xs) <= t)) {
              if(centerDist == 0
                 || (xoc.size() < pow(2, lengths[i]-t+centerDist) && (!inexact || CrossLB::fancyCrossLB(xoc) <= lengths[i]-t+ceil(centerDist)))) {
                if(pathDist != 0 && (int) solution.getPath(i).size() - 1 > lengths[i] - t) {
                  auto xx = solution.getPath(i).at(lengths[i]-t).getCrossings({u,v}, ceil(pow(2,pathDist)));
                  if(xx.size() >= pow(2,pathDist))
                    continue;
                }

                if(t == lengths[i]) { // Same target for all paths
                  targetEdges[{u,v}]++;
                }
                else if(!edgeVars.contains({{u,v},t,i})) {
                  edgeVars[{{u,v},t,i}] = ++nvars;
                }
              }
            }
          }
        }
      }
    }

    for(auto [e, val] : targetEdges) {
      if(val == lengths.size()) { // Edge in all paths
        edgeVars[{e,lengths[0],0}] = ++nvars;
        for(int i = 1; i < int(lengths.size()); i++) {
          edgeVars[{e,lengths[i],i}] = nvars; // Link to same variable
        }
      }
    }

    if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Edge variables: " << nvars << ", " << edgeVars.size() << std::endl;
  }

  void buildFlipVariables() {
    flipVars.clear();
    const int n = inst->points.size();

    std::vector<std::unordered_set<int>> flipableEdges(n);

    for(auto [eti,_] : edgeVars) {
      IEdge e = std::get<IEdge>(eti);
      flipableEdges[e[0]].insert(e[1]);
    }

    // Variables for flips
    for(short u = 0; u < n; u++) {
      for(short v :flipableEdges[u]) {
        for(short u2 = 0; u2 < n; u2++) {
          if(u2 == u || u2 == v || flipableEdges[u2].empty() || !inst->emptyTriangle(u,v,u2))
            continue;
          for(short v2 : flipableEdges[u2]) {
            if(v2 != u && v2 != v && convexQuadrilateral(u,u2,v,v2) && inst->emptyTriangle(u,v,v2)) {
              IEdge uu2 = u<u2 ? IEdge{u,u2} : IEdge{u2,u};
              IEdge vu2 = v<u2 ? IEdge{v,u2} : IEdge{u2,v};
              IEdge uv2 = u<v2 ? IEdge{u,v2} : IEdge{v2,u};
              IEdge vv2 = v<v2 ? IEdge{v,v2} : IEdge{v2,v};
              for(int i = 0; i < int(lengths.size()); i++) {
                for(int t = 0; t < lengths[i]; t++) {
                  if(edgeVars.contains({{u,v},t,i}) &&
                      edgeVars.contains({uu2,t,i}) &&
                      edgeVars.contains({vu2,t,i}) &&
                      edgeVars.contains({uv2,t,i}) &&
                      edgeVars.contains({vv2,t,i}) &&
                      edgeVars.contains({{u2,v2},t+1,i}) &&
                      edgeVars.contains({uu2,t+1,i}) &&
                      edgeVars.contains({vu2,t+1,i}) &&
                      edgeVars.contains({uv2,t+1,i}) &&
                      edgeVars.contains({vv2,t+1,i})) {
                    if(!flipVars.contains({{u,v},{u2,v2},t,i}))
                      flipVars[{{u,v},{u2,v2},t,i}] = ++nvars;
                  }
                }
              }
            }
          }
        }
      }
    }

    if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Flip variables: " << flipVars.size() << std::endl;
  }


  void buildClauses() {
    buildEdgeVariables();
    buildFlipVariables();
    clauses.clear();

    if(inexact || happyConjecture) {
      for(IEdge e : happyEdges)
        for(int i = 0; i < int(lengths.size()); i++)
          for(int t = 1; t < lengths[i] + (i==0); t++)
            clauses.push_back({edgeVars.at({e,t,i})});
    }

    // Clauses for start and target edges
    for(int i = 0; i < int(lengths.size()); i++) {
      for(IEdge e :sources[i].getEdges())
        clauses.push_back({edgeVars.at({e,0,i})});
    }

    // Clauses for flips
    for(auto [f,var] : flipVars) {
      short u = std::get<0>(f)[0];
      short v = std::get<0>(f)[1];
      IEdge uv = std::get<0>(f);
      short u2 = std::get<1>(f)[0];
      short v2 = std::get<1>(f)[1];
      IEdge u2v2 = std::get<1>(f);
      int t = std::get<2>(f);
      int i = std::get<3>(f);
      IEdge uu2 = u<u2 ? IEdge{u,u2} : IEdge{u2,u};
      IEdge vu2 = v<u2 ? IEdge{v,u2} : IEdge{u2,v};
      IEdge uv2 = u<v2 ? IEdge{u,v2} : IEdge{v2,u};
      IEdge vv2 = v<v2 ? IEdge{v,v2} : IEdge{v2,v};

      // All edges except u2v2 must be present at time t
      clauses.push_back({ -var, edgeVars.at({uv,t,i}) });
      clauses.push_back({ -var, edgeVars.at({uu2,t,i}) });
      clauses.push_back({ -var, edgeVars.at({vu2,t,i}) });
      clauses.push_back({ -var, edgeVars.at({uv2,t,i}) });
      clauses.push_back({ -var, edgeVars.at({vv2,t,i}) });

      // All edges except uv must be present at time t+1
      clauses.push_back({ -var, edgeVars.at({u2v2,t+1,i}) });
      clauses.push_back({ -var, edgeVars.at({uu2,t+1,i}) });
      clauses.push_back({ -var, edgeVars.at({vu2,t+1,i}) });
      clauses.push_back({ -var, edgeVars.at({uv2,t+1,i}) });
      clauses.push_back({ -var, edgeVars.at({vv2,t+1,i}) });

      // Edges uv is flipped into u2v2 at time t+1
      if(edgeVars.contains({uv,t+1,i}))
        clauses.push_back({ -var, -edgeVars.at({uv,t+1,i}) });
      // clauses.push_back({ -var, edgeVars.at({u2v2,t+1,i}) });
    }

    // Data structures
    std::vector<std::tuple<int, int, IEdge, IEdge, int>> indexedFlipVars;
    std::vector<std::tuple<int, int, IEdge, IEdge, int>> indexedBackFlipVars;
    for(auto [f,var] : flipVars) {
      indexedFlipVars.push_back({std::get<2>(f), std::get<3>(f), std::get<0>(f), std::get<1>(f),var});
      indexedBackFlipVars.push_back({std::get<2>(f), std::get<3>(f), std::get<1>(f), std::get<0>(f),var});
    }
    std::sort(indexedFlipVars.begin(), indexedFlipVars.end());
    std::sort(indexedBackFlipVars.begin(), indexedBackFlipVars.end());


    // Edges cannot disappear without a flip
    for(auto [e,var] : edgeVars) {
      IEdge uv = std::get<0>(e);
      int t = std::get<1>(e);
      int i = std::get<2>(e);
      if(t == lengths[i])
        continue;
      Clause c = {-var};

      if(edgeVars.contains({uv,t+1,i}))
        c.push_back(edgeVars.at({uv,t+1,i}));

      // Slow version
      // for(auto [f,varf] : flipVars) {
      //   if(std::get<2>(f) == t && std::get<3>(f) == i && std::get<0>(f) == uv) {
      //     c.push_back(varf);
      //   }
      // }

      auto it = std::upper_bound(indexedFlipVars.begin(), indexedFlipVars.end(), std::make_tuple(t,i,uv,IEdge{-1,-1},-1));
      for(; it != indexedFlipVars.end() &&
            std::get<0>(*it) == t &&
            std::get<1>(*it) == i &&
            std::get<2>(*it) == uv;    ++it)
        c.push_back(std::get<4>(*it));

      clauses.push_back(c);
    }


    // Edges cannot appear without a flip
    for(auto [e,var] : edgeVars) {
      IEdge u2v2 = std::get<0>(e);
      int tplus1 = std::get<1>(e);
      int i = std::get<2>(e);

      if(tplus1 == 0)
        continue;
      Clause c = {-var};

      if(edgeVars.contains({u2v2,tplus1-1,i}))
        c.push_back(edgeVars.at({u2v2,tplus1-1,i}));

      // Slow version:
      // for(auto [f,varf] : flipVars) {
      //   if(std::get<2>(f) == tplus1-1 && std::get<3>(f) == i && std::get<1>(f) == u2v2) {
      //     c.push_back(varf);
      //   }
      // }

      auto it = std::upper_bound(indexedBackFlipVars.begin(), indexedBackFlipVars.end(), std::make_tuple(tplus1-1,i,u2v2,IEdge{-1,-1},-1));
      for(; it != indexedBackFlipVars.end() &&
            std::get<0>(*it) == tplus1-1 &&
            std::get<1>(*it) == i &&
            std::get<2>(*it) == u2v2; ++it)
        c.push_back(std::get<4>(*it));


      clauses.push_back(c);
    }

    if(verbose >= 3) std::cout << int(elapsed()) << "s: " << "Clauses: " << clauses.size() << std::endl;
  }



protected:
  bool crossHappy(IEdge iuv) {
    auto &points = inst->points;

    if(happyEdges.size() > 10*sqrt(points.size())) {
      return sources[0].cross(iuv, happyEdges);
    }

    Segment<Number> uv = {points[iuv[0]], points[iuv[1]]};
    for(IEdge ie :happyEdges) {
      Segment<Number> e = {points[ie[0]], points[ie[1]]};
      if(uv.cross(e))
        return true;
    }

    return false;
  }

  bool emptySegment(int ia, int ib) const {
    auto &points = inst->points;
    Segment<Number> t = {points[ia], points[ib]};
    for(Point<Number> p : points) {
      if(t.containsInside(p))
        return false;
    }

    return true;
  }

  bool convexQuadrilateral(int ia, int ib, int ic, int id) const {
    auto &points = inst->points;

    Segment<Number> ac = {points[ia], points[ic]};
    Segment<Number> bd = {points[ib], points[id]};

    if(ac.cross(bd))
      return true;

    return false;
  }
};
