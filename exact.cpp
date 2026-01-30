//g++ -Wall -o main -std=c++20 -Ofast main.cpp
#include "tools.hpp"
#include "triangulation.hpp"
#include "instance.hpp"
#include "distance.hpp"
#include "center.hpp"
#include "solution.hpp"
#include "satdistance.hpp"
#include "satcenter.hpp"
#include "EvalMaxSAT.h"
#include <algorithm>
#include <cstddef>
#include <csignal>

using Number = long long int;


void printMatrix(const std::vector<std::vector<int>> &distanceMatrix) {
  for(auto &v : distanceMatrix) {
    for(int d : v)
      std::cout << d << ' ';
    std::cout << std::endl;
  }
}

std::vector<Triangulation<Number>> getHeuristicPath(Triangulation<Number> &t1, Triangulation<Number> &t2) {
  DistanceSolver<Number> heuristic1(t1, t2);
  heuristic1.setReverse();
  heuristic1.squeaky();
  DistanceSolver<Number> heuristic2(t2, t1);
  heuristic2.squeaky();
  if(heuristic1.getSolution().size() < heuristic2.getSolution().size())
    return heuristic1.getSolution();
  return heuristic2.getSolution();
}

int calculateDistance(Instance<Number> &inst, Triangulation<Number> &t1, Triangulation<Number> &t2, bool happy) {
    int ub = getHeuristicPath(t1, t2).size() - 1;
    SATDistanceSolver<Number> solver(inst, t1, t2, false, happy);
    int ret = solver.solveDecreasingDistanceSAT(ub-1);
    if(ret == 0)
      return ub;
    if(ret == -1)
      return -1;
    return solver.getLength();
}

std::vector<std::vector<int>> calculateAllDistances(Instance<Number> &inst, bool happy) {
  verbose = 0;
  int t = inst.triangulations.size();
  std::vector<std::vector<int>> distanceMatrix(t, std::vector<int>(t));

  for(int i = 0; i < t-1; i++) {
    for(int j = i+1; j < t; j++) {
      int d = calculateDistance(inst, inst.triangulations[i], inst.triangulations[j], happy);
      distanceMatrix[i][j] = d;
      distanceMatrix[j][i] = d;
      std::cout << d << ' ' << std::flush;
    }
    std::cout << std::endl;
  }

  verbose = 4;

  return distanceMatrix;
}

int calculateLowerBound(const std::vector<std::vector<int>> &distanceMatrix) {
  // Use EvalMaxSAT to solve a matching-like problem

  EvalMaxSAT solver;
  int t = distanceMatrix.size();

  std::map<std::pair<int,int>,int> vars;
  for(int i = 0; i < t; i++) {
    for(int j = 0; j < t; j++) {
      if(i != j) {
        vars[{i,j}] = solver.newVar();
      }
    }
  }

  // Hard clauses: Cannot have two edges leaving or goint to the same vertex
  for(int i = 0; i < t; i++) {
    for(int j = 0; j < t; j++) {
      if(i != j) {
        for(int j2 = j+1; j2 < t; j2++) {
          if(i != j2) {
            solver.addClause({-vars[{i,j}], -vars[{i,j2}]});
            solver.addClause({-vars[{j,i}], -vars[{j2,i}]});
          }
        }
      }
    }
  }

  // Soft clauses
  for(auto [e,var] : vars) {
    int weight = distanceMatrix[e.first][e.second];
    if(weight > 0)
      solver.addClause({var}, weight);
  }

  solver.setTargetComputationTime(10); // Target time
  // solver.setBoundRefTime(.3, 50.0);
  solver.setBoundRefTime(0.0, 0.0);
  solver.disableOptimize();
  // solver.printInfo();

  int ret = solver.solve();

  if(ret == 1 && solver.getCost() < std::numeric_limits<t_weight>::max()) {
    int value = 0;
    // std::vector<std::vector<int>> m(t, std::vector<int>(t));
    for(auto [e,var] : vars) {
      // m[e.first][e.second] = solver.getValue(var);
      if(solver.getValue(var))
        value += distanceMatrix[e.first][e.second];
    }
    // printMatrix(m);

    return (value+1)/2;
  }

  return -1;
}

void generate(int idx, int k, std::vector<int>& current, const std::vector<std::vector<int>> &distanceMatrix, std::vector<std::vector<int>> &result) {
  int n = distanceMatrix.size();
  if (idx == n - 1) {
    // Last position: its value is determined to satisfy the sum = k
    current[idx] = k;
    bool skip = false;
    for(int previous = 0; previous < idx; previous++) {
      if(current[previous] + current[idx] < distanceMatrix[previous][idx]) {
        skip = true;
        break;
      }
    }
    if(!skip) {
      result.push_back(current);
      int tot = 0;
      for(size_t i = 0; i < current.size(); i++) {
        std::cout << current[i];
        tot += current[i];
        if(i != current.size() -1)
          std::cout << " + ";
      }
      std:: cout << " = " << tot << std::endl;
    }
    return;
  }

  for (int i = 0; i <= k; i++) {
    current[idx] = i;
    bool skip = false;
    for(int previous = 0; previous < idx; previous++) {
      if(current[previous] + current[idx] < distanceMatrix[previous][idx]) {
        skip = true;
        break;
      }
    }
    if(!skip)
      generate(idx + 1, k - i, current, distanceMatrix, result);
  }
}


std::vector<std::vector<int>> calculatePossibleLengths(const std::vector<std::vector<int>> &distanceMatrix, int lb) {
  int t = distanceMatrix.size();
  std::vector<int> current(t);
  std::vector<std::vector<int>> result;
  generate(0, lb, current, distanceMatrix, result);

  return result;
}

bool solveExact(Instance<Number> &inst, std::vector<int> &lengths, bool happy) {
  SATCenterSolver<Number> solver(inst, lengths, false, happy);
  if(solver.solve()) {
    std::cout << "SOLVED!" << std::endl;
    solver.getSolution().save(inst.name+(happy?".exacthappy.sol.json":".exact.sol.json"));
    return true;
  }

  std::cout << "No solution!" << std::endl;
  return false;
}

bool solveForValue(Instance<Number> &inst, std::vector<std::vector<int>> &distanceMatrix, int value, bool happy) {
  std::cout << "Possible path lengths: " << std::endl;
  std::vector<std::vector<int>> possibleLengths = calculatePossibleLengths(distanceMatrix, value);
  for(auto &lengths : possibleLengths) {
    if(solveExact(inst, lengths, happy))
      return true;
  }
  return false;
}

void solveExact(Instance<Number> &inst, bool happy) {
  std::vector<std::vector<int>> distanceMatrix = calculateAllDistances(inst, happy);
  int lb = calculateLowerBound(distanceMatrix);
  std::cout << "Distance matrix\n";
  printMatrix(distanceMatrix);
  std::cout << "Lower bound: " << lb << std::endl;

  for(int value = lb; ; value++) {
    {
      std::ofstream f(inst.name+(happy?".lbh":".lb"));
      f << value << std::endl;
    }
    std::cout << "\nSolving for value " << value << std::endl;
    if(solveForValue(inst, distanceMatrix, value, happy))
      break;
  }
}

void solvePair(Instance<Number> &inst, bool happy) {
  Triangulation<Number> &t0 = inst.triangulations.at(0);
  Triangulation<Number> &t1 = inst.triangulations.at(1);

  verbose = 0;
  std::vector<Triangulation<Number>> path = getHeuristicPath(t0, t1);
  verbose = 4;
  int ub = path.size() - 1;

  std::cout << "Heuristic found a path of length " << ub << std::endl;

  SATDistanceSolver<Number> solver(inst, t0, t1, false, happy);
  int ret = solver.solveDecreasingDistanceSAT(ub-1);
  if(ret == 1) {
    std::cout << "SOLVED!" << std::endl;
    path = solver.getReverseSolution();
  }

  Solution<Number> sol(inst);
  sol.setPath(0, path);
  sol.setPath(1, {t1});
  sol.calculateCenter();
  sol.calculateLength();
  sol.save(inst.name+(happy?".exacthappy.sol.json":".exact.sol.json"));
}

int main(int argc, char **argv) {
  std::string fn;
  bool happy = false;

  if(argc == 2) {
    fn = argv[1];
  }
  else if(argc == 3) {
    std::cout << "Assuming happy edges conjecture" << std::endl;
    happy = true;
    fn = argv[2];
  }
  else {
    std::cout << "./exact [happy] instance.json" << std::endl;
    std::cout << "Use happy to assume the happy edges conjecture." << std::endl;
    return 1;
  }

  Instance<Number> inst(fn);
  inst.meta["command"] = concatenate(argc, argv);

  if(inst.triangulations.size() == 2) {
    std::cout << "Only 2 input triangulations, solving path directly" << std::endl;
    solvePair(inst,happy);
  }
  else {
    solveExact(inst,happy);
  }

  return 0;
}
