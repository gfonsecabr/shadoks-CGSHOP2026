//g++ -Wall -o happy -std=c++20 -Ofast happy.cpp
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

using Number = long long int;

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

bool testConjecture(Instance<Number> &inst) {
  verbose = 0;
  int t = inst.triangulations.size();
  int iter = 0;

  for(int i = 0; i < t-1; i++) {
    for(int j = i+1; j < t; j++) {
      iter++;
      std::cout << "From " << i << " to " << j << ", " << std::flush;
      int exactDistance = calculateDistance(inst, inst.triangulations[i], inst.triangulations[j], false);
      std::cout << " exact = " << exactDistance << std::flush;
      int happyDistance = calculateDistance(inst, inst.triangulations[i], inst.triangulations[j], true);
      std::cout << " happy = " << happyDistance << std::endl;
      if(happyDistance != exactDistance) {
        std::cout << "The conjecture has been falsified!\n|";
        return false;
      }
    }
  }

  std::cout << "Tested: " << iter << std::endl;

  return true;
}


int main(int argc, char **argv) {
  std::string fn;
  bool happy = false;

  if(argc == 2) {
    fn = argv[1];
  }
  else {
    std::cout << "./happy instance.json" << std::endl;
    return 1;
  }

  Instance<Number> inst(fn);

  if(!testConjecture(inst))
    return 2;

  return 0;
}
