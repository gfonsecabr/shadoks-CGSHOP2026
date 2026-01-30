#include <string>
inline struct {
  std::string fn;
  std::string outfn = "";
  std::string instancePath = "../instances/";
  double centerDistance = 2;
  int trimRadius = 0;
  int rebuildRadius = 0;
  bool rebuildOnly = false;
  bool maxsat = true;
  double pathDistance = 0;
  double satTimeout = 36000;
  bool inexact = true, happy = true;
  bool repeat = false;
} PARM;


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
#include "cxxopts.hpp"
#include "boost/pfr.hpp"

using Number = long long int;

int parse(int argc, char* argv[]) {
  try {
    std::unique_ptr<cxxopts::Options> allocated(new cxxopts::Options(argv[0], "Solution improver for CG:SHOP 2026"));
    auto& options = *allocated;
    options
      .positional_help("solution.json")
      .show_positional_help();

    options
      .set_width(75)
      .set_tab_expansion()
      // .allow_unrecognised_options()
      .add_options()
      ("s,solution",
        "Name of the input solution json file (positional argument)", cxxopts::value<std::string>(PARM.fn))
      ("o,output", "Output file name (generated automatically if empty)", cxxopts::value<std::string>(PARM.outfn)->default_value(PARM.outfn))
      ("i,instance", "Name of the instance file name or path", cxxopts::value<std::string>(PARM.instancePath)->default_value(PARM.instancePath))
      ("c,center", "Force distance to previous center, i.e. less than 2^c crossings, 0 for same center, -1 for disabled", cxxopts::value<double>(PARM.centerDistance)->default_value(std::to_string(PARM.centerDistance)))
      ("t,trim", "Trim radius, 0 for no trimming", cxxopts::value<int>(PARM.trimRadius)->default_value(std::to_string(PARM.trimRadius)))
      ("r,rebuild", "Rebuild radius, 0 for no rebuilding", cxxopts::value<int>(PARM.rebuildRadius)->default_value(std::to_string(PARM.rebuildRadius)))
      ("R,rebuildonly", "Exit after rebuild", cxxopts::value<bool>(PARM.rebuildOnly)->default_value(std::to_string(PARM.rebuildOnly)))
      ("maxsat", "Use maxsat to rebuild", cxxopts::value<bool>(PARM.maxsat)->default_value(std::to_string(PARM.maxsat)))
      ("d,distance", "Distance from previous solution path, i.e. less than 2^d crossings, 0 for disabled", cxxopts::value<double>(PARM.pathDistance)->default_value(std::to_string(PARM.pathDistance)))
      ("I,inexact", "Use inexact crossing lower bound (requires lb.txt)", cxxopts::value<bool>(PARM.inexact)->default_value(std::to_string(PARM.inexact)))
      ("H,happy", "Assume happy edges conjecture", cxxopts::value<bool>(PARM.happy)->default_value(std::to_string(PARM.happy)))
      ("l,loop", "Loop while there is improvement", cxxopts::value<bool>(PARM.repeat)->default_value(std::to_string(PARM.repeat)))
      ("timeout", "Stop SAT solver after t seconds without improvement", cxxopts::value<double>(PARM.satTimeout)->default_value(std::to_string(PARM.satTimeout)))
      ("h,help", "Print help")
    ;

    options.parse_positional("solution");

    auto result = options.parse(argc, argv);

    if (!result.count("solution") || result.count("help")) {
      std::cout << options.help({"", "Group"}) << std::endl;
      return true;
    }
  }
  catch (const cxxopts::exceptions::exception& e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    return false;
  }

  return true;
}

void rebuild(Solution<Number> &sol, int delta) {
  auto &triangulations = sol.getInstance().triangulations;

  for(unsigned i = 0; i < triangulations.size(); i++) {
    SATDistanceSolver<Number> dSolver(sol.getInstance(), triangulations[i], sol.getCenter(), PARM.inexact, PARM.happy);

    int oldLength = sol.getPath(i).size()-1;
    if(verbose >= 3) std::cout << "Rebuilding path " << i << " of length " << oldLength << std::endl;
    if(PARM.maxsat) {
      if(dSolver.solveExactDistanceIterativeMaxSAT(oldLength,delta)) {
        sol.setPath(i, dSolver.getReverseSolution());
      }
    }
    else {
      if(dSolver.solveDecreasingDistanceSAT(oldLength-1)) {
        sol.improvePath(i, dSolver.getReverseSolution());
      }
    }
  }
  if(verbose >= 2) std::cout << "= " << sol.getLength() << std::endl;
}

void printLenths(const std::vector<int> &lengths) {
  int tot = 0;
  for(size_t i = 0; i < lengths.size(); i++) {
    std::cout << lengths[i];
    tot += lengths[i];
    if(i != lengths.size() -1)
      std::cout << " + ";
  }
  std:: cout << " = " << tot << std::endl;
}

bool reduceLengths(Instance<Number> &inst, std::vector<int> &lengths, Solution<Number> &sol) {
  std::map<std::tuple<IEdge,int,int>, int> edgeVars;
  std::unordered_set<int> failed;
  std::uniform_int_distribution<int> dist(0, lengths.size()-1);
  bool improved = false;

  while(failed.size() != lengths.size()) {
    int i = dist(rgen);
    if(failed.contains(i))
      continue;

    if(edgeVars.empty()) {
      SATCenterSolver<Number> tmp(inst, lengths, PARM.centerDistance, sol, PARM.pathDistance, PARM.inexact, PARM.happy);
      std::cout << int(elapsed()) << "s: " << "Building template edge variables" << std::endl;
      tmp.buildEdgeVariables();
      edgeVars = tmp.getEdgeVars();
    }

    auto newLengths = lengths;
    newLengths[i]--;
    std::cout << lengths.size() - failed.size() << " candidates. ";
    std::cout << "Reducing term " << i << "/" << lengths.size() << " to " << newLengths[i] << std::endl;
    // printLenths(newLengths);

    SATCenterSolver<Number> solver(inst, newLengths, PARM.centerDistance, sol, PARM.pathDistance, PARM.inexact, PARM.happy);
    if(!edgeVars.empty())
      solver.updateEdgeVars(edgeVars, i);

    int status = solver.solve(PARM.satTimeout);

    if(status == 0) {
      std::cout << int(elapsed()) << "s: " << "No solution!" << std::endl;
      failed.insert(i);
    }
    else if(status < 0) {
      std::cout << int(elapsed()) << "s: " << "Timeout!" << std::endl;
      failed.insert(i);
    }
    else {
      std::cout << int(elapsed()) << "s: " << "SOLVED!" << std::endl;
      improved = true;
      if(PARM.trimRadius == 0) {
        sol = solver.getSolution();
      }
      else {
        sol.improveFromSub(solver.getSolution());
      }
      sol.save(PARM.outfn);

      lengths = newLengths;
      edgeVars.clear();
    }
  }

  return improved;
}


int main(int argc, char **argv) {
  if(!parse(argc, argv))
    return 1;
  std::stringstream parmio;
  parmio << boost::pfr::io(PARM);
  std::cout << "Parameters: " << parmio.str() << std::endl;

  verbose = 4;

  std::string uid;
  auto flips = readSolutionFlips<Number>(PARM.fn, uid);
  if(PARM.instancePath.empty() || PARM.instancePath.back() == '/')
    PARM.instancePath += uid + ".json";

  Instance<Number> inst(PARM.instancePath);
  inst.meta["command"] = concatenate(argc, argv);
  Solution<Number> sol(inst, flips);

  bool improved = false;
  do {
    std::vector<int> lengths = sol.getLengths();
    std::cout << "Original lengths: ";
    printLenths(lengths);
    if(PARM.trimRadius != 0 || PARM.rebuildRadius > 0) {
      if(PARM.rebuildRadius > 0) {
        rebuild(sol, PARM.rebuildRadius);
        sol.save(PARM.outfn, ".rebuilt" + std::to_string(PARM.rebuildRadius));

        if(PARM.rebuildOnly) {
          std::cout << "Only rebuild selected\n";
          exit(0);
        }
        verbose = 4;
      }
      if(PARM.trimRadius != 0) {
        Instance<Number> sub = sol.subInstance(PARM.trimRadius);

        for(size_t i = 0; i < lengths.size(); i++) {
          lengths[i] = std::min(PARM.trimRadius, (int)sol.getPath(i).size()-1);
        }
        std::cout << "Trimmed lengths: ";
        printLenths(lengths);
        improved = reduceLengths(sub, lengths, sol);
      }
      else {
        improved = reduceLengths(inst, lengths, sol);
      }
    }
    else {
      improved = reduceLengths(inst, lengths, sol);
    }
  } while(PARM.repeat && improved);

  return 0;
}
