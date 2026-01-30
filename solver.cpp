#include <string>
inline struct {
  std::string fn;
  int sat = 0;
  int timeout = 7200;
  double best = .75;
  int littlePoolSize = 150;
  int largePoolSize = 200;
  double minExp = .5;
  double maxExp = 4.0;
  double incExp = .5;
  int children = 16;
  int satonsave = 0;
  bool inexact = true, happy = true;
} PARM;

#include "tools.hpp"
#include "triangulation.hpp"
#include "instance.hpp"
#include "distance.hpp"
#include "center.hpp"
#include "solution.hpp"
#include "delaunay.hpp"
#include "satdistance.hpp"
#include <algorithm>
#include <cstddef>
#include <csignal>
#include "cxxopts.hpp"
#include "boost/pfr.hpp"

using Number = long long int;

int parse(int argc, char* argv[]) {
  try {
    std::unique_ptr<cxxopts::Options> allocated(new cxxopts::Options(argv[0], "Initial solution solver for CG:SHOP 2026"));
    auto& options = *allocated;
    options
      .positional_help("instance.json")
      .show_positional_help();

    options
      .set_width(75)
      .set_tab_expansion()
      // .allow_unrecognised_options()
      .add_options()
      ("i,instance",
        "Name of the instance json file (positional argument)", cxxopts::value<std::string>(PARM.fn))
      ("s,sat", "When to improve distances with MaxSAT: -1 for always, 0 for never, 1,2,3...", cxxopts::value<int>(PARM.sat)->default_value(std::to_string(PARM.sat)))
      ("t,timeout", "Stop after t seconds without improvement", cxxopts::value<int>(PARM.timeout)->default_value(std::to_string(PARM.timeout)))
      ("p,prob", "Probability of improving best solution instead of first solution in FIFO", cxxopts::value<double>(PARM.best)->default_value(std::to_string(PARM.best)))
      ("littlepool", "Target pool size after pruning", cxxopts::value<int>(PARM.littlePoolSize)->default_value(std::to_string(PARM.littlePoolSize)))
      ("largepool", "Size at which to prune pool", cxxopts::value<int>(PARM.largePoolSize)->default_value(std::to_string(PARM.largePoolSize)))
      ("minexp", "Minimum exponent for initial centers", cxxopts::value<double>(PARM.minExp)->default_value(std::to_string(PARM.minExp)))
      ("maxexp", "Maximum exponent for initial centers", cxxopts::value<double>(PARM.maxExp)->default_value(std::to_string(PARM.maxExp)))
      ("incexp", "Increment of exponent for initial centers", cxxopts::value<double>(PARM.incExp)->default_value(std::to_string(PARM.incExp)))
      ("c,children", "Maximum number of child nodes to explore", cxxopts::value<int>(PARM.children)->default_value(std::to_string(PARM.children)))
      ("x,satonsave", "Run exact sat solver when saving after n steps. 0 for disabled", cxxopts::value<int>(PARM.satonsave)->default_value(std::to_string(PARM.satonsave)))
      ("I,inexact", "Use inexact crossing lower bound (requires lb.txt)", cxxopts::value<bool>(PARM.inexact)->default_value(std::to_string(PARM.inexact)))
      ("H,happy", "Assume happy edges conjecture", cxxopts::value<bool>(PARM.happy)->default_value(std::to_string(PARM.happy)))
      ("h,help", "Print help")
    ;

    options.parse_positional("instance");

    auto result = options.parse(argc, argv);


    if (!result.count("instance") || result.count("help")) {
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

using FlipList = std::vector<std::vector<std::vector<IEdge>>>;

void improvePaths(Solution<Number> &sol) {
  auto &triangulations = sol.getInstance().triangulations;

  for(unsigned i = 0; i < triangulations.size(); i++) {
    SATDistanceSolver<Number> dSolver(sol.getInstance(), triangulations[i], sol.getCenter(), PARM.inexact, PARM.happy);

    int oldLength = sol.getPath(i).size()-1;
    if(dSolver.solveDecreasingDistanceMaxSAT(oldLength)) {
      if(sol.improvePath(i, dSolver.getReverseSolution())) {
        if(verbose >= 2) std::cout << i << ": " << oldLength << " -> " << dSolver.getLength() << " " << std::flush;
      }
    }
  }
  if(verbose >= 2) std::cout << "= " << sol.getLength() << std::endl;
}

void improvePathsSAT(Solution<Number> &sol) {
  auto &triangulations = sol.getInstance().triangulations;

  for(unsigned i = 0; i < triangulations.size(); i++) {
    SATDistanceSolver<Number> dSolver(sol.getInstance(), triangulations[i], sol.getCenter(), PARM.inexact, PARM.happy);

    int oldLength = sol.getPath(i).size()-1;
    if(dSolver.solveDecreasingDistanceSAT(oldLength-1)) {
      if(sol.improvePath(i, dSolver.getReverseSolution())) {
        if(verbose >= 2) std::cout << i << ": " << oldLength << " -> " << dSolver.getLength() << " " << std::flush;
      }
    }
  }
  if(verbose >= 2) std::cout << "= " << sol.getLength() << std::endl;
}

void setPaths(Solution<Number> &sol) {
  std::vector<std::vector<Triangulation<Number>>> paths;
  auto &triangulations = sol.getInstance().triangulations;

  for(unsigned i = 0; i < triangulations.size(); i++) {
    if(verbose >= 4) std::cout << "To " << i << ":";

    if(verbose >= 2 && !sol.getPath(i).empty()) std::cout << sol.getPath(i).size()-1 << '_';

    int goal = int(sol.getPath(i).size()) - 3;

    {
      DistanceSolver<Number> dSolver(triangulations[i], sol.getCenter());
      dSolver.setReverse();
      dSolver.squeaky(16, goal);
      sol.improvePath(i, dSolver.getSolution());
    }

    {
      // Find distance backwards
      if(verbose >= 4) std::cout << " |";
      else if (verbose >= 2) std::cout << "_";

      DistanceSolver<Number> dSolver(sol.getCenter(), triangulations[i]);
      dSolver.squeaky(16, goal);
      sol.improvePath(i, dSolver.getSolution());
    }

    if(verbose >= 4) std::cout << " = " << sol.getPath(i).size()-1 << std::endl;
    // if(verbose >= 2 && verbose < 4 ) std::cout << "_" << sol.getPath(i).size()-1 << std::flush;
    if(verbose >= 2 && verbose < 4 && i != triangulations.size()-1) std::cout << " + " << std::flush;
  }

  if(verbose >= 4) std::cout << "Length: " << sol.getLength() << std::endl << std::endl;
  if(verbose >= 2 && verbose < 4) std::cout << " = " << sol.getLength() << std::endl;
}

void pruneSols(std::deque<std::tuple<int,std::vector<int>,FlipList>> &goodSols, int limit) {
  if(verbose >= 2) std::cout << "Pruning\n";
  std::sort(goodSols.begin(), goodSols.end());
  auto it = goodSols.begin() + limit;
  if(it < goodSols.end())
    goodSols.erase(it, goodSols.end());
}

FlipList extractFirst(std::deque<std::tuple<int,std::vector<int>,FlipList>> &goodSols) {
  FlipList ret = std::get<2>(goodSols.front());
  goodSols.pop_front();
  return ret;
}

FlipList getBest(const std::deque<std::tuple<int,std::vector<int>,FlipList>> &goodSols) {
  auto it = std::min_element(goodSols.begin(), goodSols.end());
  FlipList ret = std::get<2>(*it);
  return ret;
}

FlipList extractBest(std::deque<std::tuple<int,std::vector<int>,FlipList>> &goodSols) {
  auto it = std::min_element(goodSols.begin(), goodSols.end());
  FlipList ret = std::get<2>(*it);
  goodSols.erase(it);
  return ret;
}

void improveCenter(Solution<Number> &sol, std::deque<std::tuple<int,std::vector<int>,FlipList>> goodSols, double timeOut, int runSAT) {
  double lastImprovement = elapsed();
  int bestSATLength = -1;
  auto &inst = sol.getInstance();
  int t = inst.triangulations.size();
  if(goodSols.empty())
    goodSols.push_back(std::make_tuple(sol.getLength(), sol.sortedLastFlips(), sol.getFlipVectors()));
  std::sort(goodSols.begin(), goodSols.end());
  std::unordered_set<size_t> visitedHashs;
  visitedHashs.max_load_factor(256); // To save memory on very long runs

  for(int step = 0; !goodSols.empty() && elapsed()-lastImprovement < timeOut; step++) {
    // for(auto &s :goodSols)
    //   std::cout << s.getLength() << ' ';
    // std::cout << std::endl;

    if(goodSols.size() > PARM.largePoolSize) {
      pruneSols(goodSols, PARM.littlePoolSize);
    }

    static std::bernoulli_distribution distrib(PARM.best);
    bool randBool = distrib(rgen);
    auto curSol = randBool ? Solution<Number>(inst,extractBest(goodSols)) : Solution<Number>(inst,extractFirst(goodSols));

    if(visitedHashs.insert(curSol.getCenter().hash()).second == false) {
      // Center already visited (or one with the same hash)
      step--;
      continue;
    }

    if(runSAT < 0 || curSol.getLength() <= sol.getLength() + runSAT - 1) {
      if(verbose >= 2) std::cout << "Running SAT for length " << curSol.getLength() << (randBool ? "" : "*") << ": " << std::flush;
      improvePaths(curSol);
      if(sol.improve(curSol))
        lastImprovement = elapsed();
      if(verbose >= 2) std::cout << std::endl;
    }

    auto flipCounts = curSol.decoratedSortedLastFlips();

    for(int ri = 0; ri < t && ri < PARM.children; ri++) {
      int r = flipCounts[ri].second;

      if(curSol.getPath(r).size() <= 1)
        continue;
      if(verbose >= 4) std::cout << step << " (" << goodSols.size() << "): Moving towards " << r << " from [";
      else if(verbose >= 2) {
        int flipCount = -1;
        if(curSol.getPath(r).size() >= 2)
          flipCount = (curSol.getPath(r)[0] - curSol.getPath(r)[1]).size();
        std::cout << step << ": To " << r << "(" << flipCount << ") from [";
      }
      for(int j = 0; j < t; j++)
        if(verbose >= 2) std::cout << curSol.getPath(j).size()-1 << (j != t-1 ? "+":"=");
      if(verbose >= 3) std::cout << curSol.getLength() << "]\n";
      else if(verbose >= 2) std::cout << curSol.getLength() << "]: ";
      Solution<Number> newSol = curSol;
      newSol.stepTo(r);
      setPaths(newSol);
      if(sol.improve(newSol)) {
        lastImprovement = elapsed();
        if(PARM.satonsave && step > PARM.satonsave) {
          auto tmpSol = newSol;
          if(verbose >= 2) std::cout << "Running SAT on save for length " << tmpSol.getLength() << (randBool ? "" : "*") << ": " << std::flush;
          improvePathsSAT(tmpSol);
          if(bestSATLength == -1 || tmpSol.getLength() < bestSATLength) {
            tmpSol.save("",".sat");
            bestSATLength = tmpSol.getLength();
          }
        }
      }
      goodSols.push_back(std::make_tuple(newSol.getLength(), newSol.sortedLastFlips(), newSol.getFlipVectors()));
    }
  }
}

std::deque<std::tuple<int,std::vector<int>,FlipList>> manySolutions(Instance<Number> &inst) {
  std::deque<std::tuple<int,std::vector<int>,FlipList>> goodSols;

  DelaunaySolver<Number> delSolver(inst);

  // delSolver.run();
  delSolver.solveTarget(0);

  Solution<Number> delSol = delSolver.getSolution();
  if(verbose >= 2) std::cout << "Incircle length: " << delSol.getLength() << std::endl;

  setPaths(delSol);
  if(verbose >= 2) std::cout << "Improved length to Delaunay: " << delSol.getLength() << std::endl;
  // goodSols.push_back(delSol);

  for(double exp = PARM.minExp; exp <= PARM.maxExp; exp += PARM.incExp) {
    if(verbose >= 2) std::cout << "Exponent: " << exp << std::endl;
    if(verbose >= 2) std::cout << "Finding low crossing center: ";
    auto baseTriangulation = delSolver.getSolution().getCenter();
    CenterSolver<Number> csolver(baseTriangulation,inst.triangulations);
    csolver.setExponent(exp);
    csolver.greedy();

    Solution<Number> newSol(inst);
    newSol.setCenter(csolver.getSolution());
    setPaths(newSol);
    goodSols.push_back(std::make_tuple(newSol.getLength(), newSol.sortedLastFlips(), newSol.getFlipVectors()));
  }

  pruneSols(goodSols,PARM.littlePoolSize); // To sort
  if(delSol.getLength() <= std::get<0>(goodSols.back())) {
    goodSols.push_back(std::make_tuple(delSol.getLength(), delSol.sortedLastFlips(), delSol.getFlipVectors()));
  }

  return goodSols;
}


int main(int argc, char **argv) {
  if(!parse(argc, argv))
    return 1;
  std::stringstream parmio;
  parmio << boost::pfr::io(PARM);
  std::cout << "Parameters: " << parmio.str() << std::endl;

  Instance<Number> inst(PARM.fn);
  inst.meta["command"] = concatenate(argc, argv);
  std::string parmstr = parmio.str();
  std::replace( parmstr.begin(), parmstr.end(), '\"', '\'');
  inst.meta["parameters"] = parmstr;

  auto goodSols = manySolutions(inst);
  Solution<Number> sol = Solution(inst,getBest(goodSols));
  sol.save();

  improveCenter(sol, goodSols, PARM.timeout, PARM.sat);
  std::cout << "Best length to a center: " << sol.getLength() << std::endl;

  return 0;
}
