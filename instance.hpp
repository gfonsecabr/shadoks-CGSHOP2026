#pragma once
#include "tools.hpp"
#include "point.hpp"
#include "triangulation.hpp"
#include "rapidjson/reader.h"
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"



template <class Number>
struct Instance {
  std::string name;
  mutable std::unordered_map<std::string, std::string> meta;
  std::vector<Point<Number>> points;
  std::vector<Triangulation<Number>> triangulations;
  std::unordered_map<IEdge,std::unordered_set<short>> emptyTrianglesMap;

  Instance() {}

  Instance(std::string fn) {
    meta["start_time"] = timeString();
    meta["author"] = "gdf";
    char hn[80];
    gethostname(hn, 80);
    meta["host"] = std::string(hn);
    rapidjson::Document doc = read_json(fn);
    name = doc["instance_uid"].GetString();
    points = json_points(doc);

    auto edgeSets = json_edges(doc);
    for(auto &triangulation : edgeSets) {
      triangulations.push_back(Triangulation<Number>(&points,triangulation));
    }
    std::cout << "Instance " << name << " with " << points.size() << " points and " << triangulations.size() << " triangulations\n";
  }

  static rapidjson::Document read_json(std::string filename) {
    std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
    if (!in.is_open()) {
      std::cerr << "Error reading " << filename << std::endl;
      assert(false);
      exit(EXIT_FAILURE);
    }

    rapidjson::IStreamWrapper isw {in};
    rapidjson::Document doc {};
    doc.ParseStream(isw);

    if (doc.HasParseError()) {
      std::cerr << "Error  : " << doc.GetParseError()  << std::endl;
      std::cerr << "Offset : " << doc.GetErrorOffset() << std::endl;
      exit(EXIT_FAILURE);
    }
    return doc;
  }

  static std::vector<Number> json_int_vec(rapidjson::Value &values) {
    std::vector<Number> v;
    for (auto &x : values.GetArray()) {
      double d = x.GetDouble();
      v.push_back((int)d);
    }
    return v;
  }

  static std::vector<Point<Number>> json_points(rapidjson::Value &values) {
    std::vector<Number> xs = json_int_vec(values["points_x"]);
    std::vector<Number> ys = json_int_vec(values["points_y"]);
    std::vector<Point<Number>> points;
    for(size_t i = 0; i < xs.size(); i++)
      points.push_back(Point(xs[i],ys[i]));

    return points;
  }

  std::vector<std::vector<IEdge>> json_edges(rapidjson::Value &values) {
    std::vector<std::vector<IEdge>> v;

    for (auto &triangulation : values["triangulations"].GetArray()) {
      v.push_back({});
      for (auto &edge : triangulation.GetArray())
        v.back().push_back({(short)edge[0].GetInt(), (short)edge[1].GetInt()});
    }

    return v;
  }

  bool emptyTriangle(int ia, int ib, int ic) {
    IEdge e = {short(ia),short(ib)};

    if(!emptyTrianglesMap.contains(e)) {
      emptyTrianglesMap[e] = emptyTriangles(e);
    }

    return emptyTrianglesMap.at(e).contains(ic);
  }

protected:
  std::unordered_set<short> emptyTriangles(IEdge e) {
    Segment<Number> s = {points[e[0]], points[e[1]]};

    std::array<std::vector<short>,2> lrPoints;
    for(short i = 0; i < short(points.size()); i++) {
      Point<Number> p_i = points[i];
      Number o = s.orientation(p_i);
      if(o < 0)
        lrPoints[0].push_back(i);
      else
        if(o > 0)
        lrPoints[1].push_back(i);
      // If o == 0, then the points are colinear and there is no empty triangle
    }

    Point s0 = s.source();
    auto cmp = [this,s0](short a, short b) {
      return Segment<Number>(this->points[a],this->points[b]).orientation(s0) > 0;
    };

    std::sort(lrPoints[0].rbegin(), lrPoints[0].rend(), cmp);
    std::sort(lrPoints[1].begin(), lrPoints[1].end(), cmp);

    std::unordered_set<short> ret;

    Point s1 = s.target();
    for(int lr : {0,1}) {
      Point<Number> minpt = s0;
      for(short i : lrPoints[lr]) {
        Point<Number> p_i = this->points[i];
        if(lr == 0) {
          if((minpt == s0 || Segment<Number>(p_i,minpt).orientation(s1) > 0)
             && s.orientation(p_i) != 0) {
            minpt = p_i;
            ret.insert(i);
          }
        }
        else {
          if((minpt == s0 || Segment<Number>(p_i,minpt).orientation(s1) < 0)
             && s.orientation(p_i) != 0) {
            minpt = p_i;
            ret.insert(i);
          }
        }
      }
    }

    return ret;
  }

};
