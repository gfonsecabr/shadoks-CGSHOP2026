#pragma once
#include "tools.hpp"
#include "point.hpp"
#include "segment.hpp"

// Returns first element of v1 that is not in v2
template<class T1, class T2>
auto missing(T1 v1, T2 v2) {
  for(const auto e1 : v1) {
    bool found = false;
    for(const auto e2 : v2) {
      if(e1 == e2) {
        found = true;
        break;
      }
    }
    if(!found) {
      return e1;
    }
  }
  throw std::runtime_error("Error: nothing is missing!");
}

template <class Number>
class Triangle {
  std::array<Point<Number>,3> points;
public:
  Triangle() {}

  Triangle(const Point<Number> &p1, const Point<Number> &p2, const Point<Number> &p3) :
    points{p1,p2,p3}
    {}

  Triangle<Number> operator+(Point<Number> p) const {
    return Triangle(points[0]+p,points[1]+p,points[2]+p);
  }

  Triangle<Number> operator-(Point<Number> p) const {
    return Triangle(points[0]-p,points[1]-p,points[2]-p);
  }


  auto operator<=>(const Triangle<Number> &t) const {
    return points <=> t.points;
  }

  Point<Number> operator[](int i) const {
    return points[i];
  }


  Number twiceSignedArea() const {
    return (points[0].y() - points[2].y()) * (points[1].x() - points[0].x()) -
           (points[0].x() - points[2].x()) * (points[1].y() - points[0].y());
  }

  int orientation() const {
    Number val = twiceSignedArea();
    return (val > 0) - (val < 0);
  }

  __int128  incircleValue(Point<Number> pd) const {
    auto [pa,pb,pc] = points;
    __int128 adx = pa.x() - pd.x();
    __int128 ady = pa.y() - pd.y();
    __int128 bdx = pb.x() - pd.x();
    __int128  bdy = pb.y() - pd.y();
    __int128  cdx = pc.x() - pd.x();
    __int128  cdy = pc.y() - pd.y();
    __int128  abdet = adx * bdy - bdx * ady;
    __int128  bcdet = bdx * cdy - cdx * bdy;
    __int128  cadet = cdx * ady - adx * cdy;
    __int128  alift = adx * adx + ady * ady;
    __int128  blift = bdx * bdx + bdy * bdy;
    __int128  clift = cdx * cdx + cdy * cdy;
    return alift * bcdet + blift * cadet + clift * abdet;
  }

  // Tests if pd is (1) inside, (0) on the boundary, or (-1) outside
  // the circumcircle of the triangle
  int incircle(Point<Number> pd) const {
    int o = orientation();
    if(o == 0) {
      return 0; // No idea what to return when the triangle points are colinear
    }
    __int128 ic = incircleValue(pd);
    if(ic == 0) {
      return 0;
    }

    if((o < 0 && ic > 0) || (o > 0 && ic < 0)) {
      return 1;
    }

    return -1;
  }

  bool contains(const Point<Number> &p) const {
    Segment<Number> s1 = {points[0],points[1]};
    Segment<Number> s2 = {points[1],points[2]};
    Segment<Number> s3 = {points[2],points[0]};
    int o1, o2;
    return (o1 = s1.orientation(p)) != 0 && o1 == (o2 = s2.orientation(p)) && o2 == s3.orientation(p);
  }

  Triangle<Number> sorted() const {
    Triangle<Number> ret(*this);
    std::sort(ret.points.begin(),ret.points.end());
    return ret;
  }

  bool canFlip(const Triangle<Number> &t) const {
    std::array<Point<Number>,3> pts1 = points, pts2 = t.points;
    std::sort(pts1.begin(),pts1.end());
    std::sort(pts2.begin(),pts2.end());
    std::vector<Point<Number>> isec;
    std::set_intersection(pts1.begin(), pts1.end(), pts2.begin(), pts2.end(),
                          std::back_inserter(isec));
    if(isec.size() != 2)
      return false;


    Segment oldSeg(isec[0],isec[1]);
    Segment newSeg(missing(pts1, isec), missing(pts2, isec));

    return oldSeg.cross(newSeg);
  }
};

template <class Number>
std::ostream &operator<<(std::ostream &os, Triangle<Number> const &t) {
  return os << "{" << t[0] << "," << t[1]<< "," << t[2]<< "}";
}

namespace std {
  template <class Number> struct hash<Triangle<Number>> {
    size_t operator()(const Triangle<Number> &t) const {
      size_t seed = 0;
      for(auto p : t.points) {
        boost::hash_combine(seed, std::hash<Point<Number>>{}(p));
      }
      return seed;
    }
  };
}

