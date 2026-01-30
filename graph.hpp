#pragma once
#include "tools.hpp"

template <class Vertex>
class Graph {
  std::unordered_map<Vertex, std::unordered_set<Vertex>> adj;

public:
  Graph() {
  }

  Graph(const std::vector<std::array<Vertex,2>> &edges) {
    for(auto [u,v] : edges)
      addEdge(u,v);
  }

  void addVertex(Vertex v) {
    adj[v];
  }

  void addEdge(Vertex u, Vertex v) {
    if(u != v) {
      adj[u].insert(v);
      adj[v].insert(u);
    }
  }

  bool containsVertex(Vertex v) const {
    return adj.count(v) != 0;
  }

  bool containsEdge(Vertex u, Vertex v) const {
    return containsVertex(u) && adj.at(u).count(v) != 0;
  }

  int degree(Vertex v) const {
    return containsVertex(v) ? adj.at(v).size() : -1;
  }

  int maxDegree() const {
    int ret = -1;
    for(const auto &[v, neigh] : adj) {
      ret = std::max(ret,(int)neigh.size());
    }

    return ret;
  }

  int countVertices() const {
    return adj.size();
  }

  int countEdges() const {
    int ret = 0;
    for(const auto &[v, neigh] : adj) {
      ret += neigh.size();
    }
    assert(ret % 2 == 0);
    return ret / 2;
  }

  void removeEdge(Vertex u, Vertex v) {
    if(containsEdge(u,v)) {
      adj.at(u).erase(v);
      adj.at(v).erase(u);
    }
  }

  void removeVertex(Vertex v) {
    if(containsVertex(v)) {
      std::unordered_set<Vertex> neigh(adj.at(v)); // Copy because changing the data structure invalidades the iterator
      for(Vertex u : neigh) {
        removeEdge(u,v);
      }
      adj.erase(v);
    }
  }

  void clear() {
    adj.clear();
  }

  std::vector<Vertex> vertices() const {
    std::vector<Vertex> ret;
    for(const auto &[v, neigh] : adj) {
      ret.push_back(v);
    }
    return ret;
  }

  const std::unordered_set<Vertex> & neighbors(Vertex v) const {
    return adj.at(v);
  }

  std::unordered_set<Vertex> closedNeighbors(Vertex v) const {
    std::unordered_set<Vertex> ret = adj.at(v);
    ret.insert(v);
    return ret;
  }

  std::vector<Vertex> bfs(Vertex v, int maxv = 0) const {
    std::unordered_set<Vertex> visited;
    std::vector<Vertex> ret;
    std::queue<Vertex> fifo;

    if(maxv == 0)
      maxv = countVertices();

    fifo.push(v);
    while(!fifo.empty() && ret.size() < (size_t) maxv) {
      Vertex u = fifo.front();
      fifo.pop();
      if(visited.count(u) != 0)
        continue;
      ret.push_back(u);
      visited.insert(u);
      for(Vertex w : neighbors(u))
        fifo.push(w);
    }

    return ret;
  }

  class iterator {
    Graph<Vertex> *graph;
    class std::unordered_map<Vertex,std::unordered_set<Vertex>>::iterator it;
  public:
    iterator(Graph *g = nullptr) : graph(g){
      if(g != nullptr)
        it = g->adj.begin();
    }
    Vertex operator*() {
      return it->first;
    }
    void operator++() {
      it++;
      if(it == graph->adj.end())
        graph = nullptr;
    }
    bool operator==(iterator other) {
      if(graph == nullptr && other.graph == nullptr)
        return true;
      if(graph == nullptr || other.graph == nullptr)
        return false;

      return it == other.it;
    }
    bool operator!=(iterator other) {
      return !operator==(other);
    }
  };

  Graph<Vertex>::iterator begin() {
    return iterator(this);
  }

  Graph<Vertex>::iterator end() {
    return iterator();
  }

  // List of 3-cycles
  std::vector<std::array<Vertex,3>> triangles() const {
    std::vector<std::array<Vertex,3>> ret;

    for(const auto &[u,nu] : adj) {
      for(Vertex v : nu) {
        if(u < v) {
          for(Vertex w : nu) {
            if(v < w && adj.at(w).contains(v)) {
              ret.push_back({u,v,w});
            }
          }
        }
      }
    }
    return ret;
  }

};

// // https://www.variadic.xyz/post/0120-hashing-tuples/
// template<class T>
// inline void hash_combine(std::size_t& seed, const T& val) {
//   std::hash<T> hasher;
//   seed ^= hasher(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
// }
//
// // https://www.variadic.xyz/post/0120-hashing-tuples/
// namespace std {
//   template<class S, class T>
//   struct hash<std::pair<S, T>> {
//     inline size_t operator()(const std::pair<S, T>& val) const {
//       size_t seed = 0;
//       hash_combine(seed, val.first);
//       hash_combine(seed, val.second);
//       return seed;
//     }
//   };
// }
