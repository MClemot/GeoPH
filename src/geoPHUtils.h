#pragma once

#include <array>
#include <limits>
#include <vector>

#include <boost/container_hash/hash.hpp>
#include <boost/version.hpp>
#if ((BOOST_VERSION / 100) % 1000) >= 83
#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#define GPH_CONCURRENT_HASHTABLE_AVAILABLE
#endif

#ifdef ENABLE_TBB
#include <tbb/concurrent_vector.h>
#include <tbb/global_control.h>
#endif

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

#ifdef __cpp_lib_execution
#include <execution>
#define GPH_EXECUTION_POLICY std::execution::par,
#else
#define GPH_EXECUTION_POLICY
#endif

#if defined(ENABLE_OPENMP) and defined(ENABLE_TBB) and defined(GPH_CONCURRENT_HASHTABLE_AVAILABLE)
#define GPH_PARALLEL
#endif

#include "../deps/dset.h"
#include "../deps/unordered_dense.h"

namespace gph {
  using id_t = int;
  using value_t = double;
  constexpr value_t inf = std::numeric_limits<value_t>::infinity();

  template <unsigned DIM>
  using PointD = std::conditional_t<DIM==0, std::vector<value_t>, std::array<value_t,DIM>>;

  template <unsigned DIM>
  using PointCloud = std::vector<PointD<DIM>>;

  using Simplex = std::vector<id_t>;
  using FiltratedSimplex = std::pair<Simplex, value_t>;
  using PersistencePair = std::pair<FiltratedSimplex, FiltratedSimplex>;
  using Diagram = std::vector<PersistencePair>;
  using MultidimensionalDiagram = std::vector<Diagram>;

  using Edge = std::pair<id_t, id_t>;
  struct FiltratedEdge {
    Edge e;
    value_t d;
  };

  inline FiltratedEdge max(FiltratedEdge const& a, FiltratedEdge const& b) {
    if (a.d > b.d)
      return a;
    return b;
  }

  using Facet = std::array<id_t, 3>;

  using Generator1 = std::pair<std::vector<Edge>,  std::pair<value_t,value_t>>;
  using Generator2 = std::pair<std::vector<Facet>, std::pair<value_t,value_t>>;

#ifdef GPH_CONCURRENT_HASHTABLE_AVAILABLE
  template <typename X, typename Y> using ConcurrentHashMap = boost::concurrent_flat_map<X, Y>;
  template <typename X, typename Y> using SequentialHashMap = boost::unordered_flat_map<X, Y>;
#endif

  template <typename X, typename Y> using HashMap = ankerl::unordered_dense::map<X, Y, boost::hash<X>>;
  template <typename X> using HashSet = ankerl::unordered_dense::set<X>;

  class UnionFind {
    std::vector<int> parent;
    std::vector<unsigned char> rank;
  public:
    explicit UnionFind(const unsigned n) {
      parent.resize(n);
      rank.resize(n, 0);
      for (unsigned i = 0; i < n; i++)
        parent[i] = i;
    }

    int find(const int x) {
      if (parent[x] == x)
        return x;
      return parent[x] = find(parent[x]); // path compression
    }

    void merge(const int x, const int y) {
      const int rootX = find(x);
      const int rootY = find(y);
      if (rootX != rootY) {
        if (rank[rootX] > rank[rootY])
          parent[rootY] = rootX;
        else if (rank[rootX] < rank[rootY])
          parent[rootX] = rootY;
        else {
          parent[rootY] = rootX;
          rank[rootX]++;
        }
      }
    }

    int mergeRet(const int x, const int y) {
      const int rootX = find(x);
      const int rootY = find(y);
      if (rootX != rootY) {
        if (rank[rootX] > rank[rootY])
          return parent[rootY] = rootX;
        else if (rank[rootX] < rank[rootY])
          return parent[rootX] = rootY;
        else {
          rank[rootX]++;
          return parent[rootY] = rootX;
        }
      }
      return rootX;
    }

    [[nodiscard]] bool isRoot(const int x) const {
      return parent[x] == x;
    }
  };

}