#pragma once

#include "geoPHUtils.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_id_2.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>

#include <ranges>

namespace gph {

  struct FiltratedQuadEdge {
    Edge e;
    int f1;
    int f2;
    value_t d;
  };

  class DRPersistence2 {
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Vb = CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K>;
    using Fb = CGAL::Triangulation_face_base_with_id_2<K>;
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    using Delaunay = CGAL::Delaunay_triangulation_2<K, Tds>;
    using Point = Delaunay::Point_2;
    using Traits = CGAL::Search_traits_2<K>;
    using Tree = CGAL::Kd_tree<Traits>;
    using Fuzzy_sphere = CGAL::Fuzzy_sphere<Traits>;

  public:
    explicit DRPersistence2(const std::vector<std::vector<double>> &_points) : N_p(_points.size()) {
      points.resize(N_p);
      for (unsigned i=0; i<N_p; ++i)
        points[i] = std::make_pair(Point(_points[i][0], _points[i][1]), i);
      computeDelaunay();
    }

    explicit DRPersistence2(double* data, int n) : N_p(n) {
      points.resize(N_p);
      for (unsigned i=0; i<N_p; ++i)
        points[i] = std::make_pair(Point(data[2*i], data[2*i+1]), i);
      computeDelaunay();
    }

    void compute0Persistence(Diagram &ph0, bool parallelSort = false) {
      ph0 = Diagram(0);

      //keep only the Urquhart edges, maintain polygon structure
      UnionFind UF(N_f);
      std::vector<FiltratedEdge> max_delaunay(N_f, FiltratedEdge{{-1,-1},0.});
      computeUrquhart(UF, max_delaunay, parallelSort);

      // compute EMST with Kruskal algorithm
      UnionFind UF_p(N_p);
      for (FiltratedQuadEdge const& e : urquhart) {
        if(UF_p.find(e.e.first) != UF_p.find(e.e.second)) { //we know e is a EMST edge
          UF_p.merge(e.e.first, e.e.second);
          ph0.emplace_back(FiltratedSimplex{{-1},0.}, FiltratedSimplex{{e.e.first, e.e.second},e.d});
        }
      }
      ph0.emplace_back(FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); //infinite pair
    }

    void computeDelaunayRips0And1Persistence(MultidimensionalDiagram &ph, bool parallelSort = false) {
      ph = MultidimensionalDiagram(2);

      //keep only the Urquhart edges, maintain polygon structure
      UnionFind UF(N_f);
      death_poly.resize(N_f, {{-1,-1},0.});
      computeUrquhart(UF, death_poly, parallelSort);

      // compute EMST with Kruskal algorithm
      UnionFind UF_p(N_p);
      std::vector<FiltratedQuadEdge> critical(0);
      for (FiltratedQuadEdge const& e : urquhart) {
        if(UF_p.find(e.e.first) != UF_p.find(e.e.second)) { //we know e is a EMST edge
          UF_p.merge(e.e.first, e.e.second);
          ph[0].emplace_back(FiltratedSimplex{{-1},0.}, FiltratedSimplex{{e.e.first, e.e.second},e.d});
        }
        else // we know e is a UG-EMST edge, i.e. it creates a cycle
          critical.push_back(e);
      }
      ph[0].emplace_back(FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); //infinite pair

      compute1PH(critical, UF, ph);
    }

    void computeRips0And1Persistence(MultidimensionalDiagram &ph, bool parallelSort = false, bool parallelMML = false){
      ph = MultidimensionalDiagram(2);

      //compute k-d tree
      Tree tree;
      for (auto const& p : points)
        tree.insert(p.first);

      //keep only the Urquhart edges, maintain polygon structure
      UnionFind UF(N_f);
      std::vector<FiltratedEdge> max_delaunay(N_f, FiltratedEdge{{-1,-1},0.});
      computeUrquhart(UF, max_delaunay, parallelSort);

      // compute EMST with Kruskal algorithm
      UnionFind UF_p(N_p);
      std::vector<FiltratedQuadEdge> critical(0); //RNG-EMST edges

      for (FiltratedQuadEdge const& e : urquhart) {
        if(UF_p.find(e.e.first) != UF_p.find(e.e.second)) { //we know e is a EMST edge
          UF_p.merge(e.e.first, e.e.second);
          rng.push_back(e);
          ph[0].emplace_back(FiltratedSimplex{{-1},0.}, FiltratedSimplex{{e.e.first, e.e.second},e.d});
        }
        else { //we know e is a UG-EMST edge
          //check if e is RNG
          if (isLensEmpty(points[e.e.first].first, points[e.e.second].first, tree, e.d)) { //RNG edge
            critical.push_back(e);
            rng.push_back(e);
          }
          else { //not RNG edge : merge neighboring polygons
            const int poly1 = UF.find(e.f1);
            const int poly2 = UF.find(e.f2);
            UF.merge(poly1, poly2);
            max_delaunay[UF.find(poly1)] = max(FiltratedEdge{e.e, e.d}, max(max_delaunay[poly1], max_delaunay[poly2]));
          }
        }
      }
      ph[0].emplace_back(FiltratedSimplex{{-1}, 0.}, FiltratedSimplex{{-1}, inf}); //infinite pair

      // polygon reindexation
      std::vector<int> index_polys (0);
      reindexPolygons(UF, max_delaunay, index_polys);
      for (FiltratedQuadEdge &e : rng) {
        e.f1 = index_polys[UF.find(e.f1)]; //final path compression
        e.f2 = index_polys[UF.find(e.f2)]; //final path compression
      }
      for (FiltratedQuadEdge &e : critical) {
        e.f1 = index_polys[UF.find(e.f1)]; //final path compression
        e.f2 = index_polys[UF.find(e.f2)]; //final path compression
      }

      computePolygonRipsDeath(parallelMML, UF, index_polys);

      UnionFind UF_poly (death_poly.size());
      compute1PH(critical, UF_poly, ph);
    }

    void exportRips1Generators(std::vector<Generator1> &generators) const {
      const unsigned N_polys = death_poly.size();

      std::vector<Generator1> generators_(N_polys);
      for (unsigned i=0; i<death_poly.size(); ++i)
        generators_[i].second = {birth_poly[i], death_poly[i].d};
      for (const FiltratedQuadEdge &e : rng) {
        if (e.f1 != e.f2) {
          if(death_poly[e.f1].d != inf)
            generators_[e.f1].first.push_back(e.e);
          if(death_poly[e.f2].d != inf)
            generators_[e.f2].first.push_back(e.e);
        }
      }

      generators.resize(0);
      for (auto const& g : generators_) {
        if (!g.first.empty())
          generators.push_back(g);
      }
    }

  private:
    const unsigned N_p;
    unsigned N_f {0};
    std::vector<std::pair<Point,unsigned>> points;
    Delaunay del;

    std::vector<FiltratedQuadEdge> urquhart;
    std::vector<FiltratedQuadEdge> rng;
    std::vector<FiltratedEdge> death_poly;
    std::vector<double> birth_poly;

    //common to Delaunay-Rips and Rips

    void computeDelaunay() {
      //compute Delaunay and initialize triangle IDs (including infinite triangles to deal with the convex hull)
      del = Delaunay(points.begin(), points.end());
      int k = 0;
      for(auto const& f : del.all_face_handles())
        f->id() = k++;
      N_f = k;
    }

    void computeUrquhart(UnionFind& UF, std::vector<FiltratedEdge>& max_delaunay, bool parallelSort) {
      for(Delaunay::Edge const& e : del.finite_edges()) {
        const unsigned a = e.first->vertex((e.second+1)%3)->info();
        const unsigned b = e.first->vertex((e.second+2)%3)->info();
        const double d2 = CGAL::squared_distance(points[a].first, points[b].first);
        const double d = sqrt(d2);
        const Delaunay::Edge e_m = del.mirror_edge(e);

        // check if edge is Urquhart
        bool is_urquhart = true;
        const Point p_k = e.first->vertex(e.second)->point();
        if (!del.is_infinite(e.first) && CGAL::squared_distance(points[a].first, p_k) < d2 && CGAL::squared_distance(points[b].first, p_k) < d2)
          is_urquhart = false;
        else {
          const Point p_l = e_m.first->vertex(e_m.second)->point();
          if (!del.is_infinite(e_m.first) && CGAL::squared_distance(points[a].first, p_l) < d2 && CGAL::squared_distance(points[b].first, p_l) < d2)
            is_urquhart = false;
        }

        //maintain polygon structure
        if (is_urquhart) //UG edge
          urquhart.push_back(FiltratedQuadEdge{{a,b}, e.first->id(), e_m.first->id(), d});
        else { //not UG edge: maintain UF
          const int poly1 = UF.find(e.first->id());
          const int poly2 = UF.find(e_m.first->id());
          max_delaunay[UF.mergeRet(poly1, poly2)] = max(FiltratedEdge{{a,b}, d}, max(max_delaunay[poly1], max_delaunay[poly2]));
        }

        //detect infinite polygons (going beyond convex hull)
        if (del.is_infinite(e.first))
          max_delaunay[UF.find(e.first->id())].d = inf;
        else if (del.is_infinite(e_m.first))
          max_delaunay[UF.find(e_m.first->id())].d = inf;
      }

      if (parallelSort)
        std::sort(GPH_EXECUTION_POLICY urquhart.begin(), urquhart.end(), [](FiltratedQuadEdge e1, FiltratedQuadEdge e2){return e1.d < e2.d;});
      else
        std::sort(urquhart.begin(), urquhart.end(), [](FiltratedQuadEdge e1, FiltratedQuadEdge e2){return e1.d < e2.d;});
    }

    void compute1PH(std::vector<FiltratedQuadEdge> const& critical, UnionFind &UF, MultidimensionalDiagram &ph) {
      std::vector<int> latest(death_poly.size());
      std::iota(latest.begin(), latest.end(), 0);
      birth_poly.resize(death_poly.size());

      for (FiltratedQuadEdge const& e : std::ranges::reverse_view(critical)) {
        const int v1 = UF.find(e.f1);
        const int v2 = UF.find(e.f2);
        UF.merge(v1, v2);

        const int latest1 = latest[v1];
        const int latest2 = latest[v2];

        if (death_poly[latest1].d < death_poly[latest2].d) {
          ph[1].emplace_back(FiltratedSimplex{{e.e.first, e.e.second}, e.d},
                             FiltratedSimplex{{death_poly[latest1].e.first, death_poly[latest1].e.second}, death_poly[latest1].d});
          birth_poly[latest1] = e.d;
          latest[UF.find(v1)] = latest2;
        }
        else {
          ph[1].emplace_back(FiltratedSimplex{{e.e.first, e.e.second}, e.d},
                             FiltratedSimplex{{death_poly[latest2].e.first, death_poly[latest2].e.second}, death_poly[latest2].d});
          birth_poly[latest2] = e.d;
          latest[UF.find(v1)] = latest1;
        }
      }
    }

    //specific to Rips
    [[nodiscard]] static bool isLensEmpty(Point const& p1, Point const& p2, Tree const& tree, double const& d) {
      const Fuzzy_sphere fs (p1, d);
      std::vector<Point> ball;
      tree.search(back_inserter(ball), fs);
      for (Point const& p: ball) {
        if (p != p1 && p != p2) {
          if (CGAL::squared_distance(p, p2) < d * d)
            return false;
        }
      }
      return true;
    }

    [[nodiscard]] static bool isRightSemiLensEmpty(Point const& p1, Point const& p2, Tree const& tree) {
      const double d2 = CGAL::squared_distance(p1, p2);
      const Fuzzy_sphere fs (p1, sqrt(d2));
      std::vector<Point> ball;
      tree.search(back_inserter(ball), fs);
      for (Point const& p: ball) {
        if (CGAL::squared_distance(p, p2) < d2 && CGAL::right_turn(p, p1, p2))
          return false;
      }
      return true;
    }

    void reindexPolygons(UnionFind const& UF, std::vector<FiltratedEdge> const& max_delaunay, std::vector<int>& index_polys) {
      //reindex polygons and find those that are beyond convex hull
      index_polys = std::vector<int> (N_f, -1);
      int N_polys = 0;
      for (unsigned x=0; x<N_f; ++x) {
        if (UF.isRoot(x)) {
          index_polys[x] = N_polys;
          death_poly.push_back(max_delaunay[x]);
          N_polys++;
        }
      }
    }

    void computePolygonRipsDeath(bool parallel, UnionFind &UF, std::vector<int> const& index_polys) {
      const unsigned N_polys = death_poly.size();
      std::vector<std::vector<int>> vertices(N_polys);

      //store coarsely vertices in each polygon
      for (const FiltratedQuadEdge &e : rng) {
        vertices[e.f1].push_back(e.e.first);
        vertices[e.f1].push_back(e.e.second);
        vertices[e.f2].push_back(e.e.first);
        vertices[e.f2].push_back(e.e.second);
      }
      std::vector<Delaunay::Face_handle> repr_face (N_polys);
      for (auto const& fh : del.finite_face_handles()) {
        if (UF.isRoot(fh->id()))
          repr_face[index_polys[fh->id()]] = fh;
      }

      //loop on polygons
#ifdef ENABLE_OPENMP
#pragma omp parallel for if(parallel)
#endif // ENABLE_OPENMP
      for (int poly=0; poly<N_polys; ++poly) { //int for MSVC
        if (death_poly[poly].d != inf) { // deal only with true polygons
          const double bound_max = death_poly[poly].d;
          const double bound_min = sqrt(3)/2 * bound_max;

          //find involved vertices
          sort(vertices[poly].begin(), vertices[poly].end());
          auto last = unique(vertices[poly].begin(), vertices[poly].end());
          vertices[poly].erase(last, vertices[poly].end());

          std::vector<Point> pts(0);
          std::vector<int> global(0);
          for (int const& v : vertices[poly]) {
            pts.push_back(points[v].first);
            global.push_back(v);
          }
          const unsigned N_pts = pts.size();
          const Tree loc_tree(pts.begin(), pts.end());

          //enumerate edges and do stuff
          FiltratedEdge death {{-1, -1}, inf};
          for (unsigned i = 0; i<N_pts-1; ++i) {
            for (unsigned j = i+1; j<N_pts; ++j) {
              const double d2 = CGAL::squared_distance(pts[i], pts[j]);
              const double d = sqrt(d2);
              if (bound_min <= d && d <= bound_max && d < death.d) { // edge length interval allowing to be a candidate
                const Point p_i = pts[i];
                const Point p_j = pts[j];

                const Fuzzy_sphere fs (p_i, d);
                std::vector<Point> ball;
                loc_tree.search(back_inserter(ball), fs);

                //check if edge is a 2-edge (i.e. both semi-lens are not empty)
                std::vector<Point> lens_l;
                std::vector<Point> lens_r;
                for (auto const& p: ball) {
                  if (p != p_i && p != p_j) {
                    if (CGAL::squared_distance(p, p_j) < d2) {
                      if (CGAL::right_turn(p_i, p_j, p))
                        lens_r.push_back(p);
                      else
                        lens_l.push_back(p);
                    }
                  }
                }

                //if 2-edge, check if it is expandable
                if (!lens_l.empty() && !lens_r.empty()) { // 2-edge
                  bool l_expandable = false, r_expandable = false;
                  for (Point const& p : lens_l) {
                    if (!isRightSemiLensEmpty(p_i, p, loc_tree))
                      continue;
                    else if (isRightSemiLensEmpty(p, p_j, loc_tree)) {
                      l_expandable = true;
                      break;
                    }
                  }
                  if (!l_expandable)
                    continue;
                  for (Point const& p : lens_r) {
                    if (!isRightSemiLensEmpty(p_j, p, loc_tree))
                      continue;
                    else if (isRightSemiLensEmpty(p, p_i, loc_tree)) {
                      r_expandable = true;
                      break;
                    }
                  }
                  if (r_expandable) { //expandable 2-edge, check if it is a diagonal of the current polygon (costly)
                    if (index_polys[UF.find(del.locate(CGAL::midpoint(p_i,p_j), repr_face[poly])->id())] == int(poly))
                      death = FiltratedEdge{{global[i], global[j]}, d};
                  }
                }
              }
            }
          }
          death_poly[poly] = death;
        }
      }
    }

  };

}
