#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>

#include "../src/geoPH2.h"
#include "../src/geoPH3.h"
#include "../src/geoPHd.h"

namespace nb = nanobind;
using namespace nb::literals;

using PointCloud2 = nb::ndarray<double, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>;
using PointCloud3 = nb::ndarray<double, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>;
using PointCloudD = nb::ndarray<double, nb::shape<-1,-1>, nb::c_contig, nb::device::cpu>;

using PersistenceDiagramPoint = std::array<double,2>;
using PersistenceDiagram      = std::vector<PersistenceDiagramPoint>;
using MDPersistenceDiagram    = std::vector<PersistenceDiagram>;

using ContiguousPersistenceDiagram   = nb::ndarray<double, nb::numpy, nb::shape<-1,2>>;
using MDContiguousPersistenceDiagram = std::vector<ContiguousPersistenceDiagram>;

using PersistencePair         = std::pair<gph::Simplex, gph::Simplex>;
using PersistencePairs        = std::vector<PersistencePair>;
using MDPersistencePairs      = std::vector<PersistencePairs>;

using GudhiPersistenceDiagramPoint = std::pair<unsigned,std::pair<double,double>>;
using GudhiPersistenceDiagram      = std::vector<GudhiPersistenceDiagramPoint>;

using Generators1 = std::vector<gph::Generator1>;
using Generators2 = std::vector<gph::Generator2>;

gph::Simplex sorted(gph::Simplex s) {
  std::sort(s.begin(), s.end());
  return s;
}

inline void convert(gph::MultidimensionalDiagram const& diagram, MDPersistenceDiagram &output) {
  output.resize(diagram.size());
  for (unsigned dim = 0; dim < diagram.size(); ++dim) {
    output[dim].reserve(diagram[dim].size());
    for (auto const& [b,d] : diagram[dim])
      output[dim].push_back({b.second, d.second});
    std::sort(output[dim].begin(), output[dim].end());
  }
}

inline void convert(gph::MultidimensionalDiagram const& diagram, MDContiguousPersistenceDiagram &output) {
  MDPersistenceDiagram temp;
  convert(diagram, temp);
  output.resize(diagram.size());
  for (unsigned dim = 0; dim < diagram.size(); ++dim) {
    auto *data = new PersistenceDiagram(std::move(temp[dim]));
    nb::capsule owner(data, [](void *p) noexcept {
        delete static_cast<PersistenceDiagram *>(p);
    });
    output[dim] = ContiguousPersistenceDiagram(data->empty() ? nullptr : data->data()->data(), {data->size(), 2}, owner);
  }
}

inline void convert(gph::MultidimensionalDiagram const& diagram, MDPersistencePairs &output) {
  output.resize(diagram.size());
  for (unsigned dim = 0; dim < diagram.size(); ++dim) {
    output[dim].reserve(diagram[dim].size());
    for (auto const& [b,d] : diagram[dim]) {
      if (d.second < gph::inf)
        output[dim].emplace_back(sorted(b.first), sorted(d.first));
    }
  }
}

inline void convert(gph::MultidimensionalDiagram const& diagram, GudhiPersistenceDiagram &output) {
  for (unsigned dim = 0; dim < diagram.size(); ++dim) {
    for (auto const& [b,d] : diagram[dim])
      output.emplace_back(dim, std::make_pair(b.second, d.second));
    std::sort(output.begin(), output.end(), [](GudhiPersistenceDiagramPoint const& a, GudhiPersistenceDiagramPoint const& b) {
      return a.second < b.second;
    });
  }
}

constexpr unsigned DIM_MAX = 8;
template <unsigned DIM=4>
void findDimension(PointCloudD const& pc, gph::MultidimensionalDiagram &diagram) {
  if constexpr(DIM <= DIM_MAX) {
    if (pc.shape(1) == DIM)
      gph::runDelaunayRipsPersistenceDiagram<DIM>(pc.data(), pc.shape(0), diagram);
    else
      findDimension<DIM+1>(pc, diagram);
  }
  else
    gph::runDelaunayRipsPersistenceDiagram(pc.data(), pc.shape(0), pc.shape(1), diagram);
}

class GeoPHD {
public:
  explicit GeoPHD(PointCloudD const& pc) {
    if (pc.shape(1) == 2) {
      gph::DRPersistence2 frpd(pc.data(), pc.shape(0));
      frpd.computeDelaunayRips0And1Persistence(diagram);
    }
    else if (pc.shape(1) == 3)
      gph::runDelaunayRipsPersistenceDiagram3(pc.data(), pc.shape(0), diagram);
    else
      findDimension(pc, diagram);
  }

  template <typename Diagram>
  [[nodiscard]] Diagram get() const {
    Diagram res(0);
    convert(diagram, res);
    return res;
  }

  [[nodiscard]] std::variant<MDPersistenceDiagram, GudhiPersistenceDiagram, MDContiguousPersistenceDiagram> getPersistenceDiagram(std::string const& format) const {
    if (format == "gudhi")
      return get<GudhiPersistenceDiagram>();
    if (format == "contiguous")
      return get<MDContiguousPersistenceDiagram>();
    return get<MDPersistenceDiagram>();
  }

private:
  gph::MultidimensionalDiagram diagram;
};

template <typename Diagram>
Diagram delaunayRipsPersistence(PointCloudD const& pc) {
  gph::MultidimensionalDiagram diagram;
  if (pc.shape(1) == 2) {
    gph::DRPersistence2 frpd(pc.data(), pc.shape(0));
    frpd.computeDelaunayRips0And1Persistence(diagram);
  }
  else if (pc.shape(1) == 3)
    gph::runDelaunayRipsPersistenceDiagram3(pc.data(), pc.shape(0), diagram);
  else
    findDimension(pc, diagram);
  Diagram res(0);
  convert(diagram, res);
  return res;
}

[[nodiscard]] std::variant<MDPersistenceDiagram, GudhiPersistenceDiagram, MDContiguousPersistenceDiagram> delaunayRipsPersistenceDiagram(PointCloudD const& pc, std::string const& format) {
  if (format == "gudhi")
    return delaunayRipsPersistence<GudhiPersistenceDiagram>(pc);
  if (format == "contiguous")
    return delaunayRipsPersistence<MDContiguousPersistenceDiagram>(pc);
  return delaunayRipsPersistence<MDPersistenceDiagram>(pc);
}

template <typename Diagram>
Diagram ripsPersistence2(PointCloud2 const& pc) {
  gph::MultidimensionalDiagram diagram;
  gph::DRPersistence2 frpd(pc.data(), pc.shape(0));
  frpd.computeRips0And1Persistence(diagram);
  Diagram res(0);
  convert(diagram, res);
  return res;
}

[[nodiscard]] std::variant<MDPersistenceDiagram, GudhiPersistenceDiagram, MDContiguousPersistenceDiagram> ripsPersistenceDiagram2(PointCloud2 const& pc, std::string const& format) {
  if (format == "gudhi")
    return ripsPersistence2<GudhiPersistenceDiagram>(pc);
  if (format == "contiguous")
    return ripsPersistence2<MDContiguousPersistenceDiagram>(pc);
  return ripsPersistence2<MDPersistenceDiagram>(pc);
}

Generators1 delaunayRipsPersistenceGenerators2(PointCloud2 const& pc) {
  gph::MultidimensionalDiagram diagram;
  Generators1 generators1;
  gph::DRPersistence2 frpd(pc.data(), pc.shape(0));
  frpd.computeDelaunayRips0And1Persistence(diagram);
  frpd.exportRips1Generators(generators1);
  return generators1;
}

Generators1 ripsPersistenceGenerators2(PointCloud2 const& pc) {
  gph::MultidimensionalDiagram diagram;
  Generators1 generators1;
  gph::DRPersistence2 frpd(pc.data(), pc.shape(0));
  frpd.computeRips0And1Persistence(diagram);
  frpd.exportRips1Generators(generators1);
  return generators1;
}

std::pair<Generators1, Generators2> delaunayRipsPersistenceGenerators3(PointCloud3 const& pc) {
  gph::MultidimensionalDiagram diagram;
  Generators1 generators1;
  Generators2 generators2;
  gph::runDelaunayRipsPersistenceDiagram3(pc.data(), pc.shape(0), diagram, generators1, generators2);
  return std::make_pair(generators1, generators2);
}



NB_MODULE(geoph_impl, m) {
  m.doc() = "A Delaunay-Rips persistent homology package for low-dimensional point clouds";

  m.def("ripsPersistenceDiagram2",
        &ripsPersistenceDiagram2,
        "X"_a,
        "format"_a = "default");

  m.def("ripsPersistencePairs2",
        &ripsPersistence2<MDPersistencePairs>,
        "X"_a);

  m.def("ripsPersistenceGenerators2",
        &ripsPersistenceGenerators2,
        "X"_a);

  m.def("delaunayRipsPersistenceGenerators2",
        &delaunayRipsPersistenceGenerators2,
        "X"_a);

  m.def("delaunayRipsPersistenceGenerators3",
        &delaunayRipsPersistenceGenerators3,
        "X"_a);

  m.def("delaunayRipsPersistenceDiagram",
        &delaunayRipsPersistenceDiagram,
        "X"_a,
        "format"_a = "default");

  m.def("delaunayRipsPersistencePairs",
        &delaunayRipsPersistence<MDPersistencePairs>,
        "X"_a);

  nb::class_<GeoPHD>(m, "GeoPHD")
      .def(nb::init<PointCloudD>(), "X"_a)
      .def("delaunayRipsPersistenceDiagram", &GeoPHD::getPersistenceDiagram, "format"_a="default")
      .def("delaunayRipsPersistencePairs", &GeoPHD::get<MDPersistencePairs>);
}