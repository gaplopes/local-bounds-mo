#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/optional.h>
#include "local_bounds.hpp"

namespace nb = nanobind;
using namespace local_bounds;

NB_MODULE(_core, m) {
    m.doc() = "Python bindings for the C++ local_bounds library.";

    nb::enum_<Objective>(m, "Objective")
        .value("MINIMIZE", Objective::MINIMIZE)
        .value("MAXIMIZE", Objective::MAXIMIZE)
        .export_values();

    nb::class_<Point<double>>(m, "Point")
        .def(nb::init<>())
        .def(nb::init<std::vector<double>>(), nb::arg("coordinates"))
        .def(nb::init<std::string, std::vector<double>>(), nb::arg("id"), nb::arg("coordinates"))
        .def_rw("id", &Point<double>::id)
        .def_rw("coordinates", &Point<double>::coordinates)
        .def("dimensions", &Point<double>::dimensions)
        .def("__str__", &Point<double>::to_string)
        .def("__repr__", &Point<double>::to_string)
        .def("__eq__", [](const Point<double>& self, const Point<double>& other) { return self == other; })
        .def("__ne__", [](const Point<double>& self, const Point<double>& other) { return self != other; });

    nb::class_<LocalBound<double>>(m, "LocalBound")
        .def(nb::init<>())
        .def(nb::init<std::string, std::vector<double>>(), nb::arg("id"), nb::arg("coordinates"))
        .def_rw("id", &LocalBound<double>::id)
        .def_rw("coordinates", &LocalBound<double>::coordinates)
        .def_rw("defining_points", &LocalBound<double>::defining_points)
        .def_rw("defining_point_sets", &LocalBound<double>::defining_point_sets)
        .def("dimensions", &LocalBound<double>::dimensions)
        .def("__str__", &LocalBound<double>::to_string)
        .def("__repr__", &LocalBound<double>::to_string)
        .def("__eq__", [](const LocalBound<double>& self, const LocalBound<double>& other) { return self == other; })
        .def("__ne__", [](const LocalBound<double>& self, const LocalBound<double>& other) { return self != other; });

    nb::class_<BoundSet<double, Objective::MINIMIZE>>(m, "BoundSetMinimize")
        .def(nb::init<const std::vector<double>&>(), nb::arg("reference_point"))
        .def(nb::init<const std::vector<double>&, const std::vector<double>&>(), nb::arg("reference_point"), nb::arg("anti_reference"))
        .def("update_auto", &BoundSet<double, Objective::MINIMIZE>::update_auto, nb::arg("point"))
        .def("update_re", &BoundSet<double, Objective::MINIMIZE>::update_re, nb::arg("point"))
        .def("update_re_enhanced", &BoundSet<double, Objective::MINIMIZE>::update_re_enhanced, nb::arg("point"))
        .def("update_ra_sa", &BoundSet<double, Objective::MINIMIZE>::update_ra_sa, nb::arg("point"))
        .def("update_ra", &BoundSet<double, Objective::MINIMIZE>::update_ra, nb::arg("point"))
        .def("update_naive", &BoundSet<double, Objective::MINIMIZE>::update_naive, nb::arg("point"))
        .def("size", &BoundSet<double, Objective::MINIMIZE>::size)
        .def("dimensions", &BoundSet<double, Objective::MINIMIZE>::dimensions)
        .def("is_in_search_region", &BoundSet<double, Objective::MINIMIZE>::is_in_search_region, nb::arg("point"))
        .def("find_containing_bound", &BoundSet<double, Objective::MINIMIZE>::find_containing_bound, nb::arg("point"))
        .def_prop_ro("bounds", &BoundSet<double, Objective::MINIMIZE>::bounds);

    nb::class_<BoundSet<double, Objective::MAXIMIZE>>(m, "BoundSetMaximize")
        .def(nb::init<const std::vector<double>&>(), nb::arg("reference_point"))
        .def(nb::init<const std::vector<double>&, const std::vector<double>&>(), nb::arg("reference_point"), nb::arg("anti_reference"))
        .def("update_auto", &BoundSet<double, Objective::MAXIMIZE>::update_auto, nb::arg("point"))
        .def("update_re", &BoundSet<double, Objective::MAXIMIZE>::update_re, nb::arg("point"))
        .def("update_re_enhanced", &BoundSet<double, Objective::MAXIMIZE>::update_re_enhanced, nb::arg("point"))
        .def("update_ra_sa", &BoundSet<double, Objective::MAXIMIZE>::update_ra_sa, nb::arg("point"))
        .def("update_ra", &BoundSet<double, Objective::MAXIMIZE>::update_ra, nb::arg("point"))
        .def("update_naive", &BoundSet<double, Objective::MAXIMIZE>::update_naive, nb::arg("point"))
        .def("size", &BoundSet<double, Objective::MAXIMIZE>::size)
        .def("dimensions", &BoundSet<double, Objective::MAXIMIZE>::dimensions)
        .def("is_in_search_region", &BoundSet<double, Objective::MAXIMIZE>::is_in_search_region, nb::arg("point"))
        .def("find_containing_bound", &BoundSet<double, Objective::MAXIMIZE>::find_containing_bound, nb::arg("point"))
        .def_prop_ro("bounds", &BoundSet<double, Objective::MAXIMIZE>::bounds);

    nb::class_<NeighborhoodBoundSet<double, Objective::MINIMIZE>::AdjacencyGraph>(m, "AdjacencyGraphMinimize")
        .def_ro("nodes", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::AdjacencyGraph::nodes)
        .def_ro("adjacency_list", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::AdjacencyGraph::adjacency_list)
        .def_ro("k_neighbors", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::AdjacencyGraph::k_neighbors);

    nb::class_<NeighborhoodBoundSet<double, Objective::MAXIMIZE>::AdjacencyGraph>(m, "AdjacencyGraphMaximize")
        .def_ro("nodes", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::AdjacencyGraph::nodes)
        .def_ro("adjacency_list", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::AdjacencyGraph::adjacency_list)
        .def_ro("k_neighbors", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::AdjacencyGraph::k_neighbors);

    nb::class_<NeighborhoodBoundSet<double, Objective::MINIMIZE>>(m, "NeighborhoodBoundSetMinimize")
        .def(nb::init<const std::vector<double>&, const std::vector<double>&>(), nb::arg("reference_point"), nb::arg("anti_reference"))
        .def("update", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::update, nb::arg("point"))
        .def("size", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::size)
        .def("nonredundant_size", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::nonredundant_size)
        .def("dimensions", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::dimensions)
        .def("is_in_search_region", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::is_in_search_region, nb::arg("point"))
        .def("find_containing_bound", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::find_containing_bound, nb::arg("point"))
        .def("bounds", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::bounds)
        .def("nonredundant_bounds", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::nonredundant_bounds)
        .def("get_adjacency_graph", &NeighborhoodBoundSet<double, Objective::MINIMIZE>::get_adjacency_graph);

    nb::class_<NeighborhoodBoundSet<double, Objective::MAXIMIZE>>(m, "NeighborhoodBoundSetMaximize")
        .def(nb::init<const std::vector<double>&, const std::vector<double>&>(), nb::arg("reference_point"), nb::arg("anti_reference"))
        .def("update", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::update, nb::arg("point"))
        .def("size", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::size)
        .def("nonredundant_size", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::nonredundant_size)
        .def("dimensions", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::dimensions)
        .def("is_in_search_region", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::is_in_search_region, nb::arg("point"))
        .def("find_containing_bound", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::find_containing_bound, nb::arg("point"))
        .def("bounds", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::bounds)
        .def("nonredundant_bounds", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::nonredundant_bounds)
        .def("get_adjacency_graph", &NeighborhoodBoundSet<double, Objective::MAXIMIZE>::get_adjacency_graph);

    nb::class_<BoundSetTree<double, Objective::MINIMIZE>>(m, "BoundSetTreeMinimize")
        .def(nb::init<const std::vector<double>&, size_t, size_t>(), 
             nb::arg("reference_point"), nb::arg("max_leaf_size") = 32, nb::arg("num_children") = 8)
        .def(nb::init<const std::vector<double>&, const std::vector<double>&, size_t, size_t>(), 
             nb::arg("reference_point"), nb::arg("anti_reference"), nb::arg("max_leaf_size") = 32, nb::arg("num_children") = 8)
        .def("update_auto", &BoundSetTree<double, Objective::MINIMIZE>::update_auto, nb::arg("point"))
        .def("update_re", &BoundSetTree<double, Objective::MINIMIZE>::update_re, nb::arg("point"))
        .def("update_re_enhanced", &BoundSetTree<double, Objective::MINIMIZE>::update_re_enhanced, nb::arg("point"))
        .def("update_naive", &BoundSetTree<double, Objective::MINIMIZE>::update_naive, nb::arg("point"))
        .def("size", &BoundSetTree<double, Objective::MINIMIZE>::size)
        .def("dimensions", &BoundSetTree<double, Objective::MINIMIZE>::dimensions)
        .def("is_in_search_region", &BoundSetTree<double, Objective::MINIMIZE>::is_in_search_region, nb::arg("point"))
        .def("find_containing_bound", &BoundSetTree<double, Objective::MINIMIZE>::find_containing_bound, nb::arg("point"))
        .def("bounds", &BoundSetTree<double, Objective::MINIMIZE>::bounds);

    nb::class_<BoundSetTree<double, Objective::MAXIMIZE>>(m, "BoundSetTreeMaximize")
        .def(nb::init<const std::vector<double>&, size_t, size_t>(), 
             nb::arg("reference_point"), nb::arg("max_leaf_size") = 32, nb::arg("num_children") = 8)
        .def(nb::init<const std::vector<double>&, const std::vector<double>&, size_t, size_t>(), 
             nb::arg("reference_point"), nb::arg("anti_reference"), nb::arg("max_leaf_size") = 32, nb::arg("num_children") = 8)
        .def("update_auto", &BoundSetTree<double, Objective::MAXIMIZE>::update_auto, nb::arg("point"))
        .def("update_re", &BoundSetTree<double, Objective::MAXIMIZE>::update_re, nb::arg("point"))
        .def("update_re_enhanced", &BoundSetTree<double, Objective::MAXIMIZE>::update_re_enhanced, nb::arg("point"))
        .def("update_naive", &BoundSetTree<double, Objective::MAXIMIZE>::update_naive, nb::arg("point"))
        .def("size", &BoundSetTree<double, Objective::MAXIMIZE>::size)
        .def("dimensions", &BoundSetTree<double, Objective::MAXIMIZE>::dimensions)
        .def("is_in_search_region", &BoundSetTree<double, Objective::MAXIMIZE>::is_in_search_region, nb::arg("point"))
        .def("find_containing_bound", &BoundSetTree<double, Objective::MAXIMIZE>::find_containing_bound, nb::arg("point"))
        .def("bounds", &BoundSetTree<double, Objective::MAXIMIZE>::bounds);
}
