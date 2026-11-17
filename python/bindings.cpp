#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include "local_bounds.hpp"

namespace nb = nanobind;
using namespace local_bounds;

// The C++ classes take the objective sense as a template parameter, so each is
// bound once per sense. python/local_bounds/__init__.py picks one by `sense`.
template <Objective Sense>
void bind_sense(nb::module_& m, const char* set_name, const char* tree_name,
                const char* neighborhood_name, const char* graph_name) {
    using Set = BoundSet<double, Sense>;
    using Tree = BoundSetTree<double, Sense>;
    using Neighborhood = NeighborhoodBoundSet<double, Sense>;
    using Graph = typename Neighborhood::AdjacencyGraph;
    using Coordinates = const std::vector<double>&;

    nb::class_<Set>(m, set_name)
        .def(nb::init<Coordinates>(), nb::arg("reference_point"))
        .def(nb::init<Coordinates, Coordinates>(), nb::arg("reference_point"), nb::arg("anti_reference"))
        .def("update_auto", &Set::update_auto, nb::arg("point"))
        .def("update_re", &Set::update_re, nb::arg("point"))
        .def("update_re_enhanced", &Set::update_re_enhanced, nb::arg("point"))
        .def("update_ra_sa", &Set::update_ra_sa, nb::arg("point"))
        .def("update_ra", &Set::update_ra, nb::arg("point"))
        .def("update_naive", &Set::update_naive, nb::arg("point"))
        .def("size", &Set::size)
        .def("dimensions", &Set::dimensions)
        .def("is_in_search_region", &Set::is_in_search_region, nb::arg("point"))
        .def("find_containing_bound", &Set::find_containing_bound, nb::arg("point"))
        .def("bounds", &Set::bounds)
        .def_prop_ro("sense", [](const Set&) { return Sense; });

    nb::class_<Graph>(m, graph_name)
        .def_ro("nodes", &Graph::nodes)
        .def_ro("adjacency_list", &Graph::adjacency_list)
        .def_ro("k_neighbors", &Graph::k_neighbors)
        .def_ro("quasi", &Graph::quasi);

    nb::class_<Neighborhood>(m, neighborhood_name)
        .def(nb::init<Coordinates, Coordinates>(), nb::arg("reference_point"), nb::arg("anti_reference"))
        .def("update", &Neighborhood::update, nb::arg("point"))
        .def("size", &Neighborhood::size)
        .def("nonredundant_size", &Neighborhood::nonredundant_size)
        .def("dimensions", &Neighborhood::dimensions)
        .def("is_in_search_region", &Neighborhood::is_in_search_region, nb::arg("point"))
        .def("find_containing_bound", &Neighborhood::find_containing_bound, nb::arg("point"))
        .def("bounds", &Neighborhood::bounds)
        .def("nonredundant_bounds", &Neighborhood::nonredundant_bounds)
        .def("get_adjacency_graph", &Neighborhood::get_adjacency_graph, nb::arg("include_quasi") = false)
        .def_prop_ro("sense", [](const Neighborhood&) { return Sense; });

    nb::class_<Tree>(m, tree_name)
        .def(nb::init<Coordinates, size_t, size_t>(),
             nb::arg("reference_point"), nb::arg("max_leaf_size") = 32, nb::arg("num_children") = 8)
        .def(nb::init<Coordinates, Coordinates, size_t, size_t>(),
             nb::arg("reference_point"), nb::arg("anti_reference"), nb::arg("max_leaf_size") = 32, nb::arg("num_children") = 8)
        .def("update_auto", &Tree::update_auto, nb::arg("point"))
        .def("update_re", &Tree::update_re, nb::arg("point"))
        .def("update_re_enhanced", &Tree::update_re_enhanced, nb::arg("point"))
        .def("update_naive", &Tree::update_naive, nb::arg("point"))
        .def("size", &Tree::size)
        .def("dimensions", &Tree::dimensions)
        .def("is_in_search_region", &Tree::is_in_search_region, nb::arg("point"))
        .def("find_containing_bound", &Tree::find_containing_bound, nb::arg("point"))
        .def("bounds", &Tree::bounds)
        .def_prop_ro("sense", [](const Tree&) { return Sense; });
}

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

    bind_sense<Objective::MINIMIZE>(m, "BoundSetMinimize", "BoundSetTreeMinimize",
                                    "NeighborhoodBoundSetMinimize", "AdjacencyGraphMinimize");
    bind_sense<Objective::MAXIMIZE>(m, "BoundSetMaximize", "BoundSetTreeMaximize",
                                    "NeighborhoodBoundSetMaximize", "AdjacencyGraphMaximize");
}
