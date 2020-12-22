#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "LBMKernel.hxx"
#include "CollisionOperator.hxx"
#include "StreamingOperator.hxx"
#include "BoundaryCondition.hxx"

namespace py = pybind11;

void interface_lbm(py::module& m)
{
	py::class_<LBMKernel>(m, "LBM")
        .def(py::init<>())
        .def("buildSpace", &LBMKernel::build_space,
            py::arg("x_min"),
            py::arg("x_max"),
            py::arg("nx"),
            py::arg("y_min"),
            py::arg("y_max"),
            py::arg("ny"))
        .def("buildTime", &LBMKernel::build_time,
            py::arg("t_start"),
            py::arg("t_end"),
            py::arg("nb_steps"))
        .def("setInitialState", &LBMKernel::set_initial_state,
            py::arg("rho"),
            py::arg("u"),
            py::arg("e"),
            py::arg("gamma"))
        .def("simulate", &LBMKernel::simulate,
            py::arg("collision_op"),
            py::arg("streaming_op"),
            py::arg("f_out"));

    py::class_<CollisionOperator>(m, "CollisionOperator");

    py::class_<BGK, CollisionOperator>(m, "BGK")
        .def(py::init<double, double>(),
            py::arg("nu"),
            py::arg("gamma"));

    py::class_<StreamingOperator>(m, "StreamingOperator");

    py::class_<Stream2D, StreamingOperator>(m, "Stream2D")
        .def(py::init<std::vector<BoundaryCondition>>(),
            py::arg("boundary_conditions"));

    py::class_<BoundaryCondition>(m, "BoundaryCondition");

    py::class_<D2Q9BounceBack, BoundaryCondition>(m, "D2Q9BounceBack")
        .def(py::init<std::vector<std::pair<size_t, size_t>>, std::vector<std::vector<size_t>>>(),
            py::arg("ij"),
            py::arg("dir"));
}