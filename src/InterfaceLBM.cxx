#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include "LBMKernel.hxx"
#include "LBMConfig.hxx"

namespace py = pybind11;

void interface_lbm(py::module& m)
{
	py::class_<LBMKernel>(m, "LBM")
        .def(py::init<LBMConfig>(),
            py::arg("config"))
        .def(py::init<std::string>(),
            py::arg("config_file"))
        .def("setInitialState", &LBMKernel::set_initial_state,
            py::arg("rho"),
            py::arg("u"),
            py::arg("e"))
        .def("simulate", &LBMKernel::simulate,
            py::arg("f_out"));

    py::class_<LBMConfig>(m, "LBMConfig")
        .def(py::init<>())
        .def(py::init<std::string>(),
            py::arg("config_file"));
}