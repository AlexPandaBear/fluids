/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-std=c++11', '-I/usr/local/lib/python3.7/dist-packages/pybind11-2.4.3-py3.7.egg/']
cfg['dependencies'] = ['SimManager.hxx', 'IncompressibleKernel.hxx', 'ThermalKernel.hxx']
cfg['sources'] = ['SimManager.cxx', 'IncompressibleKernel.cxx', 'ThermalKernel.cxx']
%>
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "SimManager.hxx"

namespace py = pybind11;

PYBIND11_MODULE(_vf, m)
{
    py::class_<SimManager>(m, "SM")
    	.def(py::init<>())
        .def("defineTimeParameters", &SimManager::defineTimeParameters,
            py::arg("tMax"),
            py::arg("nb_steps"))
        .def("defineGridParameters", &SimManager::defineGridParameters,
            py::arg("Lx"),
            py::arg("Ly"),
            py::arg("nx"),
            py::arg("ny"))
        .def("defineFluidProperties", &SimManager::defineFluidProperties,
            py::arg("lambda"),
            py::arg("rho"),
            py::arg("cv"),
            py::arg("mu"))
        .def("defineInitialState", &SimManager::defineInitialState,
        	py::arg("U0"),
        	py::arg("V0"),
        	py::arg("P0"),
        	py::arg("T0"))
        .def("defineBoundaryConditions", &SimManager::defineBoundaryConditions,
        	py::arg("type"),
        	py::arg("value"))
        .def("defineErrorTolerance", &SimManager::defineErrorTolerance,
        	py::arg("max_pressure_error"))
        .def("launchSimulation", &SimManager::launchSimulation)
        .def("getTemperatureAt", &SimManager::getTemperatureAt,
        	py::arg("t"),
        	py::arg("i"),
        	py::arg("j"))
        .def("getPressureAt", &SimManager::getPressureAt,
        	py::arg("t"),
        	py::arg("i"),
        	py::arg("j"))
        .def("getXVelocityAt", &SimManager::getXVelocityAt,
        	py::arg("t"),
        	py::arg("i"),
        	py::arg("j"))
        .def("getYVelocityAt", &SimManager::getYVelocityAt,
        	py::arg("t"),
        	py::arg("i"),
        	py::arg("j"))
        .def("saveData", &SimManager::saveData,
        	py::arg("file_name"))
        .def("loadData", &SimManager::loadData,
        	py::arg("file_name"));
}