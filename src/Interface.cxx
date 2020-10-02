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
        .def("defineBody", &SimManager::defineBody,
            py::arg("body"))
        .def("defineInitialState", &SimManager::defineInitialState,
        	py::arg("U0"),
        	py::arg("V0"),
        	py::arg("T0"))
        .def("defineUniformDynamicBoundaryConditions", &SimManager::defineUniformDynamicBoundaryConditions,
            py::arg("U_bc"),
            py::arg("V_bc"),
            py::arg("P_bc"))
        .def("defineDynamicBoundaryConditions", &SimManager::defineDynamicBoundaryConditions,
            py::arg("U_bc"),
            py::arg("V_bc"),
            py::arg("P_bc"))
        .def("defineThermalBoundaryConditions", &SimManager::defineThermalBoundaryConditions,
        	py::arg("type"),
        	py::arg("value"))
        .def("defineThermalIntegrationParameters", &SimManager::defineThermalIntegrationParameters,
            py::arg("theta"),
            py::arg("accuracy"))
        .def("defineFlowIntegrationParameters", &SimManager::defineFlowIntegrationParameters,
            py::arg("theta"),
            py::arg("accuracy_u"),
            py::arg("accuracy_v"),
            py::arg("accuracy_p"),
            py::arg("clean_pressure"))
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