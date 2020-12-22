#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "SimManager.hxx"

namespace py = pybind11;


void interface_lbm(py::module&);
void interface_common(py::module&);


PYBIND11_MODULE(_fluids, m)
{
    interface_lbm(m);
    interface_common(m);

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
        	py::arg("file_name"))
        .def("getTemperatureFieldAt", &SimManager::getTemperatureFieldAt,
            py::arg("t"))
        .def("getPressureFieldAt", &SimManager::getPressureFieldAt,
            py::arg("t"))
        .def("getXVelocityFieldAt", &SimManager::getXVelocityFieldAt,
            py::arg("t"))        
        .def("getYVelocityFieldAt", &SimManager::getYVelocityFieldAt,
            py::arg("t"))
        .def("getVelocityNormFieldAt", &SimManager::getVelocityNormFieldAt,
            py::arg("t"))
        .def("getVelocityDivergenceFieldAt", &SimManager::getVelocityDivergenceFieldAt,
            py::arg("t"));  
}