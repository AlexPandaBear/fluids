#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Field.hxx"

namespace py = pybind11;

void interface_common(py::module& m)
{
	py::class_<Field>(m, "Field", py::buffer_protocol())
        .def(py::init<std::vector<size_t>>(),
            py::arg("dimensions"))
        .def("reshape", &Field::reshape,
            py::arg("dimensions"))
        .def("__getitem__", [](Field& f, std::vector<size_t> index){return f(index);})
        .def("__setitem__", [](Field& f, std::vector<size_t> index, double value){f(index) = value;})
        .def_buffer([](Field& f) -> py::buffer_info
			{
				return py::buffer_info(
					f.get_ptr(),
					f.get_element_size(),
					py::format_descriptor<double>::format(),
					f.get_nb_dimensions(),
					py::detail::any_container<long int>(f.get_dimensions_sizes()),
					py::detail::any_container<long int>(f.get_adjusted_element_sizes())
				);
			});
}