// wrapper.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "./integrators.h"

namespace py = pybind11;


// modulo de python llamado "metodos_numericos"
PYBIND11_MODULE(metodos_numericos, m) {


    m.doc() = "lib interna de integradores numericos en cpp";

    //explicitos

    // exponer euler explicito
    m.def("integrate_euler", &integrate_euler, 
          "integra usando euler explicito",
          py::arg("fderiv"), py::arg("s"), py::arg("tinit"), py::arg("tend"), py::arg("dt"), py::arg("writer"));

    // exponer heun
    m.def("integrate_heun", &integrate_heun, 
          "integra usando heun",
          py::arg("fderiv"), py::arg("s"), py::arg("tinit"), py::arg("tend"), py::arg("dt"), py::arg("writer"));

    // exponer rk4
    m.def("integrate_rk4", &integrate_rk4, 
          "integra usando rk4",
          py::arg("fderiv"), py::arg("s"), py::arg("tinit"), py::arg("tend"), py::arg("dt"), py::arg("writer"));


    //implicitos

    // exponer euler implicito
    m.def("integrate_implicit_euler", &integrate_implicit_euler, 
          "integra usando euler implicito",
          py::arg("fderiv"), py::arg("s"), py::arg("tinit"), py::arg("tend"), py::arg("dt"), py::arg("writer"));

    // exponer cranck nickolson
    m.def("integrate_crank_nicolson", &integrate_crank_nicolson, 
          "integra usando crank-nicolson",
          py::arg("fderiv"), py::arg("s"), py::arg("tinit"), py::arg("tend"), py::arg("dt"), py::arg("writer"));
}