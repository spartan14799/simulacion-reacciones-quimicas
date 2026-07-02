//API
#pragma once

#include <iostream>
#include <valarray>
#include <functional>

using system_t = std::valarray<double>;
using deriv_t = std::function<void(const system_t&, system_t&, double)>; //pa decirle que es void function
using printer_t = std::function<void(const system_t&, double)>;


//explicit
void integrate_euler(deriv_t fderiv, system_t & s, double tinit, double tend, 
                    double dt, printer_t writer);

void integrate_heun(deriv_t fderiv, system_t & s, double tinit, double tend, 
                    double dt, printer_t writer); //rk2

void integrate_rk4(deriv_t fderiv, system_t & s, double tinit, double tend, 
                    double dt, printer_t writer);


//implicit
void integrate_implicit_euler (deriv_t fderiv, system_t&s, double tinit, double tend,
                                double dt, printer_t writer);

void integrate_crank_nicolson(deriv_t fderiv, system_t & s, double tinit, double tend, 
                              double dt, printer_t writer);

