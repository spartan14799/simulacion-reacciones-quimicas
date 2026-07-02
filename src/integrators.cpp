
#include "integrators.h"

//explicit

void integrate_euler(deriv_t fderiv, system_t & s, double tinit, double tend, 
                    double dt, printer_t writer)
{
    // vector to store derivs
    system_t dsdt(s.size());

    // time loop
    for(double t = tinit; t <= tend; t = t + dt) { // NOTE: Last time step not necessarily tf
        
        // write state
        writer(s, t);
        
        // compute derivs
        fderiv(s, dsdt, t); //fderiv dsdt=smthng s state vector

        // compute new state. NOTE: Not using components, assuming valarray or similar
        s = s + dt*dsdt;

      }
}

void integrate_heun(deriv_t fderiv, system_t & s, double tinit, double tend, 
                    double dt, printer_t writer) //rk2
{
    // vector to store derivs
    system_t dsdt(s.size());
    system_t k1(s.size());
    system_t k2(s.size());

    // time loop
    for(double t = tinit; t <= tend; t = t + dt) { // NOTE: Last time step not necessarily tf
        
        writer(s, t);

        fderiv(s, dsdt, t);
        k1 = dsdt;
        fderiv(s + dt*k1, dsdt, t+dt);
        k2 = dsdt;

        s = s + 0.5*dt*(k1 + k2);
      }
}

void integrate_rk4(deriv_t fderiv, system_t & s, double tinit, double tend, 
                    double dt, printer_t writer)
{
    // vector to store derivs
    system_t dsdt(s.size());
    system_t k1(s.size());
    system_t k2(s.size());
    system_t k3(s.size());
    system_t k4(s.size());

    // time loop
    for(double t = tinit; t <= tend; t = t + dt) { // NOTE: Last time step not necessarily tf
        
        writer(s, t);
        
        // compute derivs
        fderiv(s, dsdt, t);
        k1 = dsdt;
        fderiv(s + dt*k1/2, dsdt, t+dt/2);
        k2 = dsdt;
        fderiv(s+dt*k2/2, dsdt, t+dt/2);
        k3=dsdt;
        fderiv(s+dt*k3, dsdt, t+dt);
        k4=dsdt;

        s = s + (dt/6.0)*(k1 + 2*k2+2*k3+k4);
      }
}



//implicit



void integrate_implicit_euler(deriv_t fderiv, system_t & s, double tinit, double tend, 
                              double dt, printer_t writer)
{
    int dim = s.size();
    system_t dsdt(dim); //deriv en el tiempo actual para predictor de euler
    system_t s_next(dim);
    system_t s_guess(dim);
    system_t f_next(dim); 

    const int max_iter = 50;
    const double tolerancia = 1e-9;

    //buucle temporal
    for(double t = tinit; t <= tend; t += dt) {
        
        writer(s, t);
        fderiv(s, dsdt, t);

        s_guess = s + dt * dsdt; //primera predict usando euler explicito

        // implicit
        for(int ii = 0; ii < max_iter; ii++) { //iteración de punto fijo

            fderiv(s_guess, f_next, t + dt); //evalua la derivada usando la suposicion s_guess
            s_next = s + dt * f_next;  //fórmula de euler implicito

            double error = 0.0;
            for(int jj = 0; jj < dim; ++jj) {
                error += std::pow(s_next[jj] - s_guess[jj], 2);
            }
            error = std::sqrt(error); //distancia entre suposicion y resultado

            if (error < tolerancia) { //si convergió nos salimos
                break;
            }
            s_guess = s_next; //si no convergio actualizamos y repetimos
        }
        
        s = s_next;
    }
}

void integrate_crank_nicolson(deriv_t fderiv, system_t & s, double tinit, double tend, 
                              double dt, printer_t writer)
{
    int dim = s.size();
    system_t dsdt(dim);
    system_t s_next(dim);
    system_t s_guess(dim);
    system_t f_next(dim);

    const int max_iter = 50;
    const double tolerancia = 1e-9;

    // Bucle temporal
    for(double t = tinit; t <= tend; t += dt) {
        
        writer(s, t);
        fderiv(s, dsdt, t);
        //predict with euler
        s_guess = s + dt * dsdt;

        //implicit
        for(int ii = 0; ii < max_iter; ii++) {
            fderiv(s_guess, f_next, t + dt); //guarda en fnext

            s_next = s + 0.5 * dt * (dsdt + f_next);

            //cambio respecto a anterior iter despreciable?
            double error = 0.0;
            for(int jj = 0; jj < dim; ++jj) {
                error += std::pow(s_next[jj] - s_guess[jj], 2);
            }
            error = std::sqrt(error);

            //si sí nos salimos xd
            if (error < tolerancia) {
                break;
            }

            //sino, la corrección es la nueva suposición
            s_guess = s_next;
        }
        s = s_next;
    }
}
