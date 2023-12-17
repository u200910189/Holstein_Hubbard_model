//
//  ga.cpp
//
//  Created by Gia-Wei Chern on 6/20/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include <iostream>

#include "model.hpp"
#include "util.hpp"
#include "test.hpp"

void HH_simulation(Square_lattice & system, double W, int steps, int wait_steps, int par_steps) {

    int n_save = par_steps;

    system.init_quenched_disorder(W);

    for (int r=0; r < steps; r++) {

        if (r % n_save == 0)
            cout << "r = " << r << endl;
        
        system.step_NVT();

        if (r % n_save == 0 && r >= wait_steps) {
            string f_name = "c" + std::to_string(r) + ".dat";
            system.print_config(f_name);
        }
    }

}


int main(int argc, const char * argv[]) {
    
    std::cout << "Hello, World!\n";
    
    double U                = argc > 1 ? atof(argv[1]) : 1.;
    double W                = argc > 2 ? atof(argv[2]) : 1.1;
    int L                   = argc > 3 ? atoi(argv[3]) : 10;
    double filling_frac     = argc > 4 ? atof(argv[4]) : 0.5;
    
    double kT_start         = argc > 5 ? atof(argv[5]) : 0.15;
    double kT_stop          = argc > 6 ? atof(argv[6]) : 0.1;
    int n_annealing         = argc > 7 ? atoi(argv[7]) : 25;
    double r_annealing      = argc > 8 ? atof(argv[8]) : 0.95;

    int steps = argc > 9 ? atoi(argv[9]) : 40;
    int wait_steps = argc > 10 ? atoi(argv[10]) : 0;
   
    double g = argc > 11 ? atof(argv[11]) : 1.5;
    double k0 = argc > 12 ? atof(argv[12]) : 1.;
    double kappa = argc > 13 ? atof(argv[13]) : 0;  //0.18;
    double t2 = argc > 14 ? atof(argv[14]) : 0;

    int solver_type = argc > 15 ? atoi(argv[15]) : 0;
    int par_steps = argc > 16 ? atoi(argv[16]) : 20; 


    cout << "L = " << L << endl;
    cout << "U = " << U << "\t W = " << W << endl;
    
    //test_all();
    
    Square_lattice system(L);
    
    system.GA_kT_start = kT_start;
    system.GA_kT_stop = kT_stop;
    system.GA_n_annealing = n_annealing;
    system.GA_r_annealing = r_annealing;
    system.max_iter = 2000;

    system.t1 = -1;
    system.t2 = t2;
    system.filling = filling_frac;

    system.U = U;
    system.g = g;
    system.k0 = k0;
    system.kappa = kappa;
    system.mass = 5.;
    system.gamma = 0.04;
    system.dt = 0.05;

    system.solver_type = solver_type;

    //double n_data = 100;
    //system.statistics_GA(n_data, U, W);

    HH_simulation(system, W, steps, wait_steps, par_steps);



    return 0;
}
