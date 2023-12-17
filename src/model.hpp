//
//  model.hpp
//
//  Created by Gia-Wei Chern on 3/3/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#ifndef lattice_hpp
#define lattice_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <random>

#include "vec3.hpp"
#include "util.hpp"

using namespace std;

typedef mt19937 RNG;

class Square_lattice {
public:
    int L, Ns;
    int dim;
    
    double t1;
    double t2;
    double U;
    
    double filling;
    double mu;
    double kT;

    int solver_type;
    
    arma::mat Mpp, Npp, Upp;
    
    static constexpr int N_nn1 = 4;
    static constexpr int N_nn2 = 4;
    class Site {
    public:
        int idx;
        int x, y;
        
        Site * nn1[N_nn1];
        Site * nn2[N_nn2];        

        cx_double R;
        cx_double Delta;
        
        double lambda;
        
        
        Site(void) { };
        
    } *site;
    
    double *onsite_V;
    
    arma::sp_cx_mat Hamiltonian;
    arma::cx_mat Density_Mat;
    arma::cx_mat Phi;

    Square_lattice(int linear_size) {
      
        L = linear_size;
        Ns = L * L;
        dim = Ns;
        
        onsite_V = new double[Ns];
        
        site = new Site[Ns];
        
        Hamiltonian = arma::sp_cx_mat(dim, dim);
        Density_Mat = arma::cx_mat(dim, dim);

        Phi = arma::cx_mat(3, dim);
        
        init_lattice();
        init_pp_matrices();

        rng = RNG(seed());

        velocity = vector<double>(Ns, 0);
        force = vector<double>(Ns, 0);

    };
    
    ~Square_lattice(void) {
        
        delete [] site;
    };
    
    inline int index(int x, int y) {
        return L * y + x;
    };
    
    void init_lattice(void);
    void init_quenched_disorder(double W);

    void init_pp_matrices(void);

    
    void reset_GA_parameters(void) {
        for(int i=0; i<Ns; i++) {
            site[i].R = 1.0;
            site[i].lambda = U;
        }
    };

    double error;
    double err_f;
    int iter;
    int max_iter;


    
    arma::sp_cx_mat build_Hamiltonian(void);
    //arma::sp_cx_mat build_Hamiltonian(arma::cx_mat &, arma::cx_mat &);
    
    void compute_fermi_level(arma::vec & eigE);
    arma::cx_mat compute_density_matrix(arma::sp_cx_mat &);

    arma::cx_mat compute_tV_density_matrix(arma::sp_cx_mat &);
    int compute_GA_solution(void);
    
    void compute_Deltas(arma::cx_mat &);
    void compute_renormalizations(arma::cx_mat &, arma::cx_mat &);
    arma::cx_mat solve_Phi(arma::cx_mat &);
    double adjust_lambda(arma::cx_mat &, arma::cx_mat &);
    
    arma::cx_mat init_uncorrelated_Phi(arma::cx_mat &);
    
    void statistics_GA(int ndata, double U, double W);


    RNG rng;
    std::random_device seed;

    //======= MD =======
    void move_atoms();
    void compute_forces();
    void integrate_forces();
    void integrate_Langevin();
    void step_NVT();

    void print_config(string const filename);

    double g;
    double k0;
    double kappa;
    double mass;
    double gamma;
    double dt;
    
    vector<double> velocity;
    vector<double> force;

    
    
    // ========================================
    
    double GA_kT_start;
    double GA_kT_stop;
    int GA_n_annealing;
    double GA_r_annealing;
    
    double avg_R2 = 0, sgm_R2 = 0;
    double avg_nd = 0, sgm_nd = 0;
    double avg_np = 0, sgm_np = 0;
    double avg_d = 0, sgm_d = 0;
    
    double E_e = 0;
    
    void compute_electronic_energy(void);

    void analyze_data(void);
    
    void test(void);
};


#endif /* lattice_hpp */




