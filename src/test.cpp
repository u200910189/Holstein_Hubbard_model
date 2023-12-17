//
//  test.cpp
//
//  Created by Gia-Wei Chern on 3/3/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include <iostream>
#include <fstream>

#include "test.hpp"
#include "model.hpp"
#include "util.hpp"

using namespace std;

void test_all(void) {
    
    //test_sparse_matrix();
    
    //test_lattice();
    
    test_Hamiltonian();
    
    exit(1);
}

void test_lattice(void) {
    
}

void test_sparse_matrix(void) {
    int dim = 40000;
    arma::sp_cx_mat A(dim, dim);
    arma::sp_cx_mat B(dim, dim);
    arma::sp_cx_mat C(dim, dim);
    
    std::random_device seed;
    
    RNG rng = RNG(seed());
    
    std::uniform_real_distribution<double> rd;
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    
    for(int i=0; i<50000; i++) {
        int m = (int) (dim * rd(rng));
        int n = (int) (dim * rd(rng));
        A(m, n) = rn(rng);
    }
    
    
    for(int i=0; i<dim-1; i++) B(i, i+1) = 0.23 * i * (i-1);
    
    C = A * B;
    
    cout << C << endl;
}

void test_square_model(double W, double U0) {

    Square_lattice system(30);
    
    system.t1 = -1;
    system.U = U0;
    
    system.kT = 0.200001;
    system.filling = 0.5;
    
    
    system.init_quenched_disorder(W);
    
    system.compute_GA_solution();
    
}

void test_Hamiltonian(void) {
    
    Square_lattice system(4);
    
    system.t1 = -1;
    system.U = 3.0;
    
    system.kT = 0.200001;
    system.filling = 0.5;
    
    system.init_quenched_disorder(0);

    for(int i=0; i<system.Ns; i++) system.site[i].R = 0.5;
    auto Hm = system.build_Hamiltonian();
    auto Dm = system.compute_tV_density_matrix(Hm);
    
    arma::vec eigval(system.dim);
    arma::cx_mat eigvec;
    arma::cx_mat Hd(Hm);
    
    arma::eig_sym(eigval, eigvec, Hd);

    cout << Hm << endl;
    
    cout << eigval << endl;
    cout << eigvec << endl;
    
    int ix = 1;
    for(int k=0; k<system.N_nn1; k++) {
        int j = system.site[ix].nn1[k]->idx;
        cout << Dm(ix, j) << endl;
    }
    
    
    
}
