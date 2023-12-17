//
//  model.cpp
//
//  Created by Gia-Wei Chern on 3/3/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include <cassert>
#include <iomanip>
#include "model.hpp"

using namespace std;

// {s1, s2} components of pauli matrix vector,
// sigma1     sigma2     sigma3
//  0  1       0 -I       1  0
//  1  0       I  0       0 -1
const Vec3<cx_double> pauli[2][2] {
    {{0, 0, 1}, {1, -_I, 0}},
    {{1, _I, 0}, {0, 0, -1}}
};


void Square_lattice::init_pp_matrices(void) {
    Mpp = Npp = Upp = arma::mat(3, 3);
    Mpp.zeros();
    Npp.zeros();
    Upp.zeros();
    
    Mpp(0, 1) = Mpp(1, 0) = Mpp(1, 2) = Mpp(2, 1) = sqrt(2.);
    Npp(1, 1) = 0.5;
    Npp(2, 2) = 1;
    Upp(2, 2) = 1;
}


void Square_lattice::init_lattice(void) {
    
    for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
            int idx = index(x, y);
            site[idx].idx = idx;
            site[idx].x = x;
            site[idx].y = y;
    }
    
    for(int i=0; i<Ns; i++) {
        int j;
        int x = site[i].x;
        int y = site[i].y;
        
        j = index(mod(x+1, L), y);
        site[i].nn1[0] = &site[j];
        
        j = index(mod(x-1, L), y);
        site[i].nn1[1] = &site[j];
        
        j = index(x, mod(y+1, L));
        site[i].nn1[2] = &site[j];
        
        j = index(x, mod(y-1, L));
        site[i].nn1[3] = &site[j];



        j = index(mod(x+1, L), mod(y+1, L));
        site[i].nn2[0] = &site[j];
        
        j = index(mod(x-1, L), mod(y+1, L));
        site[i].nn2[1] = &site[j];
        
        j = index(mod(x+1, L), mod(y-1, L));
        site[i].nn2[2] = &site[j];
        
        j = index(mod(x-1, L), mod(y-1, L));
        site[i].nn2[3] = &site[j];
    }
}

void Square_lattice::init_quenched_disorder(double W) {
    
    std::uniform_real_distribution<double> rd(-W, W);

    for (int i=0; i<Ns; i++) 
        onsite_V[i] = rd(rng);
}

arma::sp_cx_mat Square_lattice::build_Hamiltonian(void) {
    
    arma::sp_cx_mat H(dim, dim);
    H.zeros();
    
    for (int i=0; i<Ns; i++) {
        for (int k=0; k<N_nn1; k++) {
            int j = site[i].nn1[k]->idx;
            
            H(i, j) = site[i].R * conj(site[j].R) * t1;
        }

        for (int k=0; k<N_nn2; k++) {
            int j = site[i].nn2[k]->idx;
            
            H(i, j) = site[i].R * conj(site[j].R) * t2;
        }
        
        H(i, i) = -g * onsite_V[i] + site[i].lambda;
    }
    
    return H;
}


//arma::sp_cx_mat Square_lattice::build_Hamiltonian(arma::cx_mat & Dm, arma::cx_mat & f) {
//    compute_renormalizations(Dm, f);
//    arma::sp_cx_mat H = build_Hamiltonian();
//    return H;
//}

void Square_lattice::compute_renormalizations(arma::cx_mat & Dm, arma::cx_mat & f) {
    
    for(int i=0; i<Ns; i++) {
        double n0 = real(Dm(i, i));
        //double n0 = real(conj(f(1, i)) * f(1, i) + conj(f(2, i)) * f(2, i));

        site[i].R = (conj(f(2, i)) * f(1, i) + conj(f(1, i)) * f(0, i)) / sqrt(n0 * (1. - n0));
  
//        double mix = 0.5;
//        site[i].R = mix * site[i].R + (1.-mix) * (conj(f(2, i)) * f(1, i) + conj(f(1, i)) * f(0, i)) / sqrt(n0 * (1. - n0));

    }
}

void Square_lattice::compute_Deltas(arma::cx_mat & Dm) {
    for (int i=0; i<Ns; i++) {
        site[i].Delta = 0;
        for (int k=0; k<N_nn1; k++) {
            int j = site[i].nn1[k]->idx;
            site[i].Delta += t1 * conj(site[j].R) * Dm(j, i);
        }

        for (int k=0; k<N_nn2; k++) {
            int j = site[i].nn2[k]->idx;
            site[i].Delta += t2 * conj(site[j].R) * Dm(j, i);
        }
    }
}

void Square_lattice::compute_fermi_level(arma::vec & eigE) {
    
    double x1 = eigE(0);
    double x2 = eigE(eigE.size()-1);
    
    int max_bisection = 200;
    double eps_bisection = 1.e-12;
    double eps_conv = 1.e-6;
    
    int iter = 0;
    while(iter < max_bisection || fabs(x2 - x1) > eps_bisection) {
        
        double xm = 0.5*(x1 + x2);
        double density = 0;
        for(unsigned int i=0; i<eigE.size(); i++) {
            density += fermi_density(eigE(i), kT, xm);
        }
        density /= ((double) dim);
         
        if (fabs(density - filling) < eps_conv)
            break;
   
        if (density <= filling) 
            x1 = xm;
        else 
            x2 = xm;
        
        iter++;

        //cout << iter << " " << x1 << " " << x2 << " " << fabs(x2 - x1) << " " << setprecision(15) << xm << " " << density << endl;
    }
   
    mu = 0.5*(x1 + x2);
    
}

arma::cx_mat Square_lattice::compute_density_matrix(arma::sp_cx_mat & H) {
    
    arma::cx_mat Hm(H);
    arma::cx_mat D(dim, dim);
    
    arma::vec eigval(dim);
    arma::cx_mat eigvec;
    
    arma::eig_sym(eigval, eigvec, Hm);
    
    compute_fermi_level(eigval);
    
    arma::vec fd_factor(dim);
    for(int i=0; i<dim; i++) 
        fd_factor(i) = fermi_density(eigval(i), kT, mu);

    D.zeros();
    for(int a=0; a<dim; a++)
    for(int b=a; b<dim; b++) {
    
        cx_double sum = 0;
        for(int m=0; m<dim; m++) {
            sum += fd_factor(m) * conj(eigvec(a, m)) * eigvec(b, m);
        }
        D(a, b) = sum;
        if(a != b) D(b, a) = conj(sum);
    }
    
    return D;
}

arma::cx_mat Square_lattice::compute_tV_density_matrix(arma::sp_cx_mat & H) {
    
    arma::cx_mat Hm(H);
    arma::cx_mat D(dim, dim);

    arma::vec eigval(dim);
    arma::cx_mat eigvec;
    
    arma::eig_sym(eigval, eigvec, Hm);
    
    compute_fermi_level(eigval);
    
    arma::vec fd_factor(dim);
    for(int i=0; i<dim; i++) 
        fd_factor(i) = fermi_density(eigval(i), kT, mu);
        

    D.zeros();
    for (int i=0; i<Ns; i++) {
        cx_double sum = 0;
        for (int m=0; m<dim; m++) {
            sum += fd_factor(m) * conj(eigvec(i, m)) * eigvec(i, m);
        }
        D(i, i) = sum;

        for (int k=0; k<N_nn1; k++) {
            int j = site[i].nn1[k]->idx;
            
            if (j > i) {
                cx_double sum2 = 0;
                for (int m=0; m<dim; m++) {
                    sum2 += fd_factor(m) * conj(eigvec(j, m)) * eigvec(i, m);
                }
                D(i, j) = sum2;
                D(j, i) = conj(sum2);
            }
        }

        for (int k=0; k<N_nn2; k++) {
            int j = site[i].nn2[k]->idx;
            
            if (j > i) {
                cx_double sum2 = 0;
                for (int m=0; m<dim; m++) {
                    sum2 += fd_factor(m) * conj(eigvec(j, m)) * eigvec(i, m);
                }
                D(i, j) = sum2;
                D(j, i) = conj(sum2);
            }
        }
    }

    return D;
}

arma::cx_mat Square_lattice::init_uncorrelated_Phi(arma::cx_mat & Dm) {
    
    arma::cx_mat f(3, dim);
    for(int i=0; i<Ns; i++) {
//        double n0 = real(Dm(i, i));
//        f(0, i) = 1. - n0;
//        f(1, i) = sqrt(n0 * (1. - n0));
//        f(2, i) = n0;
        
        f(0, i) = 0.5;
        f(1, i) = 1./sqrt(2.);
        f(2, i) = 0.5;
        
    }
    
    return f;
}

arma::cx_mat Square_lattice::solve_Phi(arma::cx_mat & Dm) {
    
    //compute_renormalizations(Dm, Phi);
    
    compute_Deltas(Dm);
    
    arma::cx_mat f(3, dim);
    
    for(int i=0; i<Ns; i++) {
        arma::cx_mat Hpp(3, 3);
        double n0 = real(Dm(i, i));
        //double n0 = real(conj(f(1, i)) * f(1, i) + conj(f(2, i)) * f(2, i));
        
        auto tmpx = site[i].R * site[i].Delta;
        double Lmd = 2. * real(tmpx) * (n0 - 0.5) / (n0 * (1. - n0));

        Hpp = (site[i].Delta / sqrt(n0 * (1. - n0))) * Mpp + U * Upp + (2. * Lmd - site[i].lambda) * Npp;
                
        arma::vec E(3);
        arma::cx_mat v(3, 3);
        
        arma::eig_sym(E, v, Hpp);
        
        int idx_min = 0;
        
        double e_min = E(0);
        for(int k=1; k<3; k++) {
            if(E(k) < e_min) {
                e_min = E(k);
                idx_min = k;
            }
        }
        
        f(0, i) = v(0, idx_min);
        f(1, i) = v(1, idx_min) / sqrt(2.);
        f(2, i) = v(2, idx_min);
        
    }
    
    return f;
}

double Square_lattice::adjust_lambda(arma::cx_mat & Dm, arma::cx_mat & f) {
    double sum = 0;
    for(int i=0; i<Ns; i++) {
        double n0 = real(Dm(i, i));
        double npp = real(conj(f(1, i)) * f(1, i) + conj(f(2, i)) * f(2, i));
        
//        int sgn = (n0 > npp) ? +1 : -1;
//        site[i].lambda += sgn * 0.1 * pow(fabs(n0 - npp), 0.3); // 0.01 * sgn * pow(fabs(n0 - npp), 0.8);
        
        site[i].lambda += 1. * (n0 - npp);
        
        sum += fabs(n0 - npp);
    }
    
    return sum/((double) Ns);
}

void Square_lattice::test(void) {
    
    int ix = 2;
    
    arma::cx_mat Hpp(3, 3);
    double n0 = real(Density_Mat(ix, ix));
    
    Hpp = (2. * site[ix].Delta / sqrt(n0 * (1. - n0))) * Mpp + U * Upp - site[ix].lambda * Npp;
    
    arma::vec E(3);
    arma::cx_mat v(3, 3);
    
    arma::eig_sym(E, v, Hpp);
    
    cout << Hpp << endl;

    cout << E << endl;
    cout << v << endl;
}

int Square_lattice::compute_GA_solution(void) {

    arma::cx_mat Dm(dim, dim);
    arma::cx_mat f(3, dim);         // Gutzwiller variational parameter
    arma::cx_mat f_prev(3, dim);    // Gutzwiller variational parameter

    kT = GA_kT_start;    

    error = 100;
    err_f = 100;
    iter = 0;

    double error_conv = 1.e-6;

    reset_GA_parameters();
    
    //for(int i=0; i<Ns; i++) site[i].lambda = -0.999 * onsite_V[i];
    
    Hamiltonian = build_Hamiltonian();
    Dm = compute_tV_density_matrix(Hamiltonian);
    f = init_uncorrelated_Phi(Dm);
    f_prev = f;
    
    double tag_zeroR = 10;
    
    while (((iter < max_iter && error > error_conv) || kT > GA_kT_stop) && tag_zeroR > 1.e-6) {

        compute_renormalizations(Dm, f);
        Hamiltonian = build_Hamiltonian();
        
        Dm = compute_tV_density_matrix(Hamiltonian);
        compute_Deltas(Dm);
        f = solve_Phi(Dm);
        
        auto df = f - f_prev;
        err_f = norm(df)/(3.*dim);
        f_prev = f;
        
        tag_zeroR = 0;
        for (int i=0; i<Ns; i++) 
            tag_zeroR += real(f(2, i) * conj(f(2, i)));
        tag_zeroR /= ((double) Ns);
        
        error = adjust_lambda(Dm, f);
        
        iter++;
        
        if (iter % GA_n_annealing == 0 && kT > GA_kT_stop)
            kT *= GA_r_annealing;
 
        
        /* if (iter % 20 == 0) {
            ofstream fs("t.dat");
            fs.precision(10);
            for(int i=0; i<Ns; i++) {
                fs << i << '\t' << real(site[i].R) << '\t' << real(f(0, i)) << '\t' << real(f(1, i)) << '\t' << real(f(2, i)) << '\t' << real(Dm(i, i)) 
                << '\t' << real(conj(f(1, i)) * f(1, i) + conj(f(2, i)) * f(2, i)) << '\t' << site[i].lambda << '\t' << onsite_V[i] << endl;
            }
            fs.close();
        } */
        
//        ofstream fxx("m.dat", ios::app);
//        int ix = 5;
//        arma::cx_mat Hpp(3, 3);
//        double n0 = real(Dm(ix, ix));
//        cout << "n0 = " << n0 << endl;
//        compute_Deltas(Dm);
//        double e_tnn = 0;
//        cx_double _delta = 0;
//        for(int k=0; k<N_nn1; k++) {
//            int j = site[ix].nn1[k]->idx;
//            e_tnn += t1 * real(Dm(j, ix));
//            _delta += t1 * conj(site[j].R) * Dm(j, ix);
//        }
//        cout << "e_k = " << e_tnn << "\t R = " << real(site[ix].R) << "\t Dlt = " << _delta.real() << "," << _delta.imag() << endl;
//        
//        Hpp = (site[ix].Delta / sqrt(n0 * (1. - n0))) * Mpp + U * Upp - site[ix].lambda * Npp;
//        cout << Hpp << endl;
//        fxx.close();
        
    }
    
    if ((iter >= max_iter) && ((error > error_conv) || (tag_zeroR > 1.e-6))) {
        //cout << "unable to reach convergence in this case" << endl;

        Density_Mat = compute_density_matrix(Hamiltonian);
        

        return 0;
    } 
    else {
        //cout << "reach convergence, ";
        //cout << "computing full density matrix ..." << endl;
        
        Density_Mat = compute_density_matrix(Hamiltonian);
        Phi = f;
        
        return 1;
    }
    
}

void Square_lattice::analyze_data(void) {
    double sum1_R2 = 0, sum2_R2 = 0;
    double sum1_nd = 0, sum2_nd = 0;
    double sum1_np = 0, sum2_np = 0;
    double sum1_d = 0, sum2_d = 0;
    
    int nb = 0;
    for(int i=0; i<Ns; i++) {
        
        double nr = real(Density_Mat(i, i));
        sum1_nd += nr;
        sum2_nd += pow(nr, 2);
        
        double npp = real(conj(Phi(1, i)) * Phi(1, i) + conj(Phi(2, i)) * Phi(2, i));
        sum1_np += npp;
        sum2_np += pow(npp, 2);
        
        double db = real(conj(Phi(2, i)) * Phi(2, i));
        sum1_d += db;
        sum2_d += pow(db, 2);
        
        if((site[i].x + site[i].y) % 2 == 0) {
            
            for(int k=0; k<N_nn1; k++) {
                int j = site[i].nn1[k]->idx;
                
                cx_double tmp = site[i].R * conj(site[j].R);
                
                double tmp2 = real(tmp * conj(tmp));
                
                sum1_R2 += tmp2;
                sum2_R2 += pow(tmp2, 2);
                nb++;
            }
        }
    }
    
    sum1_R2 /= ((double) nb);
    sum2_R2 /= ((double) nb);
    
    avg_R2 = sum1_R2;
    sgm_R2 = sqrt(sum2_R2 - sum1_R2 * sum1_R2);
    
    sum1_nd /= ((double) Ns);
    sum2_nd /= ((double) Ns);
    sum1_np /= ((double) Ns);
    sum2_np /= ((double) Ns);
    sum1_d /= ((double) Ns);
    sum2_d /= ((double) Ns);
    
    avg_nd = sum1_nd;
    sgm_nd = sqrt(sum2_nd - sum1_nd * sum1_nd);
    
    avg_np = sum1_np;
    sgm_np = sqrt(sum2_np - sum1_np * sum1_np);
    
    avg_d = sum1_d;
    sgm_d = sqrt(sum2_d - sum1_d * sum1_d);
}

void Square_lattice::compute_electronic_energy(void) {
    
    cx_double sum = 0;
    cx_double sum2 = 0;
    
    for(int i=0; i<Ns; i++) {
        cx_double tmp;
        
        tmp = Hamiltonian(i, i);
        sum += (tmp - site[i].lambda) * Density_Mat(i, i);
//        sum += (tmp) * Density_Mat(i, i);
        
        for (int k=0; k<N_nn1; k++) {
            int j = site[i].nn1[k]->idx;
            tmp = Hamiltonian(i, j);
            sum += tmp * Density_Mat(j, i);
        }

        for (int k=0; k<N_nn2; k++) {
            int j = site[i].nn2[k]->idx;
            tmp = Hamiltonian(i, j);
            sum += tmp * Density_Mat(j, i);
        }
        
        sum2 += Phi(2, i) * conj(Phi(2, i));
    }
    
    
    E_e = (2. * sum.real() + U * sum2.real()) / ((double) Ns);
}

void Square_lattice::statistics_GA(int ndata, double _u, double W) {
    
    U = _u;
    
    ofstream fs("r.dat");

    double accu_avg_R2 = 0, accu_sgm_R2 = 0;
    double accu_avg_d = 0, accu_sgm_d = 0;
    double accu_avg_nd = 0, accu_sgm_nd = 0;
    double accu_avg_np = 0, accu_sgm_np = 0;
    double accu_Ee = 0;
    
    fs.precision(10);
    
    int r = 0;
    for(int q=0; q<ndata; q++) {
        
        init_quenched_disorder(W);

        int converge = compute_GA_solution();
        
        if (converge == 1) {
        
            double n_tot = real( arma::trace(Density_Mat) ) / ((double) dim);
            analyze_data();
            compute_electronic_energy();
            
            fs << r << '\t' << avg_R2 << '\t' << sgm_R2 << '\t' << avg_d << '\t' << sgm_d << '\t';
            fs << avg_nd << '\t' << sgm_nd << '\t' << avg_np << '\t' << sgm_np << '\t';
            fs << E_e << '\t';
            fs << n_tot << endl;
            
            int n = r;
            accu_avg_R2 = (n * accu_avg_R2 + avg_R2) / (n + 1.);
            accu_sgm_R2 = (n * accu_sgm_R2 + sgm_R2) / (n + 1.);
            accu_avg_d  = (n * accu_avg_d  + avg_d ) / (n + 1.);
            accu_sgm_d  = (n * accu_sgm_d  + sgm_d ) / (n + 1.);
            accu_avg_nd = (n * accu_avg_nd + avg_nd) / (n + 1.);
            accu_sgm_nd = (n * accu_sgm_nd + sgm_nd) / (n + 1.);
            accu_avg_np = (n * accu_avg_np + avg_np) / (n + 1.);
            accu_sgm_np = (n * accu_sgm_np + sgm_np) / (n + 1.);
            accu_Ee     = (n * accu_Ee     + E_e)    / (n + 1.);
            
            ofstream fx("a.dat");
            fx.precision(10);
            fx << accu_avg_R2 << '\t' << accu_sgm_R2 << '\t' << accu_avg_d << '\t' << accu_sgm_d << '\t';
            fx << accu_avg_nd << '\t' << accu_sgm_nd << '\t' << accu_avg_np << '\t' << accu_sgm_np << '\t';
            fx << accu_Ee << endl;
            fx.close();
            
            fx.open("ga" + std::to_string(r) + ".dat");
            fx.precision(10);
            for(int i=0; i<Ns; i++) {
                fx << site[i].x << '\t' << site[i].y << '\t';
                fx << onsite_V[i] << '\t' << site[i].lambda << '\t';
                fx << real(Density_Mat(i, i)) << '\t';
                fx << real(conj(Phi(2, i)) * Phi(2, i)) << '\t';
                fx << real(site[i].R) << '\t'; // << imag(site[i].R) << '\t';
                fx << endl;
            }
            fx.close();
         
            r++;
        }
    }
    
    fs.close();
}

void Square_lattice::print_config(string const filename) {
    std::ofstream fx;
        
    fx.open(filename.c_str(), ios::out);
    fx.precision(8);

    for (int i=0; i<Ns; i++) {
        fx << site[i].x << '\t' << site[i].y << '\t';
        fx << onsite_V[i] << '\t' << site[i].lambda << '\t';
        fx << force[i] << '\t' << 2 * real(Density_Mat(i, i)) << '\t';
        fx << real(conj(Phi(2, i)) * Phi(2, i)) << '\t';
        fx << real(site[i].R) << '\t';
        fx << iter << '\t' << error << '\t' << err_f << '\t';
        fx << endl;
    }
    fx.close();
}


void Square_lattice::step_NVT() {
    move_atoms();

    if (solver_type == 0)
        compute_GA_solution();

    if (solver_type == 1) {
        for (int i=0; i<Ns; i++) {
            site[i].R = 1.0;
             site[i].lambda = 0;
        }
        Hamiltonian = build_Hamiltonian();
        Density_Mat = compute_density_matrix(Hamiltonian);
    }

    compute_forces();
    integrate_forces();
    integrate_Langevin();

}


void Square_lattice::move_atoms() {
    for (int i=0; i<Ns; i++) {
        double dlt = 0.5 * (force[i]/mass) * dt;
        dlt += velocity[i];
        velocity[i] = dlt;  // velocity at t + 0.5*dt
        dlt *= dt;
        onsite_V[i] += dlt;
    }
}

void Square_lattice::integrate_forces() {
    for (int i=0; i<Ns; i++)
        velocity[i] += 0.5 * (force[i]/mass) * dt;
}

void Square_lattice::integrate_Langevin() {
    std::normal_distribution<double> rd;    // default mean = 0, var = 1
    
    double alpha2 = exp(-gamma * dt);
    double sigma2 = sqrt((1 - pow(alpha2, 2)) * kT / mass);
    
    for (int i=0; i<Ns; i++)
        velocity[i] = alpha2 * velocity[i] + sigma2 * rd(rng);
}

void Square_lattice::compute_forces(void) {
    force = vector<double>(Ns, 0);

    for (int i=0; i<Ns; i++) {
        force[i] = g * (2 * real(Density_Mat(i, i)) - 1) - k0 * onsite_V[i];

        for (int k=0; k<N_nn1; k++) {
            int j = site[i].nn1[k]->idx;
            force[i] += -kappa * onsite_V[j];
        }
    }
}




