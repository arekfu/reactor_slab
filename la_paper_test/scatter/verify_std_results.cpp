/*
 *  Test of "Direct Sampling of Monte Carlo Flight Paths In Media
 *  With Continuously Varying Cross-Sections" (LA-UR-02-6530)
 *  by Forrest Brown and Bill Martin.
 *
 *  Implements solution of direct path lenght, delta tracking,
 *  and meshed delta tracking.
 *
 * */

#include"random_pcg.hpp"

#include<cmath>
#include<iomanip>
#include<iostream>
#include<memory>
#include<vector>
#include<fstream>
#include<omp.h>

using namespace lmct; // namespace for my personal PCG RNG wrapper

const double EPS = 1e-6;
const int NPART = 1e7;
const int NBIN = 5;
const int NFOMBINS = 100;
const double Fdx = 2.0/(double)NFOMBINS;
const double dx = 2.0/(double)NBIN;
const double p = 0.75; // Esamp / Emaj for full system
const double p_mshd = 0.75; // Esamp[i] / Emaj[i] for bin
const double q = 0.3;
const double q_mshd = 0.3;

const double P_abs = 0.3;
const double P_sct = 0.7;
const double P_straight_ahead = 0.5;

const double wgt_cutoff = 0.25;
const double wgt_survival = 1.0;
const double wgt_split = 2.0;
const double alpha = 1.0;

// Constants for gaussians A*exp(-a*(x-z)^2)
const double A = 2.0/std::sqrt(2.0*M_PI);
const double a_s = (1.0/0.05)*(1.0/0.05);
const double a_b = 1.0;
const double z = 1.23;

// Seed variables for PCG RNG
std::vector<uint64_t> pcg_seeds;

// Vectors to hold avg Coll rates
std::vector<double> avg_coll_rates_real;
std::vector<double> avg_coll_rates_all;
const int NTRIALS = 20;

// Tallies
double collide; // Sum of weights that collide
double collide_sqr; // Sum of weights squared that collide
double all_collide; // Sum of all weights that collide at all collisions (virt. and real)
double all_collide_sqr; // Summ of all weights squared that collide (virt. and real)
double escape; // Sum of weights that leak
double escape_sqr; // Sum of weights squared that lead
double xs_evals; // # of xs look ups or xs integrations
double wgt_chngs; // # of times particle wgt sign flipped
std::vector<std::vector<double>> coll_density; //[0] #sum coll in box,[1] sum coll sqr, [2] coll STD, [3] Coll FOM in box
std::vector<std::vector<double>> all_coll_density; // Same as above but scores at all collisions (virt. and real)

// Constants for Gaussian XS functions
const double T_2_s = A*std::sqrt(M_PI)*((std::erf(std::sqrt(a_s)*(2.0 - z)) 
                     - std::erf(std::sqrt(a_s)*(-z)))/(2.0*std::sqrt(a_s)));

const double T_2_b = A*std::sqrt(M_PI)*((std::erf(std::sqrt(a_b)*(2.0 - z)) 
                     - std::erf(std::sqrt(a_b)*(-z)))/(2.0*std::sqrt(a_b)));


// Base Cross Section Class
class XS {
    public:
        XS(double _Emax) {
            Emax = _Emax;
        }

        // All virtual methods
        virtual double Et(double /*x*/) {return -1.0;}

        double P_nc() {return Pnc;}

        // Data
        double Pnc; // Prob. of no collision
        double G; // 1 - Pnc (Prob. of collision)
        double Emax; // Maj over region
        double Em[NBIN]; // Maj in each bin
        double Esmp[NBIN]; // Sampling XS in each bin for NWDT and MNWDT
};

// Constant XS (Et = 1.0 ∀ 0<x<2, or Et=0)
class Constant : public XS {
    public:
        Constant():XS(1.0) {
            for(int b = 0; b < NBIN; b++) {
                Em[b] = 1.0;
                Esmp[b] = p_mshd;
            }
        }

        double Et(double /*x*/) {return 1.0;}

}; // Constants

// Step XS (Et = 1.5 ∀ 0<x<0.278), Et=1.0 ∀ 0.278<x<2)
class Step : public XS {
    public:
        Step():XS(1.5) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double Et(double x) {
            if((x >= 0.0) and (x < 0.278)) {return 1.5;}
            else if ((x >= 0.278) and (x < 2.0)) {return 1.0;}
            else {exit(1);}
        }

}; // Step

class Lin_Decrease : public XS {
    public:
        Lin_Decrease():XS(2.0) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }    
        }

        double Et(double x) {
            return 2.0 - x;
        }
};

class Lin_Increase : public XS {
    public:
        Lin_Increase():XS(2.0) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx + dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double Et(double x) {
            return x;
        }
};

class Exp_Decrease : public XS {
    public:
        Exp_Decrease():XS(1.0) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double Et(double x) {
            return std::exp(-3.0*x);
        }
};

class Exp_Increase : public XS {
    public:
        Exp_Increase():XS(0.1*std::exp(4.0)) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx + dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double Et(double x) {
            return 0.1*std::exp(2.0*x);
        }
};

class Gauss_Sharp : public XS {
    public:
        Gauss_Sharp():XS(A) {
            int bin_peak = std::floor(z/dx);
            double x;
            for(int b = 0; b < NBIN; b++) {
                if(b < bin_peak) {
                    x = b*dx + dx;
                }
                else if(b > bin_peak) {
                    x = b*dx;
                }
                else {
                    x = z;
                }
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double Et(double x) {
            return A*std::exp(-a_s*(x-z)*(x-z));
        }
};

class Gauss_Broad : public XS {
    public:
        Gauss_Broad():XS(A) {
            int bin_peak = std::floor(z/dx);
            double x;
            for(int b = 0; b < NBIN; b++) {
                if(b < bin_peak) {
                    x = b*dx + dx;
                }
                else if(b > bin_peak) {
                    x = b*dx;
                }
                else {
                    x = z;
                }
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double Et(double x) {
            return A*std::exp(-a_b*(x-z)*(x-z));
        }
};

void roulette(double& wgt, bool& alive, const double& xi) {
    if(std::abs(wgt) < wgt_cutoff) {
        double P_kill = 1.0 - (std::abs(wgt)/wgt_survival);
        if(xi < P_kill) alive = false;
        else {
            if(wgt > 0) wgt = wgt_survival;
            else wgt = -wgt_survival;
        }
    }
}

void score_real_collision(double& wgt, double& x) {
    #pragma omp atomic
    collide += wgt;
    #pragma omp atomic
    collide_sqr += wgt*wgt;

    int coll_bin = std::floor(x/Fdx);
    #pragma omp atomic
    coll_density[coll_bin][0] += wgt;
    #pragma omp atomic
    coll_density[coll_bin][1] += wgt*wgt;
}

void score_all_collision(double& wgt, double& x) {
    #pragma omp atomic
    all_collide += wgt;
    #pragma omp atomic
    all_collide_sqr += wgt*wgt;

    int coll_bin = std::floor(x/Fdx);
    #pragma omp atomic
    all_coll_density[coll_bin][0] += wgt;
    #pragma omp atomic
    all_coll_density[coll_bin][1] += wgt*wgt;
}

void score_escape(double& wgt) {
    #pragma omp atomic
    escape += wgt;
    #pragma omp atomic
    escape_sqr += wgt*wgt;
}

void Delta_Tracking(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Delta Tracking\n";

    int cnts_sum = 0;

    #pragma omp parallel
    {    
        PCG rng;
        int thread_id = omp_get_thread_num();
        uint64_t pcg_seed;
        #pragma omp atomic read
        pcg_seed = pcg_seeds[thread_id];
        rng.seed(pcg_seed);

        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            bool virtual_collision;
            double x = 0.0;
            double u = 1.0;
            int cnt = 0;
            double wgt = 1.0;
            bool alive = true;
            while(alive) {
                virtual_collision = true;
                while(virtual_collision) {
                    double d = -std::log(rng.rand())/(xs->Emax);
                    x += u*d;
                    if((x >= 2.0) or (x <= 0.0)) {
                        score_escape(wgt); 
                        virtual_collision = false;
                        alive = false;
                        #pragma omp atomic
                        cnts_sum += cnt;
                    } else {
                        cnt++;
                        if(rng.rand() < (xs->Et(x)/xs->Emax)) {
                            virtual_collision = false;
                        }
                        // Score every collision
                        double score = wgt*xs->Et(x) / xs->Emax;
                        score_all_collision(score,x);
                    }
                }
                
                if(alive) {
                    // Score real collision
                    score_real_collision(wgt,x); 

                    // Implicit capture
                    wgt *= 1.0 - P_abs; // Equivalent to 1.0 - (Ea/Et)

                    // Scatter
                    if(rng.rand() > P_straight_ahead) u *= -1;

                    // Russian Roulette
                    double xi = rng.rand();
                    roulette(wgt, alive, xi);
                    if(alive == false) {
                        #pragma omp atomic
                        cnts_sum += cnt;
                    }

                } // If alive for real collision
            } // While alive
        } // For all particles
        #pragma omp atomic write
        pcg_seeds[thread_id] = rng.get_seed();
    } // Parallel
    xs_evals += cnts_sum;
}

void Meshed_Delta_Tracking(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Meshed Delta Tracking\n";

    int cnts_sum = 0;
    int virtual_cnt_sum = 0;
    int bin_cnt_sum = 0;

    #pragma omp parallel
    {    
        PCG rng;
        int thread_id = omp_get_thread_num();
        uint64_t pcg_seed;
        #pragma omp atomic read
        pcg_seed = pcg_seeds[thread_id];
        rng.seed(pcg_seed);

        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            bool virtual_collision = true;
            double x = 0.0;
            double u = 1.0;
            double wgt = 1.0;
            double Emax = xs->Em[0];
            int bin = 0;
            double d,d_bin;
            int cnt = 0;
            int virtual_cnt = 0;
            int bin_cnt = 0;
            bool alive = true;
            while(alive) {
                virtual_collision = true;
                while(virtual_collision) {
                    if(u == -1.0) d_bin = x - static_cast<double>(bin)*dx;
                    else d_bin = (static_cast<double>(bin)*dx + dx) - x;
                    d = -std::log(rng.rand())/Emax;
                    if(d_bin < d) {
                        bin_cnt++;
                        d = d_bin + 1e-6;
                        x += u*d;
                        bin = std::floor(x/dx);
                        if((x >= 2.0) or (x  <= 0.0)) {
                            score_escape(wgt);
                            virtual_collision = false;
                            alive = false;
                            
                            // Score other tallies 
                            #pragma omp atomic
                            cnts_sum += cnt;
                            #pragma omp atomic
                            virtual_cnt_sum += virtual_cnt;
                            #pragma omp atomic
                            bin_cnt_sum += bin_cnt;
                        }
                        else {
                            Emax = xs->Em[bin];
                        }
                    }
                    else {
                        cnt++;
                        x += u*d;
                        double xi = rng.rand();
                        double Pr = xs->Et(x)/Emax;
                        if(xi < Pr) {
                            virtual_collision = false;
                        } else {virtual_cnt++;}
                        
                        double score = wgt*xs->Et(x) /Emax;
                        score_all_collision(score, x);
                    }
                }// While virtual

                if(alive) {
                    // Score real collision
                    score_real_collision(wgt,x);
                   
                    // Implicit capture
                    wgt *= 1.0 - P_abs;
                    
                    // Scatter
                    if(rng.rand() > P_straight_ahead) u *= -1;

                    // Russian Roulette
                    double xi = rng.rand();
                    roulette(wgt, alive, xi);
                    if(alive == false){
                        #pragma omp atomic
                        cnts_sum += cnt;
                        #pragma omp atomic
                        virtual_cnt_sum += virtual_cnt;
                        #pragma omp atomic
                        bin_cnt_sum += bin_cnt;
                    }
                }
            } // While alive
        }// For all particles
        #pragma omp atomic write
        pcg_seeds[thread_id] = rng.get_seed();
    }// Parallel
    xs_evals += cnts_sum;
}

void Negative_Weight_Delta_Tracking(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    double sign_change = 0.0;

    // Particle bank vectors
    std::vector<double> Positions;
    std::vector<double> Directions;
    std::vector<double> Wgts;
    int n_particles = NPART;
    // Generate initial source
    for(int i = 0; i < n_particles; i++) {
        Positions.push_back(0.0);
        Directions.push_back(1.0);
        Wgts.push_back(1.0);
    }
    std::vector<double> Split_Positions;
    std::vector<double> Split_Direction;
    std::vector<double> Split_Wgts;

    while(n_particles > 0) {
        //std::cout << " nparticles = " << n_particles << "\n";
        #pragma omp parallel
        {    
            PCG rng;
            int thread_id = omp_get_thread_num();
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                double Esamp = p*(xs->Emax);
                double x = Positions[n];
                double u = Directions[n];
                double wgt = Wgts[n];
                bool alive = true;
                int cnt = 0;
                while(alive) {
                    double d = -std::log(rng.rand())/Esamp;
                    x += u*d;
                    
                    if((x >= 2.0) or (x <= 0.0)) {
                        score_escape(wgt);

                        alive = false;
                        #pragma omp atomic
                        cnts_sum += cnt;
                    }
                    else if(rng.rand() < q) {
                        // Collision is real, addjust weight
                        wgt *= xs->Et(x)/(Esamp * q);
                        double score = wgt *q;
                        score_all_collision(score,x);
                        score_real_collision(wgt,x);
                        cnt++;
                        
                        // Implicit capture
                        wgt *= 1.0 - P_abs;
                        
                        // Scatter
                        if(rng.rand() > P_straight_ahead) {u *= -1;}

                        // Russian Roulette
                        double xi = rng.rand();
                        roulette(wgt,alive,xi);
                        if(alive == false) {
                            #pragma omp atomic
                            cnts_sum += cnt;
                        }
                        
                    }
                    else { // Collision is virtual
                        double dw = ((1.0 - (xs->Et(x)/Esamp))/(1.0 - q));
                        if(dw < 0.0) {
                            #pragma omp atomic
                            sign_change += 1.0;
                        }
                        double score = wgt*(xs->Et(x)/Esamp);
                        score_all_collision(score,x);
                        wgt = wgt*dw;
                        cnt++;
                    }

                    // Split if needed
                    if(alive and (std::abs(wgt) >= wgt_split)) {
                        double n_new = std::round(std::abs(wgt));
                        wgt /= n_new;
                        for(int j = 0; j < static_cast<int>(n_new-1); j++) {
                            #pragma omp critical
                            {
                                Split_Positions.push_back(x);
                                Split_Direction.push_back(u);
                                Split_Wgts.push_back(wgt);
                            }
                        }
                    }// split
                } // While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();
        }// Parallel
        n_particles = static_cast<int>(Split_Positions.size());
        Positions = Split_Positions;
        Directions = Split_Direction;
        Wgts = Split_Wgts;
        
        Split_Positions.clear();
        Split_Direction.clear();
        Split_Wgts.clear();
    } // While still split particles
    
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Meshed_Negative_Weight_Delta_Tracking(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Meshed Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    int bin_cnt_sum = 0;
    double sign_change = 0.0;

    // Particle bank vectors
    std::vector<double> Positions;
    std::vector<double> Directions;
    std::vector<double> Wgts;
    std::vector<int> Bins;
    int n_particles = NPART;
    // Generate initial source
    for(int i = 0; i < n_particles; i++) {
        Positions.push_back(0.0);
        Directions.push_back(1.0);
        Wgts.push_back(1.0);
        Bins.push_back(0);
    }
    std::vector<double> Split_Positions;
    std::vector<double> Split_Direction;
    std::vector<double> Split_Wgts;
    std::vector<int> Split_Bins;

    while(n_particles > 0) {
        //std::cout << " nparticles = " << n_particles << "\n";
        #pragma omp parallel
        {    
            PCG rng;
            int thread_id = omp_get_thread_num();
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                bool alive = true;
                double x = Positions[n];
                double u = Directions[n];
                int bin = Bins[n];
                double Esamp = xs->Esmp[bin];
                double d,d_bin;
                int cnt = 0;
                int bin_cnt = 0;
                double wgt = Wgts[n];
                while(alive) {
                    if(u == -1.0) d_bin = x - static_cast<double>(bin)*dx;
                    else d_bin = (static_cast<double>(bin)*dx + dx) - x;
                    d = -std::log(rng.rand())/Esamp;
                    if(d_bin < d) {
                        bin_cnt++;
                        d = d_bin + 1e-6;
                        x += u*d;
                        bin = std::floor(x/dx);
                        if((x >= 2.0) or (x <= 0.0)) {
                            score_escape(wgt);
                            alive = false;
                            
                            #pragma omp atomic
                            cnts_sum += cnt;
                            #pragma omp atomic
                            bin_cnt_sum += bin_cnt;
                        }
                        else {
                            Esamp = xs->Esmp[bin];
                        }
                    }
                    else {
                        x += u*d;
                        if(rng.rand() < q_mshd) { // Real collision
                            // update weight
                            wgt *= (xs->Et(x)/(Esamp*q_mshd));
                            double score = wgt*(q_mshd);
                            score_all_collision(score,x);
                            score_real_collision(wgt,x);
                            cnt++;

                            // Implicit capture
                            wgt *= 1.0 - P_abs;

                            // Scatter
                            if(rng.rand() > P_straight_ahead) u *= -1;

                            // Russian Roulette
                            double xi = rng.rand();
                            roulette(wgt,alive,xi);
                            if(alive == false) {
                                #pragma omp atomic
                                cnts_sum += cnt;
                                #pragma omp atomic
                                bin_cnt_sum += bin_cnt;
                            }
                            
                        }
                        else {
                            cnt++;
                            double dw=((1.0-(xs->Et(x)/Esamp))/(1.0-q_mshd));
                            if(dw < 0.0) {
                                #pragma omp atomic
                                sign_change += 1.0;
                            }
                            double score = wgt*(xs->Et(x)/Esamp);
                            score_all_collision(score,x);
                            wgt = wgt*dw;
                        }

                        // Split if needed
                        if(alive and (std::abs(wgt) >= wgt_split)) {
                            double n_new = std::round(std::abs(wgt));
                            wgt /= n_new;
                            for(int j = 0; j < static_cast<int>(n_new-1); j++) {
                                #pragma omp critical
                                {
                                    Split_Positions.push_back(x);
                                    Split_Direction.push_back(u);
                                    Split_Wgts.push_back(wgt);
                                    Split_Bins.push_back(bin);
                                }
                            }
                        }// split
                    }
                }// While alive
            } // For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();
        }// Parallel
        
        n_particles = static_cast<int>(Split_Positions.size());
        Positions = Split_Positions;
        Directions = Split_Direction;
        Wgts = Split_Wgts;
        Bins = Split_Bins;
        
        Split_Positions.clear();
        Split_Direction.clear();
        Split_Wgts.clear();
        Split_Bins.clear();

    }// While still split particles

    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Bomb_Transport(std::unique_ptr<XS> const &xs, double P) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Bomb Paper Transport, p = " << P << "\n";
    int cnts_sum = 0;
    double sign_change = 0.0;
    
    // Particle bank vectors
    std::vector<double> Positions;
    std::vector<double> Directions;
    std::vector<double> Wgts;
    int n_particles = NPART;
    // Generate initial source
    for(int i = 0; i < n_particles; i++) {
        Positions.push_back(0.0);
        Directions.push_back(1.0);
        Wgts.push_back(1.0);
    }
    std::vector<double> Split_Positions;
    std::vector<double> Split_Direction;
    std::vector<double> Split_Wgts;

    while(n_particles > 0) {
        //std::cout << " nparticles = " << n_particles << "\n";
        #pragma omp parallel
        {
            PCG rng;
            int thread_id = omp_get_thread_num();
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                double Esmp = P*xs->Emax;
                double x = Positions[n];
                double u = Directions[n];
                bool alive = true;
                int cnt = 0;
                double wgt = Wgts[n];
                bool real_collision = false;

                while(alive) {
                    double d = -std::log(rng.rand())/Esmp;
                    x += u*d;
                    real_collision = false;
                    
                    // Fist check for leak
                    if((x >= 2.0) or (x <= 0.0)) {
                        score_escape(wgt);
                        alive = false;
                        #pragma omp atomic
                        cnts_sum += cnt;
                    } else {
                        double E_tot = xs->Et(x);
                        cnt += 1;
                        if(E_tot > Esmp) { // First negative branch
                            //double D_alpha = alpha*E_tot / ((2.0 + alpha)*E_tot - Esmp);
                            double D = E_tot / (2*E_tot - Esmp);
                            double F = (E_tot / (D*Esmp));

                            double score = wgt*E_tot/Esmp;
                            score_all_collision(score,x);

                            //if(rand(rng) < D_alpha) {
                            if(rng.rand() < D) {
                                real_collision = true;
                                wgt *= F;
                                //wgt *=  F*(D/D_alpha);
                            }
                            else {
                                wgt *= -F;
                                //wgt *= -F*((1. - D)/(1. - D_alpha));
                                #pragma omp atomic
                                sign_change += 1.0;
                            }
                            

                        } else { // Delta tracking branch
                            double P_real = E_tot/ Esmp;
                            if(rng.rand() < P_real) {real_collision = true;}
                            double score = wgt*E_tot/Esmp;
                            score_all_collision(score,x);
                        }

                        if(real_collision) {
                            // Score real collision
                            score_real_collision(wgt,x);

                            // Implicit caputure
                            wgt *= 1.0 - P_abs;

                            // Scatter
                            if(rng.rand() > P_straight_ahead) u *= -1;

                            // Russian Roulette
                            double xi = rng.rand();
                            roulette(wgt,alive,xi);
                            if(alive == false) {
                                #pragma omp atomic
                                cnts_sum += cnt;
                            }
                            
                        }// End real coll.
                    }

                    // Split if needed
                    if(alive and (std::abs(wgt) >= wgt_split)) {
                        double n_new = std::floor(std::abs(wgt));
                        wgt /= n_new;
                        for(int j = 0; j < static_cast<int>(n_new-1); j++) {
                            #pragma omp critical
                            {
                                Split_Positions.push_back(x);
                                Split_Direction.push_back(u);
                                Split_Wgts.push_back(wgt);
                            }
                        }
                    }// split

                }// While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();
        }// Parallel

        n_particles = static_cast<int>(Split_Positions.size());
        Positions = Split_Positions;
        Directions = Split_Direction;
        Wgts = Split_Wgts;
        
        Split_Positions.clear();
        Split_Direction.clear();
        Split_Wgts.clear();
    }// While still particles

    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Meshed_Bomb_Transport(std::unique_ptr<XS> const &xs, double P) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Meshed Bomb Paper Transport, p = " << P << "\n";
    int cnts_sum = 0;
    double sign_change = 0.0;

    // Particle bank vectors
    std::vector<double> Positions;
    std::vector<double> Directions;
    std::vector<double> Wgts;
    std::vector<int> Bins;
    int n_particles = NPART;
    // Generate initial source
    for(int i = 0; i < n_particles; i++) {
        Positions.push_back(0.0);
        Directions.push_back(1.0);
        Wgts.push_back(1.0);
        Bins.push_back(0);
    }
    std::vector<double> Split_Positions;
    std::vector<double> Split_Direction;
    std::vector<double> Split_Wgts;
    std::vector<int> Split_Bins;

    while(n_particles > 0) {
        //std::cout << " nparticles = " << n_particles << "\n";
        #pragma omp parallel
        {
            PCG rng;
            int thread_id = omp_get_thread_num();
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                double x = Positions[n];
                double u = Directions[n];
                double d_bin;
                int bin = Bins[n];
                double Esmp = P*xs->Em[bin];
                bool alive = true;
                int cnt = 0;
                double wgt = Wgts[n];
                bool real_collision = false;

                while(alive) {
                    double d = -std::log(rng.rand())/Esmp;
                    if(u == -1.0) d_bin = x - static_cast<double>(bin)*dx;
                    else d_bin = (static_cast<double>(bin)*dx + dx) - x;
                    real_collision = false;
                    
                    if(d_bin < d) {
                        x += u*(d_bin+1e-6);
                        if((x >= 2.0) or (x <= 0.0)) {
                            alive = false;
                            score_escape(wgt);
                            #pragma omp atomic
                            cnts_sum += cnt;
                        } else {
                            bin = std::floor(x/dx);
                            Esmp = P*xs->Em[bin];
                        }

                    } else {
                        x += u*d;
                        // Fist check for leak
                        if((x >= 2.0) or (x <= 0.0)) {
                            alive = false;
                            score_escape(wgt);
                            #pragma omp atomic
                            cnts_sum += cnt;
                        } else {
                            double E_tot = xs->Et(x);
                            cnt += 1;
                            if(E_tot > Esmp) { // First negative branch
                                double D = E_tot / (2*E_tot - Esmp);
                                double F = E_tot / (D*Esmp);
                                double score = wgt*E_tot/Esmp;
                                score_all_collision(score,x);
                                wgt *= F;
                                if(rng.rand() < D) {real_collision = true;}
                                else {
                                    wgt *= -1.0;
                                    #pragma omp atomic
                                    sign_change += 1.0;
                                }

                            } else { // Delta tracking branch
                                double P_real = E_tot/ Esmp;
                                double score = wgt*E_tot/Esmp;
                                score_all_collision(score,x);
                                if(rng.rand() < P_real) {real_collision = true;}
                            }

                            if(real_collision) {
                                // Score real collision
                                score_real_collision(wgt,x);

                                // Implicit capture
                                wgt *= 1.0 - P_abs;

                                // Scatter
                                if(rng.rand() > P_straight_ahead) u *= -1;

                                // Russian Roulette
                                double xi = rng.rand();
                                roulette(wgt,alive,xi);
                                if(alive == false) {
                                    #pragma omp atomic
                                    cnts_sum += cnt;
                                }

                            }

                            // Split if needed
                            if(alive and (std::abs(wgt) >= wgt_split)) {
                                double n_new = std::round(std::abs(wgt));
                                wgt /= n_new;
                                for(int j = 0;j<static_cast<int>(n_new-1);j++) {
                                    #pragma omp critical
                                    {
                                        Split_Positions.push_back(x);
                                        Split_Direction.push_back(u);
                                        Split_Wgts.push_back(wgt);
                                        Split_Bins.push_back(bin);
                                    }
                                }
                            }// split
                        }
                    }
                }// While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = pcg_seed;
        }// Parallel

        n_particles = static_cast<int>(Split_Positions.size());
        Positions = Split_Positions;
        Directions = Split_Direction;
        Wgts = Split_Wgts;
        Bins = Split_Bins;
        
        Split_Positions.clear();
        Split_Direction.clear();
        Split_Wgts.clear();
        Split_Bins.clear();
    } // While split particles

    
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Improving_Meshed_Bomb_Transport(std::unique_ptr<XS> const &xs, double P) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Improving Meshed Bomb Paper Transport, p = " << P << "\n";
    int cnts_sum = 0;
    double sign_change = 0.0;

    // Particle bank vectors
    std::vector<double> Positions;
    std::vector<double> Directions;
    std::vector<double> Wgts;
    std::vector<int> Bins;
    int n_particles = NPART;
    // Generate initial source
    for(int i = 0; i < n_particles; i++) {
        Positions.push_back(0.0);
        Directions.push_back(1.0);
        Wgts.push_back(1.0);
        Bins.push_back(0);
    }
    std::vector<double> Split_Positions;
    std::vector<double> Split_Direction;
    std::vector<double> Split_Wgts;
    std::vector<int> Split_Bins;

    std::vector<double> Em;
    for(int i = 0; i < NBIN; i++) {
        Em.push_back(P*xs->Em[i]);
    }

    while(n_particles > 0) {
        #pragma omp parallel
        {
            PCG rng;
            int thread_id = omp_get_thread_num();
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                double x = Positions[n];
                double u = Directions[n];
                double d_bin;
                int bin = Bins[n];
                double Esmp = Em[bin];
                bool alive = true;
                int cnt = 0;
                double wgt = Wgts[n];
                bool real_collision = false;

                while(alive) {
                    double d = -std::log(rng.rand())/Esmp;
                    if(u == -1.0) d_bin = x - static_cast<double>(bin)*dx;
                    else d_bin = (static_cast<double>(bin)*dx + dx) - x;
                    real_collision = false;
                    
                    if(d_bin < d) {
                        x += u*(d_bin+1e-6);
                        if((x >= 2.0) or (x <= 0.0)) {
                            alive = false;
                            score_escape(wgt);
                            #pragma omp atomic
                            cnts_sum += cnt;
                        } else {
                            bin = std::floor(x/dx);
                            Esmp = xs->Em[bin];
                        }

                    } else {
                        x += u*d;
                        // Fist check for leak
                        if((x >= 2.0) or (x <= 0.0)) {
                            alive = false;
                            score_escape(wgt);
                            #pragma omp atomic
                            cnts_sum += cnt;
                        } else {
                            double E_tot = xs->Et(x);
                            cnt += 1;
                            if(E_tot > Esmp) { // First negative branch
                                double D = E_tot / (2*E_tot - Esmp);
                                double F = E_tot / (D*Esmp);
                                double score = wgt*E_tot/Esmp;
                                score_all_collision(score,x);
                                wgt *= F;
                                if(rng.rand() < D) {real_collision = true;}
                                else {
                                    wgt *= -1.0;
                                    #pragma omp atomic
                                    sign_change += 1.0;
                                }

                                

                            } else { // Delta tracking branch
                                double P_real = E_tot/ Esmp;
                                double score = wgt*E_tot/Esmp;
                                score_all_collision(score,x);
                                if(rng.rand() < P_real) {real_collision = true;}
                            }

                            if(real_collision) {
                                // Score real collision
                                score_real_collision(wgt,x);

                                // Implicit capture
                                wgt *= 1.0 - P_abs;

                                // Scatter
                                if(rng.rand() > P_straight_ahead) u *= -1;

                                // Russian Roulette
                                double xi = rng.rand();
                                roulette(wgt,alive,xi);
                                if(alive == false) {
                                    #pragma omp atomic
                                    cnts_sum += cnt;
                                }

                            }
                            
                            if(E_tot > Esmp) {
                                // Improve Esmp
                                Esmp = E_tot;
                                Em[bin] = Esmp;
                            }
                            

                            // Split if needed
                            if(alive and (std::abs(wgt) >= wgt_split)) {
                                double n_new = std::round(std::abs(wgt));
                                wgt /= n_new;
                                for(int j = 0;j<static_cast<int>(n_new-1);j++) {
                                    #pragma omp critical
                                    {
                                        Split_Positions.push_back(x);
                                        Split_Direction.push_back(u);
                                        Split_Wgts.push_back(wgt);
                                        Split_Bins.push_back(bin);
                                    }
                                }
                            }// split
                        }
                    }
                }// While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();
        }// Parallel

        n_particles = static_cast<int>(Split_Positions.size());
        Positions = Split_Positions;
        Directions = Split_Direction;
        Wgts = Split_Wgts;
        Bins = Split_Bins;
        
        Split_Positions.clear();
        Split_Direction.clear();
        Split_Wgts.clear();
        Split_Bins.clear();
    } // While split particles

    
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Previous_XS_Bomb_Transport(std::unique_ptr<XS> const &xs) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Previous XS Bomb Paper Transport\n";
    int cnts_sum = 0;
    double sign_change = 0.0;
    
    // Particle bank vectors
    std::vector<double> Positions;
    std::vector<double> Directions;
    std::vector<double> Wgts;
    std::vector<double> Esmps;
    int n_particles = NPART;
    // Generate initial source
    for(int i = 0; i < n_particles; i++) {
        Positions.push_back(0.0);
        Directions.push_back(1.0);
        Wgts.push_back(1.0);
        Esmps.push_back(3.0*xs->Et(0.0));
    }
    std::vector<double> Split_Positions;
    std::vector<double> Split_Direction;
    std::vector<double> Split_Wgts;
    std::vector<double> Split_Esmps;

    while(n_particles > 0) {
        #pragma omp parallel
        {
            PCG rng;
            int thread_id = omp_get_thread_num();
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                double x = Positions[n];
                double Esmp = Esmps[n];
                if(Esmp < 1.0) Esmp = 1.0;
                double u = Directions[n];
                bool alive = true;
                int cnt = 0;
                double wgt = Wgts[n];
                bool real_collision = false;

                while(alive) {
                    double d = -std::log(rng.rand())/Esmp;
                    x += u*d;
                    real_collision = false;
                    
                    // Fist check for leak
                    if((x >= 2.0) or (x <= 0.0)) {
                        score_escape(wgt);
                        alive = false;
                        #pragma omp atomic
                        cnts_sum += cnt;
                    } else {
                        double E_tot = xs->Et(x);
                        cnt += 1;
                        if(E_tot > Esmp) { // First negative branch
                            double D = E_tot / (2*E_tot - Esmp);
                            double F = E_tot / (D*Esmp);
                            wgt *= F;
                            if(rng.rand() < D) {real_collision = true;}
                            else {
                                wgt *= -1.0;
                                #pragma omp atomic
                                sign_change += 1.0;
                            }

                        } else { // Delta tracking branch
                            double P_real = E_tot/ Esmp;
                            if(rng.rand() < P_real) {real_collision = true;}
                        }

                        Esmp = 3.0*E_tot;
                        if(Esmp < 1.0) Esmp = 1.0;

                        if(real_collision) {
                            // Record collision
                            score_real_collision(wgt, x);

                            // Implicit caputure
                            wgt *= 1.0 - P_abs;

                            // Scatter
                            if(rng.rand() > P_straight_ahead) u *= -1;

                            // Russian Roulette
                            double xi = rng.rand();
                            roulette(wgt,alive,xi);
                            if(alive == false) {
                                #pragma omp atomic
                                cnts_sum += cnt;
                            }
                           
                        }// End real coll.
                    }

                    // Split if needed
                    if(alive and (std::abs(wgt) >= wgt_split)) {
                        double n_new = std::round(std::abs(wgt));
                        wgt /= n_new;
                        for(int j = 0; j < static_cast<int>(n_new-1); j++) {
                            #pragma omp critical
                            {
                                Split_Positions.push_back(x);
                                Split_Direction.push_back(u);
                                Split_Wgts.push_back(wgt);
                                Split_Esmps.push_back(Esmp);
                            }
                        }
                    }// split

                }// While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();
        }// Parallel

        n_particles = static_cast<int>(Split_Positions.size());
        Positions = Split_Positions;
        Directions = Split_Direction;
        Wgts = Split_Wgts;
        Esmps = Split_Esmps;
        
        Split_Positions.clear();
        Split_Direction.clear();
        Split_Wgts.clear();
        Split_Esmps.clear();
    }// While still particles

    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Old_Previous_XS_Bomb_Transport(std::unique_ptr<XS> const &xs) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Previous XS Bomb Paper Transport\n";
    int cnts_sum = 0;
    double sign_change = 0.0;

    // Particle bank vectors
    std::vector<double> Positions;
    std::vector<double> Directions;
    std::vector<double> Wgts;
    std::vector<double> Esmps;
    int n_particles = NPART;
    // Generate initial source
    for(int i = 0; i < n_particles; i++) {
        Positions.push_back(0.0);
        Directions.push_back(1.0);
        Wgts.push_back(1.0);
        Esmps.push_back(15.0*xs->Et(0.0));
    }
    std::vector<double> Split_Positions;
    std::vector<double> Split_Direction;
    std::vector<double> Split_Wgts;
    std::vector<double> Split_Esmps;

    while(n_particles > 0) {
        std::cout << " NParticles = " << n_particles << "\n";
        #pragma omp parallel
        {
            PCG rng;
            int thread_id = omp_get_thread_num();
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                double Esmp = Esmps[n];
                if(Esmp < 1.0) Esmp = 1.0;
                double x = Positions[n];
                double u = Directions[n];
                bool alive = true;
                int cnt = 1;
                double wgt = Wgts[n];
                bool real_collision = false;

                while(alive) {
                    real_collision = false;
                    double d = -std::log(rng.rand())/Esmp;
                    x += u*d;
                    
                    // Fist check for leak
                    if((x >= 2.0) or (x <= 0.0)) {
                        score_escape(wgt);
                        alive = false;
                        #pragma omp atomic
                        cnts_sum += cnt;
                    } else {
                        double E_tot = xs->Et(x);
                        cnt += 1;
                        if(E_tot > Esmp) { // First negative branch
                            double D = E_tot / (2*E_tot - Esmp);
                            double F = E_tot / (D*Esmp);
                            wgt *= F;
                            if(rng.rand() < D) {real_collision = true;}
                            else {
                                wgt *= -1.0;
                                #pragma omp atomic
                                sign_change += 1.0;
                            }

                        } else { // Delta tracking branch
                            double P_real = E_tot/ Esmp;
                            if(rng.rand() < P_real) {real_collision = true;}
                        }

                        Esmp = 15.0*E_tot; // Update sampling XS to current XS at location
                        if(Esmp < 1.0) Esmp = 1.0;

                        if(real_collision) {
                            score_real_collision(wgt,x);

                            // Implicit capture
                            wgt *= 1.0 - P_abs;

                            // Scatter
                            if(rng.rand() > P_straight_ahead) u *= -1;

                            // Russian Roulette
                            double xi = rng.rand();
                            roulette(wgt,alive,xi);
                            if(alive == false) {
                                #pragma omp atomic
                                cnts_sum += cnt;
                            }
                            
                        }

                        // Split if needed
                        if(alive and (std::abs(wgt) >= 2.0)) {
                            double n_new = std::round(std::abs(wgt));
                            wgt /= n_new;
                            for(int j = 0;j<static_cast<int>(n_new-1);j++) {
                                #pragma omp critical
                                {
                                    Split_Positions.push_back(x);
                                    Split_Direction.push_back(u);
                                    Split_Wgts.push_back(wgt);
                                    Split_Esmps.push_back(Esmp);
                                }
                            }
                        }// split
                    }
                }// While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();
        }// Parallel

        n_particles = static_cast<int>(Split_Positions.size());
        Positions = Split_Positions;
        Directions = Split_Direction;
        Wgts = Split_Wgts;
        Esmps = Split_Esmps;
        
        Split_Positions.clear();
        Split_Direction.clear();
        Split_Wgts.clear();
        Split_Esmps.clear();

    } // While split particles
    
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Output() {
    double collide_avg = collide / (double)NPART;
    double collide_sqr_avg = collide_sqr / (double)NPART;
    double collide_std = std::sqrt(std::abs(collide_avg*collide_avg - 
                          collide_sqr_avg)/((double)NPART - 1.0));

    avg_coll_rates_real.push_back(collide_avg);

    double all_collide_avg = all_collide / (double)NPART;
    double all_collide_sqr_avg = all_collide_sqr / (double)NPART;
    double all_collide_std = std::sqrt(std::abs(all_collide_avg*all_collide_avg - 
                          all_collide_sqr_avg)/((double)NPART - 1.0));

    avg_coll_rates_all.push_back(all_collide_avg);

    double escape_avg = escape / (double)NPART;
    double escape_sqr_avg = escape_sqr / (double)NPART;
    double escape_std = std::sqrt(std::abs(escape_avg*escape_avg - 
                           escape_sqr_avg)/((double)NPART - 1.0));

    double avg_xs_evals = xs_evals / (double)NPART;
    double avg_sgn_chngs = wgt_chngs / (double)NPART;

    double coll_rel_error = collide_std / collide_avg;
    double all_coll_rel_error = all_collide_std / all_collide_avg;
    double escape_rel_error = escape_std / escape_avg;
    double FOM_col = 1.0 / (avg_xs_evals * coll_rel_error * coll_rel_error);
    double FOM_all_coll = 1.0 / (avg_xs_evals * all_coll_rel_error * all_coll_rel_error);
    double FOM_esc = 1.0 / (avg_xs_evals*escape_rel_error*escape_rel_error);
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << " Colsn. Rate: " << collide_avg << " +/- " << collide_std;
    std::cout << ", Trans. Rate: " << escape_avg << " +/- " << escape_std << "\n";
    std::cout << " All Colsn. Rate: " << all_collide_avg << " +/- " << all_collide_std << "\n";
    std::cout << " Avg XS/Integration Evals: " << avg_xs_evals;
    std::cout << ", Avg Sign Changes: " << avg_sgn_chngs << "\n";
    std::cout << std::scientific;
    std::cout << " FOM_coll = " << FOM_col << ", FOM_all_coll = " << FOM_all_coll;
    std::cout << ", FOM_escp = " << FOM_esc << "\n\n";
}

std::unique_ptr<XS> make_cross_section(int type)
{
  // Determine cross section type for run
  if(type == -1) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Step\n\n";
    return std::make_unique<Step>();
  }
  else if(type == 0) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Constant\n\n";
    return std::make_unique<Constant>();
  }
  else if(type == 1) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Linearly Increasing\n\n";
    return std::make_unique<Lin_Increase>();
  }
  else if(type == 2) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Linearly Decreasing\n\n";
    return std::make_unique<Lin_Decrease>();
  }
  else if(type == 4) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Exponentially Decreasing\n\n";
    return std::make_unique<Exp_Decrease>();
  }
  else if(type == 3) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Exponentially Increasing\n\n";
    return std::make_unique<Exp_Increase>();
  }
  else if(type == 5) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Sharp Gaussian\n\n";
    return std::make_unique<Gauss_Sharp>();
  }
  else if(type == 6) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Broad Gaussian\n\n";
    return std::make_unique<Gauss_Broad>();
  }
  else {
    exit(1);
  }

}

void Zero_Values() {
    collide = 0.0;
    collide_sqr = 0.0;
    all_collide = 0.0;
    all_collide_sqr = 0.0;
    escape = 0.0;
    escape_sqr = 0.0;
    xs_evals = 0.0;
    wgt_chngs = 0.0;
    for(int i = 0; i < NFOMBINS; i++) {
        coll_density[i][0] = 0.0;
        coll_density[i][1] = 0.0;
        coll_density[i][2] = 0.0;
        coll_density[i][3] = 0.0;
        all_coll_density[i][0] = 0.0;
        all_coll_density[i][1] = 0.0;
        all_coll_density[i][2] = 0.0;
        all_coll_density[i][3] = 0.0;
    }
}

void RNG_Seeds() {
    int n_threads = omp_get_max_threads();
    for(int i = 0; i < n_threads; i++) {
        uint64_t seed = i+1;
        pcg_seeds.push_back(seed);
    }
}

int main() {
    RNG_Seeds();

    std::cout << "\n NParticles = " << NPART << ", NBins = " << NBIN << "\n\n";

    // Create and zero coll_density array
    for(int i = 0; i < NFOMBINS; i++) {
        std::vector<double> bin;
        bin.push_back(0.0);
        bin.push_back(0.0);
        bin.push_back(0.0);
        bin.push_back(0.0);
        coll_density.push_back(bin);
        all_coll_density.push_back(bin);
    }
    int xs_type = 1;
    for(int type = 1; type <= NTRIALS; type++) {
      std::unique_ptr<XS> crs = make_cross_section(xs_type);
      std::cout << " " << type << "\n";

      Zero_Values();
      Delta_Tracking(crs);
      Output();

      /*Zero_Values();
      Meshed_Delta_Tracking(crs);
      Output();

      Zero_Values();
      Negative_Weight_Delta_Tracking(crs);
      Output();

      Zero_Values();
      Meshed_Negative_Weight_Delta_Tracking(crs);
      Output();

      Zero_Values();
      Bomb_Transport(crs,0.80);
      Output();

      Zero_Values();
      Meshed_Bomb_Transport(crs,0.80);
      Output();
     
      Zero_Values();
      Improving_Meshed_Bomb_Transport(crs,0.80);
      Output();*/
      
      /*Zero_Values();
      Previous_XS_Bomb_Transport(crs);
      Output();*/

    }

    // Get avge of vectors and std
    double real_coll_sum = 0.0;
    double real_coll_sqr_sum = 0.0;
    double all_coll_sum = 0.0;
    double all_coll_sqr_sum = 0.0;
    for(int i = 0; i < static_cast<int>(avg_coll_rates_real.size()); i++) {
        real_coll_sum += avg_coll_rates_real[i];
        real_coll_sqr_sum += avg_coll_rates_real[i]*avg_coll_rates_real[i];
        all_coll_sum += avg_coll_rates_all[i];
        all_coll_sqr_sum += avg_coll_rates_all[i]*avg_coll_rates_all[i];
    }
    double n_trials = static_cast<double>(NTRIALS);
    double avg_coll_real = real_coll_sum / n_trials;
    double avg_coll_all = all_coll_sum / n_trials;
    double avg_coll_sqr_real = real_coll_sqr_sum / n_trials;
    double avg_coll_sqr_all = all_coll_sqr_sum / n_trials;
    double std_avg_real = std::sqrt(std::abs(avg_coll_real*avg_coll_real - avg_coll_sqr_real)/(n_trials - 1.0));
    double std_avg_all = std::sqrt(std::abs(avg_coll_all*avg_coll_all - avg_coll_sqr_all)/(n_trials - 1.0));
    std::cout << "\n\n RESUTLS:\n";
    std::cout << " Avg Coll Rate Real = " << avg_coll_real << " +/- " << std_avg_real << "\n";
    std::cout << " Avg Coll Rate All  = " << avg_coll_all <<  " +/- " << std_avg_all << "\n";

    pcg_seeds.clear();
    return 0;
}
