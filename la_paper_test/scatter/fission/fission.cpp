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
const int NGENS = 1;
const int NBIN = 5;
const int NFOMBINS = 100;
const double Fdx = 2.0/(double)NFOMBINS;
const double dx = 2.0/(double)NBIN;
const double p = 0.75; // Esamp / Emaj for full system
const double p_mshd = 0.75; // Esamp[i] / Emaj[i] for bin
const double q = 0.3;
const double q_mshd = 0.3;

const double P_abs = 0.3;
const double P_fis = 0.4;
const double nu    = 2.5;
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
const int NTHREADS = 40;
std::vector<uint64_t> pcg_seeds;

// Outputfile
std::ofstream File;

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

        virtual ~XS() = default;

        // All virtual methods
        virtual double Et(double /*x*/) {return -1.0;}

        double P_nc() {return Pnc;}

        // Data
        double Pnc; // Prob. of no collision
        double G; // 1 - Pnc (Prob. of collision)
        double Emax; // Maj over region
        double Em[NBIN]; // Maj in each bin
        double Em_imp[NBIN]; // p*Em[bin] for IBT
        double Esmp[NBIN]; // Sampling XS in each bin for NWDT and MNWDT
};

// Constant XS (Et = 1.0 ∀ 0<x<2, or Et=0)
class Constant : public XS {
    public:
        Constant():XS(1.0) {
            for(int b = 0; b < NBIN; b++) {
                Em[b] = 1.0;
                Em_imp[b] = Em[b];
                Esmp[b] = p_mshd;
            }
        }

        ~Constant() = default;

        double Et(double /*x*/) {return 1.0;}

}; // Constants

// Step XS (Et = 1.5 ∀ 0<x<0.278), Et=1.0 ∀ 0.278<x<2)
class Step : public XS {
    public:
        Step():XS(1.5) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
                Em_imp[b] = Em[b];
                Esmp[b] = p_mshd*Et(x);
            }
        }
        ~Step() = default;

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
                Em_imp[b] = Em[b];
                Esmp[b] = p_mshd*Et(x);
            }    
        }
        ~Lin_Decrease() = default;

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
                Em_imp[b] = Em[b];
                Esmp[b] = p_mshd*Et(x);
            }
        }
        ~Lin_Increase() = default;

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
                Em_imp[b] = Em[b];
                Esmp[b] = p_mshd*Et(x);
            }
        }
        ~Exp_Decrease() = default;

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
                Em_imp[b] = Em[b];
                Esmp[b] = p_mshd*Et(x);
            }
        }
        ~Exp_Increase() = default;

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
                Em_imp[b] = Em[b];
                Esmp[b] = p_mshd*Et(x);
            }
        }
        ~Gauss_Sharp() = default;

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
                Em_imp[b] = Em[b];
                Esmp[b] = p_mshd*Et(x);
            }
        }
        ~Gauss_Broad() = default;

        double Et(double x) {
            return A*std::exp(-a_b*(x-z)*(x-z));
        }
};

struct Particle {
    Particle(double x_=0.0, double u_=1.0, double wgt_=1.0) {
        x = x_;
        bin = std::floor(x / dx);
        u = u_;
        wgt = wgt_;
        alive = true;
    }

    void kill() {alive = false;}
    void move(double dist) {
        x += u*dist;
        bin = std::floor(x / dx);
    }
    void turn() {u *= -1.0;}
    void xs_eval() {xs_evals_cnt++;}
    void virt_coll() {virtual_colls++;}
    void bin_crs() {bin_cross++;}

    double x;
    int bin;
    double u;
    double wgt;
    double Esmp;
    bool alive;
    int xs_evals_cnt = 0;
    int virtual_colls = 0;
    int bin_cross = 0;
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

void Delta_Tracking(std::unique_ptr<XS> const &xs,
        std::vector<Particle> &bank) {
    std::cout << "\n Delta Tracking\n";

    int cnts_sum = 0;

    #pragma omp parallel
    {    
        PCG rng;
        int thread_id;
        #ifdef _OPENMP
        thread_id = omp_get_thread_num();
        #else
        thread_id = 0;
        #endif
        uint64_t pcg_seed;
        #pragma omp atomic read
        pcg_seed = pcg_seeds[thread_id];
        rng.seed(pcg_seed);

        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            bool virtual_collision;
            Particle p = bank[n];
            while(p.alive) {
                virtual_collision = true;
                while(virtual_collision) {
                    double d = -std::log(rng.rand())/(xs->Emax);
                    p.move(d);
                    if((p.x >= 2.0) or (p.x <= 0.0)) {
                        score_escape(p.wgt); 
                        virtual_collision = false;
                        p.kill();
                        #pragma omp atomic
                        cnts_sum += p.xs_evals_cnt;
                    } else {
                        p.xs_eval();
                        if(rng.rand() < (xs->Et(p.x)/xs->Emax)) {
                            virtual_collision = false;
                        }
                        // Score every collision
                        double score = p.wgt*xs->Et(p.x) / xs->Emax;
                        score_all_collision(score,p.x);
                    }
                }
                
                if(p.alive) {
                    // Score real collision
                    score_real_collision(p.wgt,p.x); 

                    // Implicit capture
                    p.wgt *= 1.0 - P_abs; // Equivalent to 1.0 - (Ea/Et)

                    // Scatter
                    if(rng.rand() > P_straight_ahead) p.turn();

                    // Russian Roulette
                    double xi = rng.rand();
                    roulette(p.wgt, p.alive, xi);
                    if(p.alive == false) {
                        #pragma omp atomic
                        cnts_sum += p.xs_evals_cnt;
                    }
                } // If alive for real collision
            } // While alive
        } // For all particles
        #pragma omp critical
        {
            pcg_seeds[thread_id] = rng.get_seed();
        }
    } // Parallel
    xs_evals += cnts_sum;
}

void Meshed_Delta_Tracking(std::unique_ptr<XS> const &xs,
        std::vector<Particle> const &bank) {
    std::cout << "\n Meshed Delta Tracking\n";

    int cnts_sum = 0;
    int virtual_cnt_sum = 0;
    int bin_cnt_sum = 0;

    #pragma omp parallel
    {    
        PCG rng;
        int thread_id;
        #ifdef _OPENMP
        thread_id = omp_get_thread_num();
        #else
        thread_id = 0;
        #endif
        uint64_t pcg_seed;
        #pragma omp atomic read
        pcg_seed = pcg_seeds[thread_id];
        rng.seed(pcg_seed);

        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            bool virtual_collision = true;
            Particle p = bank[n];
            double Emax = xs->Em[p.bin];
            double d,d_bin;
            while(p.alive) {
                virtual_collision = true;
                while(virtual_collision) {
                    if(p.u == -1.0) d_bin = p.x-static_cast<double>(p.bin)*dx;
                    else d_bin = (static_cast<double>(p.bin)*dx + dx) - p.x;
                    d = -std::log(rng.rand())/Emax;
                    if(d_bin < d) {
                        p.bin_crs();
                        d = d_bin + 1e-6;
                        p.move(d);
                        if((p.x >= 2.0) or (p.x  <= 0.0)) {
                            score_escape(p.wgt);
                            virtual_collision = false;
                            p.kill();
                            
                            // Score other tallies 
                            #pragma omp atomic
                            cnts_sum += p.xs_evals_cnt;
                            #pragma omp atomic
                            virtual_cnt_sum += p.virtual_colls;
                            #pragma omp atomic
                            bin_cnt_sum += p.bin_cross;
                        }
                        else {
                            Emax = xs->Em[p.bin];
                        }
                    }
                    else {
                        p.xs_eval();
                        p.move(d);
                        double xi = rng.rand();
                        double Pr = xs->Et(p.x)/Emax;
                        if(xi < Pr) {
                            virtual_collision = false;
                        } else {p.virt_coll();}
                        
                        double score = p.wgt*xs->Et(p.x) /Emax;
                        score_all_collision(score, p.x);
                    }
                }// While virtual

                if(p.alive) {
                    // Score real collision
                    score_real_collision(p.wgt,p.x);
                   
                    // Implicit capture
                    p.wgt *= 1.0 - P_abs;
                    
                    // Scatter
                    if(rng.rand() > P_straight_ahead) p.turn();

                    // Russian Roulette
                    double xi = rng.rand();
                    roulette(p.wgt, p.alive, xi);
                    if(p.alive == false){
                        #pragma omp atomic
                        cnts_sum += p.xs_evals_cnt;
                        #pragma omp atomic
                        virtual_cnt_sum += p.virtual_colls;
                        #pragma omp atomic
                        bin_cnt_sum += p.bin_cross;
                    }
                }
            } // While alive
        }// For all particles
        #pragma omp atomic write
        pcg_seeds[thread_id] = rng.get_seed();
    }// Parallel
    xs_evals += cnts_sum;
}

void Negative_Weight_Delta_Tracking(std::unique_ptr<XS> const &xs,
        std::vector<Particle> const &bank) {
    std::cout << "\n Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    double sign_change = 0.0;

   std::vector<Particle> Bank = bank; // Copy for re-writing after splits
   int n_particles = static_cast<int>(Bank.size());

    while(n_particles > 0) {
        //std::cout << " nparticles = " << n_particles << "\n";
        #pragma omp parallel
        {    
            PCG rng;
            int thread_id;
            #ifdef _OPENMP
            thread_id = omp_get_thread_num();
            #else
            thread_id = 0;
            #endif
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            std::vector<Particle> this_thread_Splits;

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                double Esamp = p*(xs->Emax);
                Particle p = Bank[n];
                while(p.alive) {
                    double d = -std::log(rng.rand())/Esamp;
                    p.move(d);
                    
                    if((p.x >= 2.0) or (p.x <= 0.0)) {
                        score_escape(p.wgt);

                        p.kill();
                        #pragma omp atomic
                        cnts_sum += p.xs_evals_cnt;
                    }
                    else if(rng.rand() < q) {
                        // Collision is real, addjust weight
                        p.wgt *= xs->Et(p.x)/(Esamp * q);
                        double score = p.wgt *q;
                        score_all_collision(score,p.x);
                        score_real_collision(p.wgt,p.x);
                        p.xs_eval();
                        
                        // Implicit capture
                        p.wgt *= 1.0 - P_abs;
                        
                        // Scatter
                        if(rng.rand() > P_straight_ahead) p.turn();

                        // Russian Roulette
                        double xi = rng.rand();
                        roulette(p.wgt,p.alive,xi);
                        if(p.alive == false) {
                            #pragma omp atomic
                            cnts_sum += p.xs_evals_cnt;
                        }
                        
                    }
                    else { // Collision is virtual
                        double dw = ((1.0 - (xs->Et(p.x)/Esamp))/(1.0 - q));
                        if(dw < 0.0) {
                            #pragma omp atomic
                            sign_change += 1.0;
                        }
                        double score = p.wgt*(xs->Et(p.x)/Esamp);
                        score_all_collision(score,p.x);
                        p.wgt = p.wgt*dw;
                        p.xs_eval();
                    }

                    // Split if needed
                    if(p.alive and (std::abs(p.wgt) >= wgt_split)) {
                        double n_new = std::round(std::abs(p.wgt));
                        p.wgt /= n_new;
                        for(int j = 0; j < static_cast<int>(n_new-1); j++) {
                            Particle p_daughter(p.x,p.u,p.wgt);
                            p_daughter.xs_evals_cnt = p.xs_evals_cnt;
                            p_daughter.virtual_colls = p.virtual_colls;
                            p_daughter.bin = p.bin;
                            this_thread_Splits.push_back(p_daughter);
                        }
                    }// split
                } // While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();

            // Clear Bank to accept splits
            #pragma omp single
            {
                Bank.clear();
            }
            // Add thread splits to Bank
            #pragma omp barrier
            #pragma omp critical
            {
                Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                        std::end(this_thread_Splits));
            }
        }// Parallel
        n_particles = static_cast<int>(Bank.size());

    } // While still split particles
    
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Meshed_Negative_Weight_Delta_Tracking(std::unique_ptr<XS> const &xs,
        std::vector<Particle> const &bank) {
    std::cout << "\n Meshed Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    int bin_cnt_sum = 0;
    double sign_change = 0.0;

    std::vector<Particle> Bank = bank;
    int n_particles = static_cast<int>(Bank.size()); 

    while(n_particles > 0) {
        //std::cout << " nparticles = " << n_particles << "\n";
        #pragma omp parallel
        {    
            PCG rng;
            int thread_id;
            #ifdef _OPENMP
            thread_id = omp_get_thread_num();
            #else
            thread_id = 0;
            #endif
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            std::vector<Particle> this_thread_Splits;

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                Particle p = Bank[n];
                double Esamp = xs->Esmp[p.bin];
                double d,d_bin;
                while(p.alive) {
                    if(p.u == -1.0) d_bin = p.x-static_cast<double>(p.bin)*dx;
                    else d_bin = (static_cast<double>(p.bin)*dx + dx) - p.x;
                    d = -std::log(rng.rand())/Esamp;
                    if(d_bin < d) {
                        p.bin_crs();
                        d = d_bin + 1e-6;
                        p.move(d);
                        if((p.x >= 2.0) or (p.x <= 0.0)) {
                            score_escape(p.wgt);
                            p.kill();
                            
                            #pragma omp atomic
                            cnts_sum += p.xs_evals_cnt;
                            #pragma omp atomic
                            bin_cnt_sum += p.bin_cross;
                        }
                        else {
                            Esamp = xs->Esmp[p.bin];
                        }
                    }
                    else {
                        p.move(d);
                        if(rng.rand() < q_mshd) { // Real collision
                            // update weight
                            p.wgt *= (xs->Et(p.x)/(Esamp*q_mshd));
                            double score = p.wgt*(q_mshd);
                            score_all_collision(score,p.x);
                            score_real_collision(p.wgt,p.x);
                            p.xs_eval();

                            // Implicit capture
                            p.wgt *= 1.0 - P_abs;

                            // Scatter
                            if(rng.rand() > P_straight_ahead) p.turn();

                            // Russian Roulette
                            double xi = rng.rand();
                            roulette(p.wgt,p.alive,xi);
                            if(p.alive == false) {
                                #pragma omp atomic
                                cnts_sum += p.xs_evals_cnt;
                                #pragma omp atomic
                                bin_cnt_sum += p.bin_cross;
                            }
                            
                        }
                        else {
                            p.xs_eval();
                            double dw=(1.0-(xs->Et(p.x)/Esamp))/(1.0-q_mshd);
                            if(dw < 0.0) {
                                #pragma omp atomic
                                sign_change += 1.0;
                            }
                            double score = p.wgt*(xs->Et(p.x)/Esamp);
                            score_all_collision(score,p.x);
                            p.wgt = p.wgt*dw;
                        }

                        // Split if needed
                        if(p.alive and (std::abs(p.wgt) >= wgt_split)) {
                            double n_new = std::round(std::abs(p.wgt));
                            p.wgt /= n_new;
                            for(int j=0; j < static_cast<int>(n_new-1); j++) {
                                Particle p_daughter(p.x,p.u,p.wgt);
                                p_daughter.xs_evals_cnt = p.xs_evals_cnt;
                                p_daughter.bin_cross = p.bin_cross;
                                p_daughter.bin = p.bin;
                                this_thread_Splits.push_back(p_daughter);
                            }
                        }// split
                    }
                }// While alive
            } // For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();

            // Clear bank for new particles
            #pragma omp single
            {
                Bank.clear();
            }

            // Add thread splits to Bank
            #pragma omp barrier
            #pragma omp critical
            {
                Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                        std::end(this_thread_Splits));
            }
        }// Parallel
        
        n_particles = static_cast<int>(Bank.size());

    }// While still split particles

    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Bomb_Transport(std::unique_ptr<XS> const &xs, double P,
        std::vector<Particle> const &bank) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Bomb Paper Transport, p = " << P << "\n";
    int cnts_sum = 0;
    double sign_change = 0.0;
    
    // Particle bank vectors
    std::vector<Particle> Bank = bank;
    int n_particles = static_cast<int>(Bank.size());
    
    while(n_particles > 0) {
        //std::cout << " nparticles = " << n_particles << "\n";
        #pragma omp parallel
        {
            PCG rng;
            int thread_id;
            #ifdef _OPENMP
            thread_id = omp_get_thread_num();
            #else
            thread_id = 0;
            #endif
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            std::vector<Particle> this_thread_Splits;

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                Particle p = Bank[n];
                double Esmp = P*xs->Emax;
                bool real_collision = false;

                while(p.alive) {
                    double d = -std::log(rng.rand())/Esmp;
                    p.move(d);
                    real_collision = false;
                    
                    // Fist check for leak
                    if((p.x >= 2.0) or (p.x <= 0.0)) {
                        score_escape(p.wgt);
                        p.kill();
                        #pragma omp atomic
                        cnts_sum += p.xs_evals_cnt;
                    } else {
                        double E_tot = xs->Et(p.x);
                        p.xs_eval();
                        if(E_tot > Esmp) { // First negative branch
                           //double D_alpha = alpha*E_tot / ((2.0 + alpha)*E_tot - Esmp);
                            double D = E_tot / (2*E_tot - Esmp);
                            double F = (E_tot / (D*Esmp));
                            double score = p.wgt*E_tot/Esmp;
                            score_all_collision(score,p.x);

                            //if(rand(rng) < D_alpha) {
                            if(rng.rand() < D) {
                                real_collision = true;
                                p.wgt *= F;
                                //wgt *=  F*(D/D_alpha);
                            }
                            else {
                                p.wgt *= -F;
                                //wgt *= -F*((1. - D)/(1. - D_alpha));
                                #pragma omp atomic
                                sign_change += 1.0;
                            }
                            

                        } else { // Delta tracking branch
                            double P_real = E_tot/ Esmp;
                            if(rng.rand() < P_real) {real_collision = true;}
                            double score = p.wgt*E_tot/Esmp;
                            score_all_collision(score,p.x);
                        }

                        if(real_collision) {
                            // Score real collision
                            score_real_collision(p.wgt,p.x);

                            // Implicit caputure
                            p.wgt *= 1.0 - P_abs;

                            // Scatter
                            if(rng.rand() > P_straight_ahead) p.turn();

                            // Russian Roulette
                            double xi = rng.rand();
                            roulette(p.wgt,p.alive,xi);
                            if(p.alive == false) {
                                #pragma omp atomic
                                cnts_sum += p.xs_evals_cnt;
                            }
                            
                        }// End real coll.
                    }

                    // Split if needed
                    if(p.alive and (std::abs(p.wgt) >= wgt_split)) {
                        double n_new = std::floor(std::abs(p.wgt));
                        p.wgt /= n_new;
                        for(int j = 0; j < static_cast<int>(n_new-1); j++) {
                            Particle p_daughter(p.x,p.u,p.wgt);
                            p_daughter.xs_evals_cnt = p.xs_evals_cnt;
                            p_daughter.bin_cross = p.bin_cross;
                            p_daughter.bin = p.bin;
                            this_thread_Splits.push_back(p_daughter);
                        }
                    }// split

                }// While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();

            #pragma omp single
            {
                Bank.clear();
            }
            
            #pragma omp barrier
            #pragma omp critical
            {
                Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                        std::end(this_thread_Splits));
            }
        }// Parallel

        n_particles = static_cast<int>(Bank.size());

    }// While still particles

    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Meshed_Bomb_Transport(std::unique_ptr<XS> const &xs, double P,
        std::vector<Particle> const &bank) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Meshed Bomb Paper Transport, p = " << P << "\n";
    int cnts_sum = 0;
    double sign_change = 0.0;

    // Particle bank vectors
    std::vector<Particle> Bank = bank;
    int n_particles = static_cast<int>(Bank.size());
    
    while(n_particles > 0) {
        //std::cout << " nparticles = " << n_particles << "\n";
        #pragma omp parallel
        {
            PCG rng;
            int thread_id;
            #ifdef _OPENMP
            thread_id = omp_get_thread_num();
            #else
            thread_id = 0;
            #endif
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];

            std::vector<Particle> this_thread_Splits;

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                Particle p = Bank[n];
                double d_bin;
                double Esmp = P*xs->Em[p.bin];
                bool real_collision = false;

                while(p.alive) {
                    double d = -std::log(rng.rand())/Esmp;
                    if(p.u == -1.0) d_bin = p.x-static_cast<double>(p.bin)*dx;
                    else d_bin = (static_cast<double>(p.bin)*dx + dx) - p.x;
                    real_collision = false;
                    
                    if(d_bin < d) {
                        p.move(d_bin+1e-6);
                        if((p.x >= 2.0) or (p.x <= 0.0)) {
                            p.kill();
                            score_escape(p.wgt);
                            #pragma omp atomic
                            cnts_sum += p.xs_evals_cnt;
                        } else {
                            Esmp = P*xs->Em[p.bin];
                        }

                    } else {
                        p.move(d);
                        // Fist check for leak
                        if((p.x >= 2.0) or (p.x <= 0.0)) {
                            p.kill();
                            score_escape(p.wgt);
                            #pragma omp atomic
                            cnts_sum += p.xs_evals_cnt;
                        } else {
                            double E_tot = xs->Et(p.x);
                            p.xs_eval();
                            if(E_tot > Esmp) { // First negative branch
                                double D = E_tot / (2*E_tot - Esmp);
                                double F = E_tot / (D*Esmp);
                                double score = p.wgt*E_tot/Esmp;
                                score_all_collision(score,p.x);
                                p.wgt *= F;
                                if(rng.rand() < D) {real_collision = true;}
                                else {
                                    p.wgt *= -1.0;
                                    #pragma omp atomic
                                    sign_change += 1.0;
                                }

                            } else { // Delta tracking branch
                                double P_real = E_tot/ Esmp;
                                double score = p.wgt*E_tot/Esmp;
                                score_all_collision(score,p.x);
                                if(rng.rand() < P_real) {real_collision=true;}
                            }

                            if(real_collision) {
                                // Score real collision
                                score_real_collision(p.wgt,p.x);

                                // Implicit capture
                                p.wgt *= 1.0 - P_abs;

                                // Scatter
                                if(rng.rand() > P_straight_ahead) p.turn();

                                // Russian Roulette
                                double xi = rng.rand();
                                roulette(p.wgt,p.alive,xi);
                                if(p.alive == false) {
                                    #pragma omp atomic
                                    cnts_sum += p.xs_evals_cnt;
                                }

                            }

                            // Split if needed
                            if(p.alive and (std::abs(p.wgt) >= wgt_split)) {
                                double n_new = std::round(std::abs(p.wgt));
                                p.wgt /= n_new;
                                for(int j=0;j<static_cast<int>(n_new-1);j++) {
                                    Particle p_daughter(p.x,p.u,p.wgt);
                                    p_daughter.xs_evals_cnt = p.xs_evals_cnt;
                                    p_daughter.bin_cross = p.bin_cross;
                                    p_daughter.bin = p.bin;
                                    this_thread_Splits.push_back(p_daughter);
                                }
                            }// split
                        }
                    }
                }// While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = pcg_seed;

            #pragma omp single
            {
                Bank.clear();
            }
            
            #pragma omp barrier
            #pragma omp critical
            {
                Bank.insert(std::end(Bank),std::begin(this_thread_Splits),
                        std::end(this_thread_Splits));
            }
        }// Parallel

        n_particles = static_cast<int>(Bank.size());

    } // While split particles
    
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Improving_Meshed_Bomb_Transport(std::unique_ptr<XS> const &xs, double P,
        std::vector<Particle> const &bank) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Improving Meshed Bomb Paper Transport, p = " << P << "\n";
    int cnts_sum = 0;
    double sign_change = 0.0;

    // Particle bank vectors
    std::vector<Particle> Bank = bank;
    int n_particles = static_cast<int>(Bank.size());
    
    while(n_particles > 0) {
        #pragma omp parallel
        {
            PCG rng;
            int thread_id;
            #ifdef _OPENMP
            thread_id = omp_get_thread_num();
            #else
            thread_id = 0;
            #endif
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            std::vector<Particle> this_thread_Splits;

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                Particle p = Bank[n];
                double d_bin;
                double Esmp = xs->Em_imp[p.bin];
                bool real_collision = false;

                while(p.alive) {
                    double d = -std::log(rng.rand())/Esmp;
                    if(p.u == -1.0) d_bin = p.x-static_cast<double>(p.bin)*dx;
                    else d_bin = (static_cast<double>(p.bin)*dx + dx) - p.x;
                    real_collision = false;
                    
                    if(d_bin < d) {
                        p.move(d_bin+1e-6);
                        if((p.x >= 2.0) or (p.x <= 0.0)) {
                            p.kill();
                            score_escape(p.wgt);
                            #pragma omp atomic
                            cnts_sum += p.xs_evals_cnt;
                        } else {
                            Esmp = xs->Em_imp[p.bin];
                        }

                    } else {
                        p.move(d);
                        // Fist check for leak
                        if((p.x >= 2.0) or (p.x <= 0.0)) {
                            p.kill();
                            score_escape(p.wgt);
                            #pragma omp atomic
                            cnts_sum += p.xs_evals_cnt;
                        } else {
                            double E_tot = xs->Et(p.x);
                            p.xs_eval();
                            if(E_tot > Esmp) { // First negative branch
                                double D = E_tot / (2*E_tot - Esmp);
                                double F = E_tot / (D*Esmp);
                                double score = p.wgt*E_tot/Esmp;
                                score_all_collision(score,p.x);
                                p.wgt *= F;
                                if(rng.rand() < D) {real_collision = true;}
                                else {
                                    p.wgt *= -1.0;
                                    #pragma omp atomic
                                    sign_change += 1.0;
                                }

                                

                            } else { // Delta tracking branch
                                double P_real = E_tot/ Esmp;
                                double score = p.wgt*E_tot/Esmp;
                                score_all_collision(score,p.x);
                                if(rng.rand() < P_real) {real_collision = true;}
                            }

                            if(real_collision) {
                                // Score real collision
                                score_real_collision(p.wgt,p.x);

                                // Implicit capture
                                p.wgt *= 1.0 - P_abs;

                                // Scatter
                                if(rng.rand() > P_straight_ahead) p.turn();

                                // Russian Roulette
                                double xi = rng.rand();
                                roulette(p.wgt,p.alive,xi);
                                if(p.alive == false) {
                                    #pragma omp atomic
                                    cnts_sum += p.xs_evals_cnt;
                                }

                            }
                            
                            if(E_tot > Esmp) {
                                // Improve Esmp
                                Esmp = E_tot;
                                #pragma omp atomic write
                                xs->Em_imp[p.bin] = Esmp;
                            }
                            

                            // Split if needed
                            if(p.alive and (std::abs(p.wgt) >= wgt_split)) {
                                double n_new = std::round(std::abs(p.wgt));
                                p.wgt /= n_new;
                                for(int j=0;j<static_cast<int>(n_new-1);j++) {
                                    Particle p_daughter(p.x,p.u,p.wgt);
                                    p_daughter.bin = p.bin;
                                    p_daughter.xs_evals_cnt = p.xs_evals_cnt;
                                    p_daughter.bin_cross = p.bin_cross;
                                    this_thread_Splits.push_back(p_daughter);
                                }
                            }// split
                        }
                    }
                }// While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();

            #pragma omp single
            {
                Bank.clear();
            }
            
            #pragma omp barrier
            #pragma omp critical
            {
                Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                        std::end(this_thread_Splits));
            }
        }// Parallel

        n_particles = static_cast<int>(Bank.size());

    } // While split particles

    
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Previous_XS_Bomb_Transport(std::unique_ptr<XS> const &xs,
        std::vector<Particle> const &bank) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Previous XS Bomb Paper Transport\n";
    int cnts_sum = 0;
    double sign_change = 0.0;
    
    // Particle bank vectors
    std::vector<Particle> Bank = bank;
    int n_particles = static_cast<int>(Bank.size());
    for(int i = 0; i < n_particles; i++) {
        Bank[i].Esmp = xs->Et(Bank[i].x);
    }

    while(n_particles > 0) {
        #pragma omp parallel
        {
            PCG rng;
            int thread_id;
            #ifdef _OPENMP
            thread_id = omp_get_thread_num();
            #else
            thread_id = 0;
            #endif
            uint64_t pcg_seed;
            #pragma omp atomic read
            pcg_seed = pcg_seeds[thread_id];
            rng.seed(pcg_seed);

            std::vector<Particle> this_thread_Splits;

            #pragma omp for
            for(int n = 0; n < n_particles; n++) {
                Particle p = Bank[n];
                double Esmp = p.Esmp;
                if(Esmp < 1.0) Esmp = 1.0;
                bool real_collision = false;

                while(p.alive) {
                    double d = -std::log(rng.rand())/Esmp;
                    p.move(d);
                    real_collision = false;
                    
                    // Fist check for leak
                    if((p.x >= 2.0) or (p.x <= 0.0)) {
                        score_escape(p.wgt);
                        p.kill();
                        #pragma omp atomic
                        cnts_sum += p.xs_evals_cnt;
                    } else {
                        double E_tot = xs->Et(p.x);
                        p.xs_eval();
                        if(E_tot > Esmp) { // First negative branch
                            double D = E_tot / (2*E_tot - Esmp);
                            double F = E_tot / (D*Esmp);
                            p.wgt *= F;
                            if(rng.rand() < D) {real_collision = true;}
                            else {
                                p.wgt *= -1.0;
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
                            score_real_collision(p.wgt, p.x);

                            // Implicit caputure
                            p.wgt *= 1.0 - P_abs;

                            // Scatter
                            if(rng.rand() > P_straight_ahead) p.turn();

                            // Russian Roulette
                            double xi = rng.rand();
                            roulette(p.wgt,p.alive,xi);
                            if(p.alive == false) {
                                #pragma omp atomic
                                cnts_sum += p.xs_evals_cnt;
                            }
                           
                        }// End real coll.
                    }

                    // Split if needed
                    if(p.alive and (std::abs(p.wgt) >= wgt_split)) {
                        double n_new = std::round(std::abs(p.wgt));
                        p.wgt /= n_new;
                        for(int j = 0; j < static_cast<int>(n_new-1); j++) {
                            Particle p_daughter(p.x,p.u,p.wgt);
                            p_daughter.bin = p.bin;
                            p_daughter.xs_evals_cnt = p.xs_evals_cnt;
                            p_daughter.bin_cross = p.bin_cross;
                            this_thread_Splits.push_back(p_daughter); 
                        }
                    }// split

                }// While alive
            }// For all particles
            #pragma omp atomic write
            pcg_seeds[thread_id] = rng.get_seed();
            
            #pragma omp barrier
            #pragma omp single
            {
                Bank.clear();
            }

            #pragma omp critical
            {
                Bank.insert(std::end(Bank), std::begin(this_thread_Splits),
                        std::end(this_thread_Splits));
            }
        }// Parallel

        n_particles = static_cast<int>(Bank.size());

    }// While still particles

    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Output() {
    double collide_avg = collide / (double)NPART;
    double collide_sqr_avg = collide_sqr / (double)NPART;
    double collide_std = std::sqrt(std::abs(collide_avg*collide_avg - 
                          collide_sqr_avg)/((double)NPART - 1.0));

    double all_collide_avg = all_collide / (double)NPART;
    double all_collide_sqr_avg = all_collide_sqr / (double)NPART;
    double all_collide_std = std::sqrt(std::abs(all_collide_avg*all_collide_avg - 
                          all_collide_sqr_avg)/((double)NPART - 1.0));

    double escape_avg = escape / (double)NPART;
    double escape_sqr_avg = escape_sqr / (double)NPART;
    double escape_std = std::sqrt(std::abs(escape_avg*escape_avg - 
                           escape_sqr_avg)/((double)NPART - 1.0));

    double avg_xs_evals = xs_evals / (double)NPART;
    double avg_sgn_chngs = wgt_chngs / (double)NPART;

    // Calculations and output for real collision density profile
    for(int i = 0; i < NFOMBINS; i++) {
        // Get avg for bin
        double coll_avg = coll_density[i][0] / static_cast<double>(NPART);
        double coll_sqr_avg = coll_density[i][1] / static_cast<double>(NPART);
        double coll_sig = std::sqrt(std::abs(coll_avg*coll_avg - coll_sqr_avg)
                /(static_cast<double>(NPART) - 1.0));
        coll_density[i][2] = coll_sig;
        double rel_error = coll_sig / coll_avg;
        coll_density[i][3] = 1.0 / (avg_xs_evals * rel_error * rel_error);

        // Output avg real coll desnity in bin
        if(i == 0) {File << coll_avg;}
        else {File << "," << coll_avg;}
    }
    File << "\n";
    // Output std of real coll density
    for(int i = 0; i < NFOMBINS; i++) {
        if(i == 0) {File << coll_density[i][2];}
        else {File << "," << coll_density[i][2];}
    }
    File << "\n";
    // Output FOM of real coll density
    for(int i = 0; i < NFOMBINS; i++) {
        if(i == 0) {File << coll_density[i][3];}
        else {File << "," << coll_density[i][3];}
    }
    File << "\n";

    // Calculations and output for all collision density profile
    for(int i = 0; i < NFOMBINS; i++) {
        // Get avg for bin
        double all_coll_avg = all_coll_density[i][0] / static_cast<double>(NPART);
        double all_coll_sqr_avg = all_coll_density[i][1] / static_cast<double>(NPART);
        double all_coll_sig = std::sqrt(std::abs(all_coll_avg*all_coll_avg - all_coll_sqr_avg)
                /(static_cast<double>(NPART) - 1.0));
        all_coll_density[i][2] = all_coll_sig;
        double all_rel_error = all_coll_sig / all_coll_avg;
        all_coll_density[i][3] = 1.0 / (avg_xs_evals * all_rel_error * all_rel_error);

        // Output avg all coll desnity in bin
        if(i == 0) {File << all_coll_avg;}
        else {File << "," << all_coll_avg;}
    }
    File << "\n";
    // Output std of all coll density
    for(int i = 0; i < NFOMBINS; i++) {
        if(i == 0) {File << all_coll_density[i][2];}
        else {File << "," << all_coll_density[i][2];}
    }
    File << "\n";
    // Output FOM of all coll density
    for(int i = 0; i < NFOMBINS; i++) {
        if(i == 0) {File << all_coll_density[i][3];}
        else {File << "," << all_coll_density[i][3];}
    }
    File << "\n\n";

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
    File << "#XS,S\n";
    return std::make_unique<Step>();
  }
  else if(type == 0) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Constant\n\n";
    File << "#XS,C\n";
    return std::make_unique<Constant>();
  }
  else if(type == 1) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Linearly Increasing\n\n";
    File << "#XS,LI\n";
    return std::make_unique<Lin_Increase>();
  }
  else if(type == 2) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Linearly Decreasing\n\n";
    File << "#XS,LD\n"; 
    return std::make_unique<Lin_Decrease>();
  }
  else if(type == 4) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Exponentially Decreasing\n\n";
    File << "#XS,ED\n";
    return std::make_unique<Exp_Decrease>();
  }
  else if(type == 3) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Exponentially Increasing\n\n";
    File << "#XS,EI\n";
    return std::make_unique<Exp_Increase>();
  }
  else if(type == 5) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Sharp Gaussian\n\n";
    File << "#XS,SG\n";
    return std::make_unique<Gauss_Sharp>();
  }
  else if(type == 6) {
    std::cout << "\n------------------------------------------------------";
    std::cout << "\n Broad Gaussian\n\n";
    File << "#XS,BG\n";
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
    pcg_seeds.resize(NTHREADS);
    for(int i = 0; i < NTHREADS; i++) {
        uint64_t seed = i+1;
        pcg_seeds[seed];
    }
}

int main() {
#ifdef _OPENMP
    omp_set_num_threads(NTHREADS);
#endif

    std::cout << "\n NParticles = " << NPART << ", NBins = " << NBIN << "\n\n";

    RNG_Seeds();

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

    File.open("Coll_Densities.txt");

    int xs_type = 1;
    //int trk_type = 1;
    
    std::unique_ptr<XS> crs = make_cross_section(xs_type);
    
    // Make initial source distribution
    std::vector<Particle> particle_bank;
    for(int i = 0; i < NPART; i++) {
        particle_bank.push_back(Particle(0.0,1.0,1.0));
        // Initials source, particles start at x = 0.0 in bin 0,
        // in forward direction, with wgt = 1.0
    }

    Zero_Values();
    File << "#TM,DT\n";
    Delta_Tracking(crs, particle_bank);
    Output();

    Zero_Values();
    File << "#TM,MDT\n";
    Meshed_Delta_Tracking(crs, particle_bank);
    Output();

    Zero_Values();
    File << "#TM,NWDT\n";
    Negative_Weight_Delta_Tracking(crs, particle_bank);
    Output();

    Zero_Values();
    File << "#TM,MNWDT\n";
    Meshed_Negative_Weight_Delta_Tracking(crs,particle_bank);
    Output();

    Zero_Values();
    File << "#TM,BT\n";
    Bomb_Transport(crs,0.8,particle_bank);
    Output();

    Zero_Values();
    File << "#TM,MBT\n";
    Meshed_Bomb_Transport(crs,0.8,particle_bank);
    Output();

    Zero_Values();
    File << "#TM,IBT\n";
    Improving_Meshed_Bomb_Transport(crs,0.8,particle_bank);
    Output();
    
    File.close();

    return 0;
}
