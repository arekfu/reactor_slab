/*
 *  Test of "Direct Sampling of Monte Carlo Flight Paths In Media
 *  With Continuously Varying Cross-Sections" (LA-UR-02-6530)
 *  by Forrest Brown and Bill Martin.
 *
 *  Implements solution of direct path lenght, delta tracking,
 *  and meshed delta tracking.
 *
 * */

#include"pcg_random.hpp"

#include<cmath>
#include<iomanip>
#include<iostream>
#include<memory>
#include<vector>
#include<fstream>

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

// Constants for gaussians A*exp(-a*(x-z)^2)
const double A = 2.0/std::sqrt(2.0*M_PI);
const double a_s = (1.0/0.05)*(1.0/0.05);
const double a_b = 1.0;
const double z = 1.23;

// Outputfile
std::ofstream File;

// Tallies
double collide; // Sum of weights that collide
double collide_sqr; // Sum of weights squared that collide
double escape; // Sum of weights that leak
double escape_sqr; // Sum of weights squared that lead
double xs_evals; // # of xs look ups or xs integrations
double wgt_chngs; // # of times particle wgt sign flipped
std::vector<std::vector<double>> coll_density; //[0] #sum coll in box,[1] sum coll sqr, [2] coll STD, [3] Coll FOM in box

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
        virtual double Et(double x) {return -1.0;}

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

        double Et(double x) {return 1.0;}

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

// Random number generation function
double rand(pcg64_unique& rng) {
    uint32_t range = 1000000000;
    uint32_t mask = ~uint32_t(0);
    --range;
    mask >>= __builtin_clz(range|1);
    uint32_t x;
    do {
        x = rng() & mask;
    } while(x > range);
    return double(x)/1000000000.0;
}

void Delta_Tracking(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Delta Tracking\n";

    int cnts_sum = 0;

    #pragma omp parallel
    {    
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            bool virtual_collision;
            double x = 0.0;
            double u = 1.0;
            int cnt = 0;
            bool alive = true;
            while(alive) {
                virtual_collision = true;
                while(virtual_collision) {
                    double d = -std::log(rand(rng))/(xs->Emax);
                    cnt++;
                    x += u*d;
                    if((x >= 2.0) or (x <= 0.0)) {
                        #pragma omp atomic
                        escape += 1.0;
                        #pragma omp atomic
                        escape_sqr += 1.0;
                        virtual_collision = false;
                        alive = false;
                        #pragma omp atomic
                        cnts_sum += cnt;
                    } else if(rand(rng) < (xs->Et(x)/xs->Emax)) {
                        virtual_collision = false;
                    }
                }
                
                if(alive) {
                    // Real collision
                    #pragma omp atomic
                    collide += 1.0;
                    #pragma omp atomic
                    collide_sqr += 1.0;

                    int coll_bin = std::floor(x/Fdx);
                    #pragma omp atomic
                    coll_density[coll_bin][0] += 1.0;
                    #pragma omp atomic
                    coll_density[coll_bin][1] += 1.0;

                    // Get collision type
                    if(rand(rng) < P_abs) {
                        alive = false;
                        #pragma omp atomic
                        cnts_sum += cnt;
                    } else {
                        if(rand(rng) > P_straight_ahead) u *= -1;
                    }
                } // If alive for real collision
            } // While alive
        } // For all particles
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
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            bool virtual_collision = true;
            double x = 0.0;
            double u = 1.0;
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
                    d = -std::log(rand(rng))/Emax;
                    if(d_bin < d) {
                        bin_cnt++;
                        d = d_bin + 1e-6;
                        x += u*d;
                        bin = std::floor(x/dx);
                        if((x >= 2.0) or (x  <= 0.0)) {
                            virtual_collision = false;
                            alive = false;
                            #pragma omp atomic
                            escape += 1.0;
                            #pragma omp atomic
                            escape_sqr += 1.0;
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
                        double xi = rand(rng);
                        double Pr = xs->Et(x)/Emax;
                        if(xi < Pr) {
                            virtual_collision = false;
                        } else {virtual_cnt++;}
                    }
                }// While virtual

                if(alive) {
                    #pragma omp atomic
                    collide += 1.0;
                    #pragma omp atomic
                    collide_sqr += 1.0;
                    int coll_bin = std::floor(x/Fdx);
                    #pragma omp atomic
                    coll_density[coll_bin][0] += 1.0;
                    #pragma omp atomic
                    coll_density[coll_bin][1] += 1.0;

                    // Dertmine collision type
                    if(rand(rng) < P_abs) {
                        alive = false;
                        #pragma omp atomic
                        cnts_sum += cnt;
                        #pragma omp atomic
                        virtual_cnt_sum += virtual_cnt;
                        #pragma omp atomic 
                        bin_cnt_sum += bin_cnt;
                    } else { // Scatter
                        if(rand(rng) > P_straight_ahead) u *= -1;
                    }
                }
            } // While alive
        }// For all particles
    }// Parallel
    xs_evals += cnts_sum;
}

void Negative_Weight_Delta_Tracking(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    double sign_change = 0.0;
    #pragma omp parallel
    {    
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double Esamp = p*(xs->Emax);
            double x = 0.0;
            double u = 1.0;
            double w = 1.0;
            bool alive = true;
            int cnt = 0;
            while(alive) {
                double d = -std::log(rand(rng))/Esamp;
                x += u*d;
                
                if((x >= 2.0) or (x <= 0.0)) {
                    #pragma omp atomic
                    escape += w;
                    #pragma omp atomic
                    escape_sqr += w*w;
                    alive = false;
                    #pragma omp atomic
                    cnts_sum += cnt;
                }
                else if(rand(rng) < q) {
                    // Collision is real, addjust weight
                    w *= xs->Et(x)/(Esamp * q);
                    #pragma omp atomic
                    collide += w;
                    #pragma omp atomic
                    collide_sqr += w*w;
                    int coll_bin = std::floor(x/Fdx);
                    #pragma omp atomic
                    coll_density[coll_bin][0] += w;
                    #pragma omp atomic
                    coll_density[coll_bin][1] += w*w;

                    // Determine collision type
                    if(rand(rng) < P_abs) {
                        #pragma omp atomic
                        cnts_sum += cnt;
                        alive = false;
                    } else {
                        if(rand(rng) > P_straight_ahead) {u *= -1;}
                    }
                    
                }
                else { // Collision is virtual
                    double dw = ((1.0 - (xs->Et(x)/Esamp))/(1.0 - q));
                    if(dw < 0.0) {
                        #pragma omp atomic
                        sign_change += 1.0;
                    }
                    w = w*dw;
                    cnt++;
                }
            }
        }
    }
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Meshed_Negative_Weight_Delta_Tracking(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Meshed Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    int bin_cnt_sum = 0;
    double sign_change = 0.0;
    #pragma omp parallel
    {    
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            bool alive = true;
            double x = 0.0;
            double u = 1.0;
            int bin = 0;
            double Esamp = xs->Esmp[bin];
            double d,d_bin;
            int cnt = 0;
            int bin_cnt = 0;
            double w = 1.0;
            while(alive) {
                if(u == -1.0) d_bin = x - static_cast<double>(bin)*dx;
                else d_bin = (static_cast<double>(bin)*dx + dx) - x;
                d = -std::log(rand(rng))/Esamp;
                if(d_bin < d) {
                    bin_cnt++;
                    d = d_bin + 1e-6;
                    x += u*d;
                    bin = std::floor(x/dx);
                    if((x >= 2.0) or (x <= 0.0)) {
                        alive = false;
                        #pragma omp atomic
                        escape += w;
                        #pragma omp atomic
                        escape_sqr += w*w;
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
                    if(rand(rng) < q_mshd) {
                        // update weight
                        w *= (xs->Et(x)/(Esamp*q_mshd));
                        #pragma omp atomic
                        collide += w;
                        #pragma omp atomic
                        collide_sqr += w*w;
                        int coll_bin = std::floor(x/Fdx);
                        #pragma omp atomic
                        coll_density[coll_bin][0] += w;
                        #pragma omp atomic
                        coll_density[coll_bin][1] += w*w;

                        if(rand(rng) < P_abs) {
                            alive = false;
                            #pragma omp atomic
                            cnts_sum += cnt;
                            #pragma omp atomic 
                            bin_cnt_sum += bin_cnt;
                        } else {
                            if(rand(rng) > P_straight_ahead) u *= -1;
                        }
                        
                    }
                    else {
                        cnt++;
                        double dw=((1.0-(xs->Et(x)/Esamp))/(1.0-q_mshd));
                        if(dw < 0.0) {
                            #pragma omp atomic
                            sign_change += 1.0;
                        }
                        w = w*dw;
                    }
                }
            }
        }
    }
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Bomb_Transport(std::unique_ptr<XS> const &xs, double P) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Bomb Paper Transport, p = " << P << "\n";
    int cnts_sum = 0;
    double sign_change = 0.0;
    #pragma omp parallel
    {
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double Esmp = P*xs->Emax;
            double x = 0.0;
            double u = 1.0;
            bool alive = true;
            int cnt = 0;
            double w = 1.0;
            bool real_collision = false;

            while(alive) {
                double d = -std::log(rand(rng))/Esmp;
                x += u*d;
                real_collision = false;
                
                // Fist check for leak
                if((x >= 2.0) or (x <= 0.0)) {
                    alive = false;
                    #pragma omp atomic
                    escape += w;
                    #pragma omp atomic
                    escape_sqr += w*w;
                    #pragma omp atomic
                    cnts_sum += cnt;
                } else {
                    double E_tot = xs->Et(x);
                    cnt += 1;
                    if(E_tot > Esmp) { // First negative branch
                        double D = E_tot / (2*E_tot - Esmp);
                        double F = E_tot / (D*Esmp);
                        w *= F;
                        if(rand(rng) < D) {real_collision = true;}
                        else {
                            w *= -1.0;
                            #pragma omp atomic
                            sign_change += 1.0;
                        }

                    } else { // Delta tracking branch
                        double P_real = E_tot/ Esmp;
                        if(rand(rng) < P_real) {real_collision = true;}
                    }

                    if(real_collision) {
                        // Record collision
                        #pragma omp atomic
                        collide += w;
                        #pragma omp atomic
                        collide_sqr += w*w;
                        int coll_bin = std::floor(x/Fdx);
                        #pragma omp atomic
                        coll_density[coll_bin][0] += w;
                        #pragma omp atomic
                        coll_density[coll_bin][1] += w*w;

                        // Determine collision type
                        if(rand(rng) < P_abs) {
                            alive = false;
                            #pragma omp atomic
                            cnts_sum += cnt;
                        } else {
                            if(rand(rng) > P_straight_ahead) u *= -1;
                        }
                    }
                }
            }// While alive
        }// For all particles
    }// Parallel
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Meshed_Bomb_Transport(std::unique_ptr<XS> const &xs, double P) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Meshed Bomb Paper Transport, p = " << P << "\n";
    int cnts_sum = 0;
    double sign_change = 0.0;
    #pragma omp parallel
    {
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double Esmp = P*xs->Em[0];
            double x = 0.0;
            double u = 1.0;
            double d_bin;
            int bin = 0;
            bool alive = true;
            int cnt = 0;
            double w = 1.0;
            bool real_collision = false;

            while(alive) {
                double d = -std::log(rand(rng))/Esmp;
                if(u == -1.0) d_bin = x - static_cast<double>(bin)*dx;
                else d_bin = (static_cast<double>(bin)*dx + dx) - x;
                real_collision = false;
                
                if(d_bin < d) {
                    x += u*(d_bin+1e-6);
                    if((x >= 2.0) or (x <= 0.0)) {
                        alive = false;
                        #pragma omp atomic
                        escape += w;
                        #pragma omp atomic
                        escape_sqr += w*w;
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
                        #pragma omp atomic
                        escape += w;
                        #pragma omp atomic
                        escape_sqr += w*w;
                        #pragma omp atomic
                        cnts_sum += cnt;
                    } else {
                        double E_tot = xs->Et(x);
                        cnt += 1;
                        if(E_tot > Esmp) { // First negative branch
                            double D = E_tot / (2*E_tot - Esmp);
                            double F = E_tot / (D*Esmp);
                            w *= F;
                            if(rand(rng) < D) {real_collision = true;}
                            else {
                                w *= -1.0;
                                #pragma omp atomic
                                sign_change += 1.0;
                            }

                        } else { // Delta tracking branch
                            double P_real = E_tot/ Esmp;
                            if(rand(rng) < P_real) {real_collision = true;}
                        }

                        if(real_collision) {
                            // Record collision
                            #pragma omp atomic
                            collide += w;
                            #pragma omp atomic
                            collide_sqr += w*w;
                            int coll_bin = std::floor(x/Fdx);
                            #pragma omp atomic
                            coll_density[coll_bin][0] += w;
                            #pragma omp atomic
                            coll_density[coll_bin][1] += w*w;

                            // Determine collision type
                            if(rand(rng) < P_abs) {
                                alive = false;
                                #pragma omp atomic
                                cnts_sum += cnt;
                            } else {
                                if(rand(rng) > P_straight_ahead) u *= -1;
                            }
                        }
                    }
                }
            }// While alive
        }// For all particles
    }// Parallel
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Previous_XS_Bomb_Transport(std::unique_ptr<XS> const &xs) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\n Previous XS Bomb Paper Transport\n";
    int cnts_sum = 0;
    double sign_change = 0.0;
    #pragma omp parallel
    {
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double Esmp = xs->Et(0.0); // Initial sampling XS is XS at birth place
            if(Esmp < 0.1) Esmp = 0.1;
            double x = 0.0;
            double u = 1.0;
            bool alive = true;
            int cnt = 0;
            double w = 1.0;
            bool real_collision = false;

            while(alive) {
                real_collision = false;
                double d = -std::log(rand(rng))/Esmp;
                x += u*d;
                
                // Fist check for leak
                if((x >= 2.0) or (x <= 0.0)) {
                    alive = false;
                    #pragma omp atomic
                    escape += w;
                    #pragma omp atomic
                    escape_sqr += w*w;
                    #pragma omp atomic
                    cnts_sum += cnt;
                } else {
                    double E_tot = xs->Et(x);
                    cnt += 1;
                    if(E_tot > Esmp) { // First negative branch
                        double D = E_tot / (2*E_tot - Esmp);
                        double F = E_tot / (D*Esmp);
                        w *= F;
                        if(rand(rng) < D) {real_collision = true;}
                        else {
                            w *= -1.0;
                            #pragma omp atomic
                            sign_change += 1.0;
                        }

                    } else { // Delta tracking branch
                        double P_real = E_tot/ Esmp;
                        if(rand(rng) < P_real) {real_collision = true;}
                    }

                    Esmp = E_tot; // Update sampling XS to current XS at location
                    if(Esmp < 0.1) Esmp = 0.1;

                    if(real_collision) {
                        #pragma omp atomic
                        collide += w;
                        #pragma omp atomic
                        collide_sqr += w*w;
                        int coll_bin = std::floor(x/Fdx);
                        #pragma omp atomic
                        coll_density[coll_bin][0] += w;
                        #pragma omp atomic
                        coll_density[coll_bin][1] += w*w;

                        // Determine collision type
                        if(rand(rng) < P_abs) {
                            alive = false;
                            #pragma omp atomic
                            cnts_sum += cnt;
                        } else {
                            if(rand(rng) > P_straight_ahead) u *= -1;
                        }
                    }
                }
            }// While alive
        }// For all particles
    }// Parallel
    xs_evals += cnts_sum;
    wgt_chngs += sign_change;
}

void Output() {
    double collide_avg = collide / (double)NPART;
    double collide_sqr_avg = collide_sqr / (double)NPART;
    double collide_std = std::sqrt(std::abs(collide_avg*collide_avg - 
                          collide_sqr_avg)/((double)NPART - 1.0));
    double escape_avg = escape / (double)NPART;
    double escape_sqr_avg = escape_sqr / (double)NPART;
    double escape_std = std::sqrt(std::abs(escape_avg*escape_avg - 
                           escape_sqr_avg)/((double)NPART - 1.0));

    double avg_xs_evals = xs_evals / (double)NPART;
    double avg_sgn_chngs = wgt_chngs / (double)NPART;

    //  Calculations and output for collision density profile
    for(int i = 0; i < NFOMBINS; i++) {
        // Get avg for bin
        double coll_avg = coll_density[i][0] / static_cast<double>(NPART);
        double coll_sqr_avg = coll_density[i][1] / static_cast<double>(NPART);
        double coll_sig = std::sqrt(std::abs(coll_avg*coll_avg - coll_sqr_avg)
                /(static_cast<double>(NPART) - 1.0));
        coll_density[i][2] = coll_sig;
        double rel_error = coll_sig / coll_avg;
        coll_density[i][3] = 1.0 / (avg_xs_evals * rel_error * rel_error);

        // Output avg coll desnity in bin
        if(i == 0) {File << coll_avg;}
        else {File << "," << coll_avg;}
    }
    File << "\n";
    // Output std of coll density
    for(int i = 0; i < NFOMBINS; i++) {
        if(i == 0) {File << coll_density[i][2];}
        else {File << "," << coll_density[i][2];}
    }
    File << "\n";
    // Output FOM of coll density
    for(int i = 0; i < NFOMBINS; i++) {
        if(i == 0) {File << coll_density[i][3];}
        else {File << "," << coll_density[i][3];}
    }
    File << "\n\n";

    double coll_rel_error = collide_std / collide_avg;
    double escape_rel_error = escape_std / escape_avg;
    double FOM_col = 1.0 / (avg_xs_evals * coll_rel_error * coll_rel_error);
    double FOM_esc = 1.0 / (avg_xs_evals*escape_rel_error*escape_rel_error);
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << " Colsn. Rate: " << collide_avg << " +/- " << collide_std;
    std::cout << ", Trans. Rate: " << escape_avg << " +/- " << escape_std << "\n";
    std::cout << " Avg XS/Integration Evals: " << avg_xs_evals;
    std::cout << ", Avg Sign Changes: " << avg_sgn_chngs << "\n";
    std::cout << std::scientific;
    std::cout << " FOM_coll = " << FOM_col << ", FOM_escp = " << FOM_esc << "\n\n";
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

void Zero_Bins() {
    for(int i = 0; i < NFOMBINS; i++) {
        coll_density[i][0] = 0.0;
        coll_density[i][1] = 0.0;
        coll_density[i][2] = 0.0;
        coll_density[i][3] = 0.0;
    }
}

int main() {
    std::cout << "\n NParticles = " << NPART << ", NBins = " << NBIN << "\n\n";

     // Create and zero coll_density array
    for(int i = 0; i < NFOMBINS; i++) {
        std::vector<double> bin;
        bin.push_back(0.0);
        bin.push_back(0.0);
        bin.push_back(0.0);
        bin.push_back(0.0);
        coll_density.push_back(bin);
    }

    File.open("Coll_Densities.txt");
    
    for(int type = 1; type <= 6; type++) {
      std::unique_ptr<XS> crs = make_cross_section(type);

      collide = 0.0;
      collide_sqr = 0.0;
      escape = 0.0;
      escape_sqr = 0.0;
      xs_evals = 0.0;
      wgt_chngs = 0.0;
      Zero_Bins();
      File << "#TM,DT\n";
      Delta_Tracking(crs);
      Output();

      collide = 0.0;
      collide_sqr = 0.0;
      escape = 0.0;
      escape_sqr = 0.0;
      xs_evals = 0.0;
      wgt_chngs = 0.0;
      Zero_Bins();
      File << "#TM,MDT\n";
      Meshed_Delta_Tracking(crs);
      Output();

      collide = 0.0;
      collide_sqr = 0.0;
      escape = 0.0;
      escape_sqr = 0.0;
      xs_evals = 0.0;
      wgt_chngs = 0.0;
      Zero_Bins();
      File << "#TM,NWDT\n";
      Negative_Weight_Delta_Tracking(crs);
      Output();

      collide = 0.0;
      collide_sqr = 0.0;
      escape = 0.0;
      escape_sqr = 0.0;
      xs_evals = 0.0;
      wgt_chngs = 0.0;
      Zero_Bins();
      File << "#TM,MNWDT\n";
      Meshed_Negative_Weight_Delta_Tracking(crs);
      Output();

      collide = 0.0;
      collide_sqr = 0.0;
      escape = 0.0;
      escape_sqr = 0.0;
      xs_evals = 0.0;
      wgt_chngs = 0.0;
      Zero_Bins();
      File << "#TM,BT\n";
      Bomb_Transport(crs,0.95);
      Output();

      collide = 0.0;
      collide_sqr = 0.0;
      escape = 0.0;
      escape_sqr = 0.0;
      xs_evals = 0.0;
      wgt_chngs = 0.0;
      Zero_Bins();
      File << "#TM,MBT\n";
      Meshed_Bomb_Transport(crs,0.95);
      Output();
      
      collide = 0.0;
      collide_sqr = 0.0;
      escape = 0.0;
      escape_sqr = 0.0;
      xs_evals = 0.0;
      wgt_chngs = 0.0;
      Zero_Bins();
      File << "#TM,PBT\n";
      Previous_XS_Bomb_Transport(crs);
      Output();
    }
    File.close();
    return 0;
}
