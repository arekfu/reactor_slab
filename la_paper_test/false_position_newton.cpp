/*
 *  Test of "Direct Sampling of Monte Carlo Flight Paths In Media
 *  With Continuously Varying Cross-Sections" (LA-UR-02-6530)
 *  by Forrest Brown and Bill Martin.
 *
 *  Implements solution of direct path lenght, delta tracking,
 *  and meshed delta tracking.
 *
 * */

#include"random/pcg_random.hpp"

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
const double Pnc_s = std::exp(-T_2_s);

const double T_2_b = A*std::sqrt(M_PI)*((std::erf(std::sqrt(a_b)*(2.0 - z)) 
                     - std::erf(std::sqrt(a_b)*(-z)))/(2.0*std::sqrt(a_b)));
const double Pnc_b = std::exp(-T_2_b);


// Base Cross Section Class
class XS {
    public:
        XS(double _Pnc, double _Emax) {
            Pnc = _Pnc;
            G = 1.0 - Pnc;
            Emax = _Emax;
        }

        // All virtual methods
        virtual double T(double x) {return -1.0;} // Calc. optical depth from 0 to x>0
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
        Constant():XS(std::exp(-2), 1.0) {
            for(int b = 0; b < NBIN; b++) {
                Em[b] = 1.0;
                Esmp[b] = p_mshd;
            }
        }

        double T(double x) {
            if((x >= 0.0) and (x <= 2.0)) {return x;}
            else if(x > 2.0) {return 2.0;}
            else {exit(1);}
        }

        double Et(double x) {return 1.0;}

}; // Constants

// Step XS (Et = 1.5 ∀ 0<x<0.278), Et=1.0 ∀ 0.278<x<2)
class Step : public XS {
    public:
        Step():XS(std::exp(-2.139)  , 1.5) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double T(double x) {
            if((x >= 0.0) and (x < 0.278)) {return 1.5*x;}
            else if((x >= 0.278) and x <= 2.0) {return 1.5*0.278 + (x - 0.278);}
            else if(x > 2.0) {return T(2.0);}
            else {
                #pragma omp critical
                {
                    std::cout << x<<"\n";
                }
                exit(1);}
        }

        double Et(double x) {
            if((x >= 0.0) and (x < 0.278)) {return 1.5;}
            else if ((x >= 0.278) and (x < 2.0)) {return 1.0;}
            else {exit(1);}
        }

}; // Step

class Lin_Decrease : public XS {
    public:
        Lin_Decrease():XS(std::exp(-2.0), 2.0) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }    
        }

        double T(double x) {
            if((x >= 0.0) and (x <= 2.0)) {
                return (2.0*x - ((x*x)/2.0));
            }
            else if(x > 2) {return T(2.0);}
            else {exit(1);}
        }

        double Et(double x) {
            return 2.0 - x;
        }
};

class Lin_Increase : public XS {
    public:
        Lin_Increase():XS(std::exp(-2.0), 2.0) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx + dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double T(double x) {
            if((x >= 0.0) and (x <= 2.0)) {
                return (x*x)/2.0 ;
            }
            else if(x > 2.0) {return T(2.0);}
            else {exit(1);}
        }

        double Et(double x) {
            return x;
        }
};

class Exp_Decrease : public XS {
    public:
        Exp_Decrease():XS(std::exp(-(1.0/3.0)*(1.0 - std::exp(-6.0))), 1.0) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double T(double x) {
            if((x >= 0.0) and (x <= 2.0)) {
                return ((1.0/3.0) - (std::exp(-3.0*x)/3.0));
            }
            else if(x > 2.0) {return T(2.0);}
            else {exit(1);}
        }

        double Et(double x) {
            return std::exp(-3.0*x);
        }
};

class Exp_Increase : public XS {
    public:
        Exp_Increase():XS(std::exp(-0.05*(std::exp(4.0)-1.0)), 0.1*std::exp(4.0)) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx + dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double T(double x) {
            if((x >= 0.0) and (x <= 2.0)) {
                return 0.05*(std::exp(2.0*x) - 1.0);
            }
            else if(x > 2.0) {return T(2.0);}
            else {exit(1);}
        }

        double Et(double x) {
            return 0.1*std::exp(2.0*x);
        }
};

class Gauss_Sharp : public XS {
    public:
        Gauss_Sharp():XS(Pnc_s, A) {
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

        double T(double x) {
            if((x >= 0.0) and (x <= 2.0)) {
                return A*std::sqrt(M_PI)*(std::erf(std::sqrt(a_s)*(x - z)) 
                       - std::erf(std::sqrt(a_s)*(-z)))/(2.0*std::sqrt(a_s));
            }
            else if(x > 2.0) {return T(2.0);}
            else {exit(1);}
        }

        double Et(double x) {
            return A*std::exp(-a_s*(x-z)*(x-z));
        }
};

class Gauss_Broad : public XS {
    public:
        Gauss_Broad():XS(Pnc_b, A) {
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

        double T(double x) {
            if((x >= 0.0) and (x <= 2.0)) {
                return A*std::sqrt(M_PI)*(std::erf(std::sqrt(a_b)*(x - z)) 
                       - std::erf(std::sqrt(a_b)*(-z)))/(2.0*std::sqrt(a_b));
            }
            else if(x > 2.0) {return T(2.0);}
            else {exit(1);}
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

double False_Position(std::unique_ptr<XS> const &xs, double T_hat, double& x0, double& x1, 
                      double eps, int& counter) {
    // x0 is starting location
    // x1 is point where particle would cross to new region
    
    double F1 = T_hat - xs->T(x1);
    double F0 = T_hat - xs->T(x0);
    counter++;
    counter++;
    double m = (F1 - F0)/(x1 - x0);
    double x = x0;
    double x_old = x1;
    if( m >= 0.0 ) {
        std::cout << " ERROR WITH INITIAL POINTS!!\n";
    } else {
        while(std::abs(x_old - x) > eps) {
            x_old = x;
            x = x0 - (F0/m);

            if(std::abs(x_old - x) < eps) {
                break;
            } else {
                double Fx = T_hat - xs->T(x);
                counter++;
                if(Fx > 0.0) {
                    x0 = x;
                    F0 = Fx;
                }
                else {
                    x1 = x;
                    F1 = Fx;
                }
                
                m = (F1 - F0)/(x1 - x0);
            }
        }
    }
    return x;
}

double Newton(std::unique_ptr<XS> const &xs, double T_hat, double x0, double x_low, 
              double x_hi, int& counter) {

    //if(T_hat == x0) {std::cout << " SCREAM\n";} 
    double x = x0;
    x0 += 0.05;
    double g = 0.0;
    double gp = 0.0;
    while(std::abs(x - x0) > EPS) {
        counter++;
        counter++;
        x0 = x;
        g = T_hat - xs->T(x0);
        gp = xs->Et(x0);
        x = x0 - (g/gp);
        if ((x > x_hi) or (x < x_low) or (x < 0.0)) {
            x = False_Position(xs, T_hat, x_low, x_hi, EPS, counter);
            break;           
        }
    }
    return x;
}

void Direct_Sampling(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Direct Sampling (False Position/Newton's Method)\n";
    double Pnc = xs->Pnc;
    #pragma omp parallel
    {
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double xi = rand(rng);
            int cnt = 0;
            if(xi < Pnc) { // No collision
                #pragma omp atomic
                escape += 1.0;    
                #pragma omp atomic
                escape_sqr += 1.0;
            } else { // Collision will occur
                xi = rand(rng); 
                double T_hat = -std::log(1.0 - (xs->G)*xi);
                double x_low = 0.0;
                double x_hi = 2.0;
                double eps = 0.01;
                double x0 = False_Position(xs, T_hat, x_low, x_hi, eps, cnt);
                double x1 = Newton(xs, T_hat, x0, x_low, x_hi, cnt);
                int coll_bin = std::floor(x1/Fdx);
                #pragma omp atomic
                coll_density[coll_bin][0] += 1.0;
                #pragma omp atomic
                coll_density[coll_bin][1] += 1.0;
                #pragma omp atomic
                collide += 1.0;
                #pragma omp atomic
                collide_sqr += 1.0;
                
            }
            #pragma omp atomic
            xs_evals += cnt;
        }
    }
}

void Delta_Tracking(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Delta Tracking\n";

    int cnts_sum = 0;

    #pragma omp parallel
    {    
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            bool virtual_collision = true;
            double x = 0.0;
            int cnt = 0;
            while(virtual_collision) {
                cnt++;
                double d = -std::log(rand(rng))/(xs->Emax);
                x += d;
                if(x >= 2.0) {
                    #pragma omp atomic
                    escape += 1.0;
                    #pragma omp atomic
                    escape_sqr += 1.0;
                    virtual_collision = false;
                    #pragma omp atomic
                    cnts_sum += cnt;
                }
                else if(rand(rng) < (xs->Et(x)/xs->Emax)) {
                    // Collision is real
                    #pragma omp atomic
                    collide += 1.0;
                    #pragma omp atomic
                    collide_sqr += 1.0;
                    #pragma omp atomic
                    cnts_sum += cnt;
                    int coll_bin = std::floor(x/Fdx);
                    #pragma omp atomic
                    coll_density[coll_bin][0] += 1.0;
                    #pragma omp atomic
                    coll_density[coll_bin][1] += 1.0;
                    virtual_collision = false;
                }
            }
        }
    }
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
            double Emax = xs->Em[0];
            int bin = 0;
            double d,d_bin;
            int cnt = 0;
            int virtual_cnt = 0;
            int bin_cnt = 0;
            while(virtual_collision) {
                d_bin = ((double)bin*dx + dx) - x;
                d = -std::log(rand(rng))/Emax;
                if(d_bin < d) {
                    bin_cnt++;
                    d = d_bin + 1e-6;
                    x += d;
                    bin = std::floor(x/dx);
                    if(x >= 2.0) {
                        virtual_collision = false;
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
                    x += d;
                    double xi = rand(rng);
                    double Pr = xs->Et(x)/Emax;
                    if(xi < Pr) {
                        virtual_collision = false;
                        #pragma omp atomic
                        collide += 1.0;
                        #pragma omp atomic
                        collide_sqr += 1.0;
                        int coll_bin = std::floor(x/Fdx);
                        #pragma omp atomic
                        coll_density[coll_bin][0] += 1.0;
                        #pragma omp atomic
                        coll_density[coll_bin][1] += 1.0;
                        #pragma omp atomic
                        cnts_sum += cnt;
                        #pragma omp atomic
                        virtual_cnt_sum += virtual_cnt;
                        #pragma omp atomic 
                        bin_cnt_sum += bin_cnt;
                    }
                    else {virtual_cnt++;}
                }
            }
        }
    }
    xs_evals += cnts_sum;
}

void Negative_Weight_Delta_Tracking(std::unique_ptr<XS> const &xs) {
    std::cout << "\n Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    double sign_change = 0.0;
    double Esamp = p*(xs->Emax);
    #pragma omp parallel
    {    
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double x = 0.0;
            double w = 1.0;
            bool alive = true;
            int cnt = 0;
            while(alive) {
                double d = -std::log(rand(rng))/Esamp;
                x += d;
                
                if(x >= 2.0) {
                    #pragma omp atomic
                    escape += w;
                    #pragma omp atomic
                    escape_sqr += w*w;
                    alive = false;
                    #pragma omp atomic
                    cnts_sum += cnt;
                }
                else if(rand(rng) < q) {
                    // Collision is real
                    #pragma omp atomic
                    collide += w*xs->Et(x)/(Esamp * q);
                    #pragma omp atomic
                    collide_sqr += (w*xs->Et(x)/(Esamp * q))*(w*xs->Et(x)/(Esamp * q));
                    int coll_bin = std::floor(x/Fdx);
                    #pragma omp atomic
                    coll_density[coll_bin][0] += w*xs->Et(x)/(Esamp * q);
                    #pragma omp atomic
                    coll_density[coll_bin][1] += (w*xs->Et(x)/(Esamp * q))*(w*xs->Et(x)/(Esamp * q));
                    #pragma omp atomic
                    cnts_sum += cnt;
                    alive = false;
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
            int bin = 0;
            double Esamp = xs->Esmp[bin];
            double d,d_bin;
            int cnt = 0;
            int bin_cnt = 0;
            double w = 1.0;
            while(alive) {
                d_bin = ((double)bin*dx + dx) - x;
                d = -std::log(rand(rng))/Esamp;
                if(d_bin < d) {
                    bin_cnt++;
                    d = d_bin + 1e-6;
                    x += d;
                    bin = std::floor(x/dx);
                    if(x >= 2.0) {
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
                    x += d;
                    if(rand(rng) < q_mshd) {
                        alive = false;
                        #pragma omp atomic
                        collide += w*(xs->Et(x)/(Esamp*q_mshd));
                        #pragma omp atomic
                        collide_sqr += (w*(xs->Et(x)/(Esamp*q_mshd)))*(w*(xs->Et(x)/(Esamp*q_mshd)));
                        int coll_bin = std::floor(x/Fdx);
                        #pragma omp atomic
                        coll_density[coll_bin][0] += w*(xs->Et(x)/(Esamp*q_mshd));
                        #pragma omp atomic
                        coll_density[coll_bin][1] += (w*(xs->Et(x)/(Esamp*q_mshd)))*(w*(xs->Et(x)/(Esamp*q_mshd)));
                        #pragma omp atomic
                        cnts_sum += cnt;
                        #pragma omp atomic 
                        bin_cnt_sum += bin_cnt;
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
    double Esmp = P*xs->Emax;
    int cnts_sum = 0;
    double sign_change = 0.0;
    #pragma omp parallel
    {
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double x = 0.0;
            bool alive = true;
            int cnt = 0;
            double w = 1.0;
            bool real_collision = false;

            while(alive) {
                double d = -std::log(rand(rng))/Esmp;
                x += d;
                
                // Fist check for leak
                if(x > 2.0) {
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
                        alive = false;
                        #pragma omp atomic
                        collide += w;
                        #pragma omp atomic
                        collide_sqr += w*w;
                        int coll_bin = std::floor(x/Fdx);
                        #pragma omp atomic
                        coll_density[coll_bin][0] += w;
                        #pragma omp atomic
                        coll_density[coll_bin][1] += w*w;
                        #pragma omp atomic
                        cnts_sum += cnt;
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
    double Esmp = xs->Et(0.0); // Initial sampling XS is XS at birth place
    if(Esmp < 0.1) Esmp = 0.1;
    int cnts_sum = 0;
    double sign_change = 0.0;
    #pragma omp parallel
    {
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double x = 0.0;
            bool alive = true;
            int cnt = 0;
            double w = 1.0;
            bool real_collision = false;

            while(alive) {
                double d = -std::log(rand(rng))/Esmp;
                x += d;
                
                // Fist check for leak
                if(x > 2.0) {
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
                        alive = false;
                        #pragma omp atomic
                        collide += w;
                        #pragma omp atomic
                        collide_sqr += w*w;
                        int coll_bin = std::floor(x/Fdx);
                        #pragma omp atomic
                        coll_density[coll_bin][0] += w;
                        #pragma omp atomic
                        coll_density[coll_bin][1] += w*w;
                        #pragma omp atomic
                        cnts_sum += cnt;
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
    std::cout << "\n Meshed Bomb Paper Transport, p = " << P  << "\n";
    int cnts_sum = 0;
    double sign_change = 0.0;
    #pragma omp parallel
    {
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double x = 0.0;
            int bin = 0;
            double Esmp = P*xs->Em[bin];
            bool alive = true;
            int cnt = 0;
            double w = 1.0;
            bool real_collision = false;

            while(alive) {
                double d = -std::log(rand(rng))/Esmp;
                double d_bin = ((double)bin*dx + dx) - x;
                
                // Fist check for bin change
                if(d_bin < d) {
                    x += d_bin; // Move to new bin
                    bin = std::floor(x/dx);
                    if(x > 2.0) {
                        alive = false;
                        #pragma omp atomic
                        escape += w;
                        #pragma omp atomic
                        escape_sqr += w*w;
                        #pragma omp atomic
                        cnts_sum += cnt;
                    } else {
                        Esmp = P*xs->Em[bin];
                    }
                    
                } else if (x + d > 2.0) {
                    alive = false;
                    #pragma omp atomic
                    escape += w;
                    #pragma omp atomic
                    escape_sqr += w*w;
                    #pragma omp atomic
                    cnts_sum += cnt;
                } else {
                    x += d;
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
                        alive = false;
                        #pragma omp atomic
                        collide += w;
                        #pragma omp atomic
                        collide_sqr += w*w;
                        int coll_bin = std::floor(x/Fdx);
                        #pragma omp atomic
                        coll_density[coll_bin][0] += w;
                        #pragma omp atomic
                        coll_density[coll_bin][1] += w*w;
                        #pragma omp atomic
                        cnts_sum += cnt;
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
      File << "#TM,DS\n";
      Direct_Sampling(crs);
      Output();

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
