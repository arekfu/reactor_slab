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

const double EPS = 1e-6;
const int NPART = 1e7;
const int NBIN = 5;
const double dx = 2.0/(double)NBIN;
const double p = 0.7;
const double p_mshd = 0.7;
const double q = 0.3;
const double q_mshd = 0.3;

// Constants for gaussians A*exp(-a*(x-z)^2)
const double A = 2.0/std::sqrt(2.0*M_PI);
const double a_s = (1.0/0.05)*(1.0/0.05);
const double a_b = 1.0;
const double z = 1.23;

// Tallies
double collide;
double collide_sqr;
double escape;
double escape_sqr;
double xs_evals; // # of xs look ups or xs integrations

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
        virtual double T(double x) {return -1.0;}
        virtual double Et(double x) {return -1.0;}
        virtual double dEt(double x) {return -1000000.0;}

        double P_nc() {return Pnc;}

        // Data
        double Pnc;
        double G;
        double Emax;
        double Em[NBIN];
        double Esmp[NBIN];
};

// Constant XS
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

        double dEt(double x) {return 0.0;}

}; // Constants

// Step XS
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

        double dEt(double x) {return 0.0;}
}; // Step

// Derived classes for each type of cross section
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

        double dEt(double x) {
            return -1.0;
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

        double dEt(double x) {
            return 1.0;
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

        double dEt(double x) {
            return -3.0*std::exp(-3.0*x);
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

        double dEt(double x) {
            return 0.2*std::exp(2.0*x);
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

        double dEt(double x) {
            double ar = (x - 1.0)/0.05;
            return (2.0/std::sqrt(2.0*M_PI))*(-2*ar*(1.0/0.05))*std::exp(-ar*ar);
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

        double dEt(double x) {
            double ar = (x - 1.0);
            return (2.0/std::sqrt(2.0*M_PI))*(-2*ar)*std::exp(-ar*ar);
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

double False_Position(XS* xs, double T_hat, double& x0, double& x1, 
                      double eps, int& counter) {
    // x0 is starting location
    // x1 is point where particle would cross to new region
    
    // Lambda function for T_hat - T(x)
    auto F = [&](double y) {
        return (T_hat - xs->T(y));
    };
    
    double F1 = F(x1);
    double F0 = F(x0);
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

            if((x > 2.0) or (x < 0.0)) {
                #pragma omp critical
                {
                std::cout << " PROBS\n";
                std::cout << " T_hat = " << T_hat << "\n";
                std::cout << " x0 = " << x0 << "\n" << " x1 = " << x1;
                std::cout << "\n x = " << x << " m = " << m << "\n";
                std::cout << " F0 = " << F0 << " F1 = " << F1 << "\n";
                std::cin >> x;
                }
            }

            if(std::abs(x_old - x) > eps) {
                break;
            } else {
                counter++;
                if(F(x) > 0.0) {
                    x0 = x;
                    F0 = F(x0);
                }
                else {
                    x1 = x;
                    F1 = F(x1);
                }
                
                m = (F1 - F0)/(x1 - x0);
            }
        }
    }
    return x;
}

double Newton(XS* xs, double T_hat, double x0, double x_low, 
              double x_hi, int& counter) {
    // Lambda function for T_hat - T(x)
    auto F = [&](double y) {
        return (T_hat - xs->T(y));
    };

    //if(T_hat == x0) {std::cout << " SCREAM\n";} 
    double x = x0;
    x0 += 0.05;
    double g = 0.0;
    double gp = 0.0;
    while(std::abs(x - x0) > EPS) {
        counter++;
        counter++;
        x0 = x;
        g = F(x0);
        gp = xs->Et(x0);
        x = x0 - (g/gp);
        if ((x > x_hi) or (x < x_low) or (x < 0.0)) {
            x = False_Position(xs, T_hat, x_low, x_hi, EPS, counter);
            break;           
        }
    }
    return x;
}

void Direct_Sampling(XS* xs) {
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
                double T_hat = -std::log(1.0 - (xs->G)*xi);
                double x_low = 0.0;
                double x_hi = 2.0;
                double eps = 0.01;
                double x0 = False_Position(xs, T_hat, x_low, x_hi, eps, cnt);
                double x1 = Newton(xs, T_hat, x0, x_low, x_hi, cnt);

                if(x1 > 2.0) {
                    escape += 1.0;
                    #pragma omp critical
                    {
                    std::cout << " Problem with Newton\n";
                    }
                } else { 
                    #pragma omp atomic
                    collide += 1.0;
                    #pragma omp atomic
                    collide_sqr += 1.0;
                }
            }
            #pragma omp atomic
            xs_evals += cnt;
        }
    }
}

void Delta_Tracking(XS* xs) {
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
                    virtual_collision = false;
                }
            }
        }
    }
    xs_evals += cnts_sum;
}

void Meshed_Delta_Tracking(XS* xs) {
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

void Negative_Weight_Delta_Tracking(XS* xs) {
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
}

void Meshed_Negative_Weight_Delta_Tracking(XS* xs) {
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

    double FOM_col = 1.0 / (avg_xs_evals * collide_std * collide_std);
    double FOM_esc = 1.0 / (avg_xs_evals * escape_std * escape_std);
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << " Colsn. Rate: " << collide_avg << " +/- " << collide_std;
    std::cout << ", Trans. Rate: " << escape_avg << " +/- " << escape_std << "\n";
    std::cout << " Average XS/Integration Evaluations: " << avg_xs_evals << "\n";
    std::cout << std::scientific;
    std::cout << " FOM_coll = " << FOM_col << ", FOM_escp = " << FOM_esc << "\n\n";
}

int main() {
    std::cout << "\n NParticles = " << NPART << ", NBins = " << NBIN << "\n\n";
    
    for(int type = 1; type <= 6; type++) {
        XS* crs;
        
        // Determine cross section type for run
        if(type == -1) {
            std::cout << "\n------------------------------------------------------";
            std::cout << "\n Step\n\n";
            Step xs = Step();
            crs = &xs;
        }
        else if(type == 0) {
            std::cout << "\n------------------------------------------------------";
            std::cout << "\n Constant\n\n";
            Constant xs = Constant();
            crs = &xs;
        }
        else if(type == 1) {
            std::cout << "\n------------------------------------------------------";
            std::cout << "\n Linearly Increasing\n\n";
            Lin_Increase xs = Lin_Increase();
            crs = &xs;
        }
        else if(type == 2) {
            std::cout << "\n------------------------------------------------------";
            std::cout << "\n Linearly Decreasing\n\n";
            Lin_Decrease xs = Lin_Decrease();
            crs = &xs;
        }
        else if(type == 4) {
            std::cout << "\n------------------------------------------------------";
            std::cout << "\n Exponentially Decreasing\n\n";
            Exp_Decrease xs = Exp_Decrease();
            crs = &xs;
        }
        else if(type == 3) {
            std::cout << "\n------------------------------------------------------";
            std::cout << "\n Exponentially Increasing\n\n";
            Exp_Increase xs = Exp_Increase();
            crs = &xs;
        }
        else if(type == 5) {
            std::cout << "\n------------------------------------------------------";
            std::cout << "\n Sharp Gaussian\n\n";
            Gauss_Sharp xs = Gauss_Sharp();
            crs = &xs;
        }
        else if(type == 6) {
            std::cout << "\n------------------------------------------------------";
            std::cout << "\n Broad Gaussian\n\n";
            Gauss_Broad xs = Gauss_Broad();
            crs = &xs;
        }
        else {exit(1);}
        
        collide = 0.0;
        collide_sqr = 0.0;
        escape = 0.0;
        escape_sqr = 0.0;
        xs_evals = 0.0;

        Direct_Sampling(crs);
        Output();
                
        collide = 0.0;
        collide_sqr = 0.0;
        escape = 0.0;
        escape_sqr = 0.0;
        xs_evals = 0.0;

        Delta_Tracking(crs);
        Output();
        
        collide = 0.0;
        collide_sqr = 0.0;
        escape = 0.0;
        escape_sqr = 0.0;
        xs_evals = 0.0;

        Meshed_Delta_Tracking(crs);
        Output();
        
        collide = 0.0;
        collide_sqr = 0.0;
        escape = 0.0;
        escape_sqr = 0.0;
        xs_evals = 0.0;

        Negative_Weight_Delta_Tracking(crs);
        Output();
        
        collide = 0.0;
        collide_sqr = 0.0;
        escape = 0.0;
        escape_sqr = 0.0;
        xs_evals = 0.0;

        Meshed_Negative_Weight_Delta_Tracking(crs);
        Output();

    }
    
    return 0;
}
