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
#include<iostream>

const double EPS = 1e-6;
const int NPART = 1e7;
const int NBIN = 2;
const double dx = 2.0/(double)NBIN;
const double p = 0.7;
const double p_mshd = 0.5;
const double q = 0.3;
const double q_mshd = 0.05;

// Constants for gaussians A*exp(-a*(x-z)^2)
const double A = 2.0/std::sqrt(2.0*M_PI);
const double a_s = (1.0/0.05)*(1.0/0.05);
const double a_b = 1.0;
const double z = 1.23;

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

        // Data
        double Pnc;
        double G;
        double Emax;
        double Em[NBIN];
        double Esmp[NBIN];
};

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
        if ((x > x_hi) or (x < x_low)) {
            x = False_Position(xs, T_hat, x_low, x_hi, EPS, counter);
            break;           
        }
    }
    return x;
}

void Direct_Sampling(XS* xs) {
    std::cout << "\n Direct Sampling (Newton's Method)\n";
    double escape = 0.0;
    double collide = 0.0;
    int cnts_sum = 0;

    #pragma omp parallel
    {
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double xi = rand(rng);
            int cnt = 0;
            if(xi < xs->Pnc) { // No collision
                #pragma omp atomic
                escape += 1.0;    
            } else { // Collision will occur 
                double T_hat = -std::log(1.0 - (xs->G)*xi);
                double x_low = 0.0;
                double x_hi = 2.0;
                double eps = 0.01;
                double x0 = False_Position(xs, T_hat, x_low, x_hi, eps, cnt);
                double x1 = Newton(xs, T_hat, x0, x_low, x_hi, cnt);

                if(x1 > 2.0) {
                    escape += 1;
                    #pragma omp critical
                    {
                    std::cout << " Problem with Newton\n";
                    }
                } else { 
                    #pragma omp atomic
                    collide += 1.0;
                }
            }
            #pragma omp atomic
            cnts_sum += cnt;
        }
    }

    double trans = escape / (double)NPART;
    double avg_cnt = (double)cnts_sum/(double)NPART;

    std::cout << " Collisions: " << (int)collide << ", Transmission: ";
    std::cout << trans << ", Average Counts: " << avg_cnt << "\n\n";
}

void Delta_Tracking(XS* xs) {
    std::cout << "\n Delta Tracking\n";

    int cnts_sum = 0;
    double escape = 0.0;
    double collide = 0.0;

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
                    virtual_collision = false;
                    #pragma omp atomic
                    cnts_sum += cnt;
                }
                else if(rand(rng) < (xs->Et(x)/xs->Emax)) {
                    // Collision is real
                    #pragma omp atomic
                    collide += 1.0;
                    #pragma omp atomic
                    cnts_sum += cnt;
                    virtual_collision = false;
                
                }
            }
        }
    }
    
    double trans = escape / (double)NPART;
    double avg_cnt = (double)cnts_sum/(double)NPART;

    std::cout << " Collisions: " << (int)collide << ", Transmission: ";
    std::cout << trans << ", Average Counts: " << avg_cnt << "\n\n";
}

void Meshed_Delta_Tracking(XS* xs) {
    std::cout << "\n Meshed Delta Tracking\n";

    int cnts_sum = 0;
    double escape = 0.0;
    double collide = 0.0;
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
    
    double trans = escape / (double)NPART;
    double avg_cnt = (double)cnts_sum/(double)NPART;
    double avg_bin_cnt = (double)bin_cnt_sum/(double)NPART;
    double avg_virtual_cnt = (double)virtual_cnt_sum/(double)NPART;

    std::cout << " Collisions: " << (int)collide << ", Transmission: ";
    std::cout << trans << ", Average Counts: " << avg_cnt << "\n";
    std::cout << " Avg Bin Cnts = " << avg_bin_cnt << ", Avg Virt. Cnts = ";
    std::cout << avg_virtual_cnt << "\n\n";
}

void Negative_Weight_Delta_Tracking(XS* xs) {
    std::cout << "\n Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    double escape = 0.0;
    double collide = 0.0;
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
                    alive = false;
                    #pragma omp atomic
                    cnts_sum += cnt;
                }
                else if(rand(rng) < q) {
                    // Collision is real
                    #pragma omp atomic
                    collide += w;
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
    
    double trans = escape / (double)NPART;
    double avg_cnt = (double)cnts_sum/(double)NPART;
    double avg_sign_chng = sign_change / (double)NPART;

    std::cout << " Collisions: " << (int)std::round(collide) << ", Transmission: ";
    std::cout << trans << ", Average Counts: " << avg_cnt;
    std::cout << "\n Avg Sign Chng = " << avg_sign_chng << "\n\n";
}

void Meshed_Negative_Weight_Delta_Tracking(XS* xs) {
    std::cout << "\n Meshed Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    double escape = 0.0;
    double collide = 0.0;
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
                        collide += w;
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
    
    double trans = escape / (double)NPART;
    double avg_cnt = (double)cnts_sum/(double)NPART;
    double avg_bin_cnt = (double)bin_cnt_sum/(double)NPART;
    double avg_sign_chng = sign_change / (double)NPART;

    std::cout << " Collisions: " << (int)collide << ", Transmission: ";
    std::cout << trans << ", Average Counts: " << avg_cnt << "\n";
    std::cout << " Avg Bin Cnts = " << avg_bin_cnt;
    std::cout << "\n Avg Sign Chng = " << avg_sign_chng << "\n\n";
}

int main() {
    std::cout << "\n NParticles = " << NPART << ", NBins = " << NBIN << "\n\n";
    
    for(int type = 1; type <= 6; type++) {
        XS* crs;
        
        // Determine cross section type for run
        if(type == 1) {
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

        Direct_Sampling(crs);

        Delta_Tracking(crs);

        Meshed_Delta_Tracking(crs);

        Negative_Weight_Delta_Tracking(crs);

        Meshed_Negative_Weight_Delta_Tracking(crs);
    }
    
    return 0;
}
