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
const double Pscatt = 0.5;
const double Pabs = 0.5;

// Constants for gaussians A*exp(-a*(x-z)^2)
const double A = 2.0/std::sqrt(2.0*M_PI);
const double a_s = (1.0/0.05)*(1.0/0.05);
const double a_b = 1.0;
const double z = 1.23;

// Base Cross Section Class
class XS {
    public:
        XS(double _Emax) {
            Emax = _Emax;
        }

        // All virtual methods
        virtual double T(double x0, double x) {return -1.0;}
        virtual double Et(double x) {return -1.0;}

        // Data
        double Emax;
        double Em[NBIN];
        double Esmp[NBIN];
};

// Derived classes for each type of cross section
class Lin_Decrease : public XS {
    public:
        Lin_Decrease():XS(2.0) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }    
        }

        double T(double x0, double x) {
            double higher = 2*x - (x*x/2.0);
            double lower = 2*x0 - (x0*x0/2.0);
            if( x0 < x ) {return higher - lower;}
            else {return lower - higher;}
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

        double T(double x0, double x) {
            double higher = x*x/2.0;
            double lower = x0*x0/2.0;
            if(x0 < x) {return higher - lower;}
            else {return lower - higher;}
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

        double T(double x0, double x) {
            double higher = Et(x)/(-3.0);
            double lower = Et(x0)/(-3.0);
            if(x0 < x) {return higher - lower;}
            else {return lower - higher;}
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

        double T(double x0, double x) {
            double higher = Et(x)/(2.0);
            double lower = Et(x0)/(2.0);
            if(x0 < x) {return higher - lower;}
            else {return lower - higher;}
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

        double T(double x0, double x) {
            double higher = std::sqrt(M_PI)*A*std::erf(std::sqrt(a_s*(x-z)))/(2.0*std::sqrt(a_s));
            double lower = std::sqrt(M_PI)*A*std::erf(std::sqrt(a_s*(x0-z)))/(2.0*std::sqrt(a_s));
            if(x0 < x) {return higher - lower;}
            else {return lower - higher;}
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

        double T(double x0, double x) {
            double higher = std::sqrt(M_PI)*A*std::erf(std::sqrt(a_b*(x-z)))/(2.0*std::sqrt(a_b));
            double lower = std::sqrt(M_PI)*A*std::erf(std::sqrt(a_b*(x0-z)))/(2.0*std::sqrt(a_b));
            if(x0 < x) {return higher - lower;}
            else {return lower - higher;}
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

double False_Position(XS* xs, double x0, double T_hat, double& x_low, double& x_hi, 
                      double eps, int& counter) {
    // x_low is lower bound (initialy current x position)
    // x_hi is point where particle would cross to new region
    
    double x_bound = x_hi;
    bool mu_positive = x0 < x_hi;
    
    // Lambda function for T_hat - T(x)
    auto F = [&](double y) {
        return (T_hat - xs->T(x0,y));
    };
    
    double F_hi = F(x_hi);
    double F_low = F(x_low);
    counter++;
    counter++;
    double m = (F_hi - F_low)/(x_hi - x_low);
    double x = x_low;
    double x_old = x_hi;
    while(std::abs(x_old - x) > eps) {
        x_old = x;
        x = x_low - (F_low/m);

        if((x > 2.0) or (x < 0.0)) {
            #pragma omp critical
            {
            std::cout << " PROBS\n";
            std::cout << " xi = " << xi << " G = " << G << "\n";
            std::cout << " Pnc = " << Pnc << "\n";
            std::cout << " T_hat = " << T_hat << "\n";
            std::cout << " Mu positive = " << mu_positive << "\n";
            std::cout << " x0 = " << x0 << " xb = " << x_bound << "\n";
            std::cout << " x_low = " << x_low << "\n" << " x_hi = " << x_hi;
            std::cout << "\n x = " << x << " m = " << m << "\n";
            std::cout << " F_low = " << F_low << " F_hi = " << F_hi << "\n";
            std::cin >> x;
            }
        }

        if(std::abs(x_old - x) > eps) {
            break;
        } else {
            counter++;
            
            if(mu_positive) {
                if(F(x) > 0.0) {
                    x_low = x;
                    F_low = F(x_low);
                }
                else {
                    x_hi = x;
                    F_hi = F(x_hi);
                }
            } else {
                if(F(x) < 0.0) {
                    x_low = x;
                    F_low = F(x_low);
                }
                else {
                    x_hi = x;
                    F_hi = F(x_hi);
                }
            }
            m = (F_hi - F_low)/(x_hi - x_low);
        }
    }
    return x;
}

double Newton(XS* xs, double x_orig, double xi, double Pnc, double G, double T_hat, double x0, double xg, double x_low, 
              double x_hi, int& counter) {
    // xg is initial guess from False_Position
    // x0 is current position
    // Lambda function for T_hat - T(x)
    auto F = [&](double y) {
        return (T_hat - xs->T(x0,y));
    };

    //if(T_hat == x0) {std::cout << " SCREAM\n";} 
    double x = xg;
    xg += 0.05;
    double g = 0.0;
    double gp = 0.0;
    while(std::abs(x - xg) > EPS) {
        counter++;
        counter++;
        xg = x;
        g = F(xg);
        gp = xs->Et(xg);
        x = xg - (g/gp);
        if ((x > x_hi) or (x < x_low)) {
            x = False_Position(xs, x_orig, xi, Pnc, G, T_hat, x_low, x_hi, EPS, counter);
            break;           
        }
    }
    return x;
}

double Distance_To_Surf(double x, double mu) {
    if(mu > 0) {return 2.0 - x;}
    else {return x;}
}

void Direct_Sampling(XS* xs) {
    std::cout << "\n Direct Sampling (Newton's Method)\n";
    double escape = 0.0;
    int cnts_sum = 0;

    #pragma omp parallel
    {
        pcg64_unique rng;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double x = 1.0; // All particles start in middle
            double mu = 1.0;
            //if(rand(rng) < 0.5) {mu = -1.0;}
            bool alive = true;
            int cnt = 0;

            while(alive) {
                double ds = Distance_To_Surf(x, mu);
                double Pnc = std::exp(-(xs->T(x,x + (mu*ds)))); // Probability of no collision
                cnt++;
                double G = 1.0 - Pnc;
                double xi = rand(rng);
                if(xi < Pnc) { // Particle escapes
                    alive = false;
                    #pragma omp atomic
                    escape += 1.0;
                } else { // Particle has collision
                    double T_hat = -std::log(1.0 - G*xi);
                    double x_low = x;
                    double x_hi = x + mu*ds;
                    double eps = 0.01;
                    double xg = False_Position(xs, x, xi, Pnc, G, T_hat, x_low, x_hi, eps, cnt);
                    x = Newton(xs, x, xi, Pnc, G, T_hat, x, xg, x_low, x_hi, cnt);

                    // Determine if scatter or absorption
                    if(rand(rng) < Pscatt) { // Scatter occurs
                        if(rand(rng) > 0.5) {mu = 1.0;}
                        else {mu = -1.0;}
                    } else { // Particle absorbed
                        alive = false;
                    }
                }
            }

            #pragma omp atomic
            cnts_sum += cnt;
        }
    }

    double trans = escape / (double)NPART;
    double avg_cnt = (double)cnts_sum/(double)NPART;

    std::cout << " Escaped: " << (int)escape << ", Transmission: ";
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
    
    for(int type = 1; type <= 1; type++) {
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
        //double Tau = crs->T(1.9999,2.0);
        //double Pnc = std::exp(-Tau);
        //double G = 1.0 - Pnc;
        //std::cout << " T = " << Tau << "\n";
        //std::cout << " Pnc = " << Pnc << "\n";
        //std::cout << " G = " << G << "\n";

        //Delta_Tracking(crs);

        //Meshed_Delta_Tracking(crs);

        //Negative_Weight_Delta_Tracking(crs);

        //Meshed_Negative_Weight_Delta_Tracking(crs);
    }
    
    return 0;
}
