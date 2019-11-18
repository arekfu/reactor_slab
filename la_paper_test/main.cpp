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

const double EPS = 10e-6;
const int NPART = 10e6;
const double dx = 0.5;
const int NBIN = 9;

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

};

// Derived classes for each type of cross section
class Lin_Decrease : public XS {
    public:
        Lin_Decrease():XS(std::exp(-2.0), 2.0) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
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
        Gauss_Sharp():XS(0.9317314, 2.0/std::sqrt(2.0*M_PI)) {
            int bin_peak = std::floor(1.0/dx);
            double x;
            for(int b = 0; b < NBIN; b++) {
                if(b < bin_peak) {
                    x = b*dx + dx;
                }
                else if(b > bin_peak) {
                    x = b*dx;
                }
                else {
                    x = 1.0;
                }
                Em[b] = Et(x);
            }
        }

        double T(double x) {
            if((x >= 0.0) and (x <= 2.0)) {
                return 0.0353553*(std::erf(0.025*(800.0*x - 800.0)) - 
                                 std::erf(-20.0));
            }
            else if(x > 2.0) {return T(2.0);}
            else {exit(1);}
        }

        double Et(double x) {
            double ar = (x - 1.0)/0.05;
            return (2.0/std::sqrt(2.0*M_PI))*std::exp(-ar*ar);
        }

        double dEt(double x) {
            double ar = (x - 1.0)/0.05;
            return (2.0/std::sqrt(2.0*M_PI))*(-2*ar*(1.0/0.05))*std::exp(-ar*ar);
        }
};

class Gauss_Broad : public XS {
    public:
        Gauss_Broad():XS(0.303686, 2.0/std::sqrt(2.0*M_PI)) {
            int bin_peak = std::floor(1.0/dx);
            double x;
            for(int b = 0; b < NBIN; b++) {
                if(b < bin_peak) {
                    x = b*dx + dx;
                }
                else if(b > bin_peak) {
                    x = b*dx;
                }
                else {
                    x = 1.0;
                }
                Em[b] = Et(x);
            }
        }

        double T(double x) {
            if(x <= 2.0) {
                return (std::erf(1.0) - std::erf(1.0 - x))/std::sqrt(2.0);
            }
            else {return T(2.0);}
        }

        double Et(double x) {
            double ar = (x - 1.0);
            return (2.0/std::sqrt(2.0*M_PI))*std::exp(-ar*ar);
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

double Newton(XS* xs, pcg64_unique& rng, int& counter) {
    double xi = rand(rng);
    if(xi < xs->Pnc){
        counter = 1;
        return 2.1;
    }
    else {
        double T_hat = -std::log(1.0 - (xs->G)*xi);
        //std::cout << " Pnc = " << xs->Pnc << "\n";
        //std::cout << " G = " << xs->G << "\n";
        //std::cout << " T_hat = " << T_hat << "\n";
        int n = 0;
        double s2;
        if(xs->Et(0.0) == 0.0) {s2 = 2.1;}
        else {
            s2 = T_hat/(xs->Et(0.0));
            if(s2 > 2.0) {s2 = 1.0;}
        }
        double s1 = s2 + 20;
        double g, gp;
        //std::cout << " s1 = " << s1 << "\n s2 = " << s2 << "\n";
        while(std::abs(s1 - s2) > EPS) {
            n++;
            s1 = s2;
            g = T_hat - (xs->T(s1));
            gp = -(xs->Et(s1));
            s2 = s1 - (g/gp);
            //std::cout << " s1 = " << s1 << "\n s2 = " << s2 << "\n";
        }
        counter = n;
        return s2;
    }
}

void Direct_Sampling(XS* xs) {
    std::cout << "\n Direct Sampling (Newton's Method)\n";
    double escape = 0.0;
    double collide = 0.0;
    int cnts_sum = 0;

    #pragma omp parallel
    {
        pcg64_unique rng;
        int cnt;
        #pragma omp for
        for(int n = 0; n < NPART; n++) {
            double x = Newton(xs, rng, cnt);
            #pragma omp atomic
            cnts_sum += cnt;
            if(x > 2.0) {
                #pragma omp atomic
                escape += 1.0;
            }
            else {
                #pragma omp atomic
                collide += 1.0;
            }
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
                cnt++;
                if(d_bin < d) {
                    bin_cnt++;
                    d = d_bin + 10e-5;
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
    }
    
    return 0;
}
