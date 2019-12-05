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
#include<chrono>
#include<iomanip>
#include<vector>
#include<fstream>
#include<iostream>

const double EPS = 1e-6;
const int NPART = 1e7;
const int NBIN = 100;
const double dx = 10.0/(double)NBIN;
const double p = 0.9;
const double p_mshd = 0.9;
const double q = 0.992;
const double q_mshd = 0.992;

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

// FOM vectors
std::vector<std::vector<double>> Coll_bin;
std::vector<std::vector<double>> FOM_DT;
std::vector<std::vector<double>> FOM_NDT;

std::ofstream File;

// Timers for FOM
std::chrono::high_resolution_clock::time_point t_start;
std::chrono::high_resolution_clock::time_point t_end;

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
        Step():XS(std::exp(-2.139)  , 0.48) {
            for(int b = 0; b < NBIN; b++) {
                double x = b*dx;
                Em[b] = Et(x);
                Esmp[b] = p_mshd*Et(x);
            }
        }

        double T(double x) {
            if((x >= 0.0) and (x < 0.1)) {return 0.48*x;}
            else if((x >= 0.1) and x <= 10.0) {return 0.48*0.1 + 0.32*(x - 0.1);}
            else if(x > 10.0) {return T(10.0);}
            else {
                #pragma omp critical
                {
                    std::cout << x<<"\n";
                }
                exit(1);}
        }

        double Et(double x) {
            if((x >= 0.0) and (x < 0.1)) {return 0.48;}
            else if ((x >= 0.1) and (x < 10.0)) {return 0.32;}
            else {exit(1);}
        }

        double dEt(double x) {return 0.0;}
}; // Step

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

void Delta_Tracking(XS* xs) {
    std::cout << "\n Delta Tracking\n";

    int cnts_sum = 0;
    t_start = std::chrono::high_resolution_clock::now();
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

                    int bin = std::floor(x/dx);
                    #pragma omp atomic
                    Coll_bin[bin][0] += 1.0;
                    #pragma omp atomic
                    Coll_bin[bin][1] += 1.0;
                    virtual_collision = false;
                }
            }
        }
    }
    xs_evals += cnts_sum;
    t_end = std::chrono::high_resolution_clock::now();
}

void Meshed_Delta_Tracking(XS* xs) {
    std::cout << "\n Meshed Delta Tracking\n";

    int cnts_sum = 0;
    int virtual_cnt_sum = 0;
    int bin_cnt_sum = 0;
    
    t_start = std::chrono::high_resolution_clock::now();
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
    t_end = std::chrono::high_resolution_clock::now();
}

void Negative_Weight_Delta_Tracking(XS* xs) {
    std::cout << "\n Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    double sign_change = 0.0;
    double Esamp = p*(xs->Emax);
    t_start = std::chrono::high_resolution_clock::now();
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

                    int bin = std::floor(x/dx);
                    #pragma omp atomic
                    Coll_bin[bin][0] += 1.0;
                    #pragma omp atomic
                    Coll_bin[bin][1] += 1.0;

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
    t_end = std::chrono::high_resolution_clock::now();
}

void Meshed_Negative_Weight_Delta_Tracking(XS* xs) {
    std::cout << "\n Meshed Negative Weight Delta Tracking\n";

    int cnts_sum = 0;
    int bin_cnt_sum = 0;
    double sign_change = 0.0;
    t_start = std::chrono::high_resolution_clock::now();
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
    t_end = std::chrono::high_resolution_clock::now();
}

void Output(int type) {
    double collide_avg = collide / (double)NPART;
    double collide_sqr_avg = collide_sqr / (double)NPART;
    double collide_std = std::sqrt(std::abs(collide_avg*collide_avg - 
                          collide_sqr_avg)/((double)NPART - 1.0));
    double escape_avg = escape / (double)NPART;
    double escape_sqr_avg = escape_sqr / (double)NPART;
    double escape_std = std::sqrt(std::abs(escape_avg*escape_avg - 
                           escape_sqr_avg)/((double)NPART - 1.0));

    double avg_xs_evals = xs_evals / (double)NPART;

    std::chrono::duration<double, std::milli> delta_T = t_end - t_start; 
    double T_ms = (double)delta_T.count() * 0.001;

    // For output file
    for(int i = 0; i < NBIN; i++) {
        double coll_bin_avg = Coll_bin[i][0] / (double)NPART;
        double coll_bin_sqr_avg = Coll_bin[i][1] / (double)NPART;
        double coll_bin_std = std::sqrt(std::abs(coll_bin_avg*coll_bin_avg - 
                    coll_bin_sqr_avg)/((double)NPART - 1.0));
        double FOM_cnts = 1.0 / (avg_xs_evals * coll_bin_std * coll_bin_std);
        double FOM_T = 1.0 / (T_ms * coll_bin_std * coll_bin_std);

        if(type == 1) {
            FOM_DT[i][0] = FOM_cnts;
            FOM_DT[i][1] = FOM_T;
        } else {
            FOM_NDT[i][0] = FOM_cnts;
            FOM_NDT[i][1] = FOM_T;
        }
    }

    double FOM_col = 1.0 / (avg_xs_evals * collide_std * collide_std);
    double FOM_esc = 1.0 / (avg_xs_evals * escape_std * escape_std);

    double FOM_T_col = 1.0 / (T_ms * collide_std * collide_std);
    double FOM_T_esc = 1.0 / (T_ms * escape_std * escape_std);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << " Colsn. Rate: " << collide_avg << " +/- " << collide_std;
    std::cout << ", Trans. Rate: " << escape_avg << " +/- " << escape_std << "\n";
    std::cout << " Average XS/Integration Evaluations: " << avg_xs_evals << "\n";
    std::cout << std::scientific;
    std::cout << " FOM_coll = " << FOM_col << ", FOM_escp = " << FOM_esc << "\n";
    std::cout << " FOM_T_col = " << FOM_T_col << ", FOM_T_escp = " << FOM_T_esc << "\n\n";
}

int main() {
   /* std::cout << "\n NParticles = " << NPART << ", NBins = " << NBIN << "\n\n";
    
    for(int type = -1; type <= 6; type++) {
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

    }*/

    Step xs = Step();
    XS* crs = &xs;

    for(int i = 0; i < NBIN; i++) {
        std::vector<double> box;
        box.push_back(0.0);
        box.push_back(0.0);
        Coll_bin.push_back(box);
        FOM_DT.push_back(box);
        FOM_NDT.push_back(box);
    }

    collide = 0.0;
    collide_sqr = 0.0;
    escape = 0.0;
    escape_sqr = 0.0;
    xs_evals = 0.0;
    Delta_Tracking(crs);
    Output(1);

    for(int i = 0; i < NBIN; i++) {
        Coll_bin[i][0] = 0.0;
        Coll_bin[i][1] = 0.0;
    }

    collide = 0.0;
    collide_sqr = 0.0;
    escape = 0.0;
    escape_sqr = 0.0;
    xs_evals = 0.0;
    Negative_Weight_Delta_Tracking(crs);
    Output(0);

    File.open("FOM_output.txt");
    File << p << " , " << q << "\n";

    File << FOM_DT[0][0];
    for(int i = 1; i < NBIN; i++) {
        File << " , " << FOM_DT[i][0];
    }
    File << "\n";

    File << FOM_DT[0][1];
    for(int i = 1; i < NBIN; i++) {
        File << " , " << FOM_DT[i][1];
    }
    File << "\n";
    
    File << FOM_NDT[0][0];
    for(int i = 1; i < NBIN; i++) {
        File << " , " << FOM_NDT[i][0];
    }
    File << "\n";

    File << FOM_NDT[0][1];
    for(int i = 1; i < NBIN; i++) {
        File << " , " << FOM_NDT[i][1];
    }
    
    File.close();
    return 0;
}
