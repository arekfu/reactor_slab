/*
 *  Test of "Woodcock tracking with arbitrary sampling cross
 *  section using negative weights".
 *
 * */

#include"random/pcg_random.hpp"

#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>

const double L = 10.0;
const double dx = 10.0/100.0;
const int N = L/dx;
const double p = 0.7;
const double q = 0.3;
const int NPART = 10000000;

// Random number generation functions
double rand(pcg64_unique& rng) {
    uint32_t range = 1000000000;
    uint32_t mask = ~uint32_t(0);
    --range;
    mask >>= __builtin_clz(range|1);
    uint32_t x;
    do {
        x = rng() & mask;
    } while (x > range);
    return double(x)/1000000000.0;
}

int main() {

    // Cross section in each area
    std::vector<double> E;
    for(int n = 0; n < N; n++) {
        if(n == 0) {E.push_back(0.48);}
        else {E.push_back(0.32);}
    }

    // Theoretical transmission
    double E_avg = 0.0;
    for(int n = 0; n < N; n++) {
        E_avg += E[n];
    }
    E_avg /= (double)N;
    std::cout << " Theoretical Transmission = " << std::exp(-E_avg*L) << "\n";

    double Esamp = p*E[0];
    
    // Make empty bins for tally
    std::vector<double> NVC, Flux;
    for(int n = 0; n < N; n++) {
        NVC.push_back(0.0);
        Flux.push_back(0.0);
    }

    // Transport!
    double escape = 0.0;
    pcg64_unique rng;
    for(int p = 0; p < NPART; p++) {
        bool alive = true;
        double x = 0.0;
        double w = 1.0;
        while(alive) {
            double d = -std::log(rand(rng))/Esamp;
            x += d;
           
            if(x > L) {
                alive = false;
                escape += w;
            } else {
                int bin = std::floor(x/dx);
                // Check if real
                if(rand(rng) < q) {
                    alive = false;
                    Flux[bin] += w;
                } else {
                    w = w*((1.0 - (E[bin]/Esamp))/(1.0 - q));
                    NVC[bin] += 1.0;
                }
            }
        }
    }

    std::cout << " Transmission = " << escape/(double)NPART << "\n";

    std::ofstream output("data.csv");
    output << Flux[0];
    for(int n = 1; n < N; n++) {
        output << "," << Flux[n];
    }
    output << "\n";
    output << NVC[0];
    for(int n = 1; n < N; n++) {
        output << "," << NVC[n];
    }
    output.close();

    return 0;
}
