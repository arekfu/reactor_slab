#include"random/pcg_random.hpp"
#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<chrono>
#include<omp.h>

using namespace std::chrono;

// General simulation settings
int nparticles = 1000000;
const int nparticles_initial = nparticles;
const int array_len = 2*nparticles;
const int ngens = 105;
const int ninactive = 5;
const int ngroups = 2;
const int nmaterials = 2;

// Entropy mesh
const double Hdx = 0.2;

// Flux tally mesh
const double Fdx = 0.1;

// Survival Biasing parameters
const double w_c = 0.25;
const double w_s = 1.0;

// Surface locations
const double Surfaces[4] = {0.0, 30.0, 70.0, 100.0};

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

// Allocate and assign cross sections
void Cross_Sections(double** &Et, double** &Ea, double** &Ef,
                    double*** &Es, double*** &Es1, double** &v) {

    // Allocate memmory for cross sections
    Et = new double*[2];
    Ea = new double*[2];
    Ef = new double*[2];
     v = new double*[2];
    for(int i = 0; i < ngroups; i++) {
        Et[i] = new double[2];
        Ea[i] = new double[2];
        Ef[i] = new double[2];
         v[i] = new double[2];
    }

    Es = new double**[2];
    Es1 = new double**[2];
    for(int i = 0; i < ngroups; i++) {
        Es[i] = new double*[2];
        Es1[i] = new double*[2];
        for(int j = 0; j < ngroups; j++) {
            Es[i][j] = new double[2];
            Es1[i][j] = new double[2];
        }
    }

    /* Notes on cross section variables
       
       First index = material : 0 = reflector, 1 = core
       Second index = energy group : 0 = fast, 1 = slow

       Scattering Matricies
       First index = material
       Second index = incoming group
       Third index = outgoing group
    */

    // Set cross sections for reflector material
    Et[0][0] = 0.295;
    Et[0][1] = 2.1;

    Ef[0][0] = 0.0;
    Ef[0][1] = 0.0;

    Ea[0][0] = 0.0004;
    Ea[0][1] = 0.02;

    Es[0][0][0] = 0.2456;
    Es[0][0][1] = 0.049;
    Es[0][1][0] = 0.0;
    Es[0][1][1] = 2.08;

    v[0][0] = 0.0;
    v[0][1] = 0.0;

    // Set cross sections for core material
    Et[1][0] = 0.276;
    Et[1][1] = 1.063;

    Ef[1][0] = 0.00339;
    Ef[1][1] = 0.074;

    Ea[1][0] = 0.012;
    Ea[1][1] = 0.121;

    Es[1][0][0] = 0.25;
    Es[1][0][1] = 0.014;
    Es[1][1][0] = 0.0;
    Es[1][1][1] = 0.942;

    v[1][0] = 2.5;
    v[1][1] = 2.5;
}

// Deallocate cross sections
void Deallocate(double** &Et, double** &Ea, double** &Ef, double*** &Es, double*** &Es1, 
                double** &v, double* &x, double* &x_F) {
    for(int i = 0; i < nmaterials; i++) {
        delete[] Et[i];
        delete[] Ea[i];
        delete[] Ef[i];
        delete[]  v[i];

        for(int j = 0; j < ngroups; j++) {
            delete[] Es[i][j];
            delete[] Es1[i][j];
        }
        delete[] Es[i];
        delete[] Es1[i];
    }

    delete[] Et;
    delete[] Ea;
    delete[] Ef;
    delete[] v;
    delete[] Es;
    delete[] Es1;
    delete[] x;
    delete[] x_F;
}

// Allocate position arrays and define initial source
void Generate_Source(double* &x, double* &x_F) {
    x = new double[array_len];
    x_F = new double[array_len];
    for(int i = 0; i < array_len; i++) {
        x[i] = 0.0;
        x_F[i] = 0.0;
    }
    pcg64_unique rnd;
    for(int i = 0; i < nparticles; i++) {
        x[i] = (Surfaces[2]-Surfaces[1])*rand(rnd) + Surfaces[1];
    }
}

// Get distance to closes surface
double Distance_to_Surface(double x, double mu) {
    double d = 10000000000000.0;
    // Loop through all for surfaces
    for(int s = 0; s < 4; s++) {
        double d_to_surf = (Surfaces[s] - x)/mu;
        if( (d_to_surf > 0) && (d_to_surf < d) ) {
            d = d_to_surf;
        }
    }
    return d;
}

// Get material where particle is located
int Get_Material(double x) {
    if( (x >= Surfaces[0]) and (x < Surfaces[1]) ) { // In left reflector
        return 0;
    }
    else if( (x >= Surfaces[1]) and (x < Surfaces[2]) ) { // In core
        return 1;
    }
    else if( (x >= Surfaces[2])  and (x < Surfaces[3]) ) { // Right reflector
        return 0;
    }
    else { // Left to vacuum
        return -1;
    }
}

// Get initial eigenvalue of original source
double Initial_Eigenvalue(double Emaj[2], double** &Et, double** &Ea, double** &Ef,
                          double*** &Es, double*** &Es1, double** &v, double* &pos) {
    // Tally for single generation analog game
    double tally = 0.0;

    #pragma omp parallel
    {
        // Random number generator for each thread
        pcg64_unique rnd;

        // Interate through all particles
        #pragma omp for
        for(int i = 0; i < nparticles; i++) {
            //rnd.advance(125);
            
            // Initialize particle as alive and fast group, in core material
            double x = pos[i];
            bool alive = true;
            int energy = 0;
            int material = 1;

            while(alive) {
                double mu = 2.0*rand(rnd) - 1.0; // Isotropic scattering
                bool virtual_collision = true;
                while(virtual_collision) {
                    double d = -std::log(rand(rnd))/Emaj[energy]; // Distance to collision
                    int new_material = Get_Material(x + (d*mu)); // Material where particle lands
                    if(new_material != -1) { // Particle hasn't left simulation
                        if(rand(rnd) < (Et[new_material][energy]/Emaj[energy])) { // Virtual collision accepted
                            virtual_collision = false;
                            x += d*mu;
                            material = new_material;
                        }
                    }
                    else { // Particle has left simulation, kill it
                        virtual_collision = false;
                        alive = false;
                    }
                }

                if (alive == true) {
                    // Particle now has a collision, what type?
                    if(rand(rnd) < Ea[material][energy]/Et[material][energy]) { // Particle absorbed
                        #pragma omp atomic
                        tally += v[material][energy]*Ef[material][energy]/Ea[material][energy];
                        alive = false;
                    }
                    else { // Particle scattered, to what group?
                        double scatter_total = Es[material][energy][0] + Es[material][energy][1];
                        if(rand(rnd) < Es[material][energy][0]/scatter_total) { energy = 0; }
                        else { energy = 1; }
                    }
                }
            }
        }
    }

    return tally / (double)nparticles;
}

// Write flux to csv file
void Write_Flux(double** &Flux_avg, double** &Flux_sqr_avg, int &nFbins) {
    std::ofstream output;
    output.open("nvc_delta_tracking_flux.csv");

    // Write central x value of each bin
    double x = 0.5*Fdx;
    output << x;
    for(int i = 1; i < nFbins; i++) {
        output << "," << (double)i*Fdx + 0.5*Fdx ;
    }
    output << "\n";

    // Write average fast and thermal fluxes
    output << Flux_avg[0][0];
    for(int i = 1; i < nFbins; i++) {
        output << "," << Flux_avg[0][i];
    }
    output << "\n";

    output << Flux_avg[1][0];
    for(int i = 1; i < nFbins; i++) {
        output << "," << Flux_avg[1][i];
    }
    output << "\n";

    // Write standard detiation for fast and thermal fluxes
    double std = std::sqrt( std::abs( Flux_sqr_avg[0][0] - (Flux_avg[0][0]*Flux_avg[0][0]) ) / (double)(ngens - ninactive - 1) );
    output << std;
    for(int i = 1; i < nFbins; i++) {
        std = std::sqrt( std::abs( Flux_sqr_avg[0][i] - (Flux_avg[0][i]*Flux_avg[0][i]) ) / (double)(ngens - ninactive - 1) );
        output << "," << std;
    }
    output << "\n";

    std = std::sqrt( std::abs( Flux_sqr_avg[1][0] - (Flux_avg[1][0]*Flux_avg[1][0]) ) / (double)(ngens - ninactive - 1) );
    output << std;
    for(int i = 1; i < nFbins; i++) {
        std = std::sqrt( std::abs( Flux_sqr_avg[1][i] - (Flux_avg[1][i]*Flux_avg[1][i]) ) / (double)(ngens - ninactive - 1) );
        output << "," << std;
    }

    output.close();
}

// Perform power series for acctual eigenvalue
void Eigenvalue(double Emaj[2], double** &Et, double** &Ea, double** &Ef, double*** &Es,
                double*** &Es1, double** &v, double* &pos, double* &pos_future, double &k) {
    
    double k_sum = 0.0;
    double k_sqr_sum = 0.0;

    high_resolution_clock::time_point t0; // Timers for FOM
    high_resolution_clock::time_point t;

    // Entropy bin allocation
    int nHbins = std::round((Surfaces[2] - Surfaces[1])/Hdx);
    double* Entropy;
    Entropy = new double[nHbins];
    for(int i = 0; i < nHbins; i++) {
        Entropy[i] = 0.0;
    }

    // Flux bin allocation
    int nFbins = std::round((Surfaces[3] - Surfaces[0])/Fdx);
    double** Flux;
    double** Flux_sum;
    double** Flux_sqr_sum;
    Flux = new double*[2];
    Flux_sum = new double*[2];
    Flux_sqr_sum = new double*[2];
    Flux[0] = new double[nFbins];
    Flux[1] = new double[nFbins];
    Flux_sum[0] = new double[nFbins];
    Flux_sum[1] = new double[nFbins];
    Flux_sqr_sum[0] = new double[nFbins];
    Flux_sqr_sum[1] = new double[nFbins];
    for(int i = 0; i < nFbins; i++) {
        Flux[0][i] = 0.0;
        Flux[1][i] = 0.0;
        Flux_sum[0][i] = 0.0;
        Flux_sum[1][i] = 0.0;
        Flux_sqr_sum[0][i] = 0.0;
        Flux_sqr_sum[1][i] = 0.0;
    }    
    
    // Loop to iterate through all generations
    for(int g = 1; g <= ngens; g++) {
        // Tally for single generation analog game
        double tally = 0.0;
        int nparticles_future = 0;
        //double w_initial = (double)nparticles_initial/(double)nparticles;
        double w_initial = 1.0;

        // Zero entropy bins
        for(int i = 0; i < nHbins; i++) {
            Entropy[i] = 0.0;
        }

        // Zero flux bins
        for(int i = 0; i < nFbins; i++) {
            Flux[0][i] = 0.0;
            Flux[1][i] = 0.0;
        }

        #pragma omp parallel
        {
            // Random number generator for each thread
            pcg64_unique rnd;

            // Interate through all particles
            #pragma omp for
            for(int i = 0; i < nparticles; i++) {
                //rnd.advance(125);

                // Initialize particle as alive and fast group, in core material
                double x = pos[i];
                bool alive = true;
                int energy = 0;
                int material = 1;
                double w = w_initial;

                // Record particle entropy
                int bin = std::floor((x-Surfaces[1])/Hdx);
                #pragma omp atomic
                Entropy[bin] += 1.0/(double)nparticles;

                while(alive) {
                    double mu = 2.0*rand(rnd) - 1.0; // Isotropic scattering
                    // Determine collision distance
                    double d = -std::log(rand(rnd))/Emaj[energy];
                    x += d*mu; // Position updated
                    material = Get_Material(x); // Material updated
                    if(material == -1) { // Check is particle left simulation
                        alive = false;
                    }
                    
                    // Original delta tracking virtual collision block
                    /*bool virtual_collision = true;
                    while(virtual_collision) {
                        double d = -std::log(rand(rnd))/Emaj[energy]; // Distance to collision
                        int new_material = Get_Material(x + (d*mu)); // Material where particle lands
                        if(new_material != -1) { // Particle hasn't left simulation
                            if(rand(rnd) < (Et[new_material][energy]/Emaj[energy])) { // Virtual collision accepted
                                virtual_collision = false;
                                x += d*mu;
                                material = new_material;
                            }
                        }
                        else { // Particle has left simulation, kill it
                            virtual_collision = false;
                            alive = false;
                        }
                    }*/

                    if(alive == true) {
                        #pragma omp atomic
                        //tally += w*v[material][energy]*Ef[material][energy]/Et[material][energy];
                        tally += w*((Et[material][energy]/Emaj[energy]))*v[material][energy]*Ef[material][energy]/Et[material][energy];

                        // Tally flux
                        if(g > ninactive) {
                            bin = std::floor(x/Fdx);
                            #pragma omp atomic
                            //Flux[energy][bin] += w/(Et[material][energy]);
                            Flux[energy][bin] += w*((Et[material][energy]/Emaj[energy]))/(Et[material][energy]);
                        }

                        //int new_particles = floor((w/k)*v[material][energy]*Ef[material][energy]/Et[material][energy] + rand(rnd));
                        int new_particles = floor((w*((Et[material][energy]/Emaj[energy]))/k)*v[material][energy]*Ef[material][energy]/Et[material][energy] + rand(rnd));
                        for(int n = 0; n < new_particles; n++) {
                            #pragma omp critical
                            {
                                pos_future[nparticles_future] = x;
                                nparticles_future += 1;
                            }
                        }
                        
                        // New weight after partial fission
                        // w = w*(1.0 - (Et[material][energy]/Emaj[energy]))*(1.0 - (Ea[material][energy]/Et[material][energy]));
                        w = w*(1.0 - Et[material][energy]/Emaj[energy])*(1.0 - (Ea[material][energy]/Et[material][energy]));

                        // Russian Roulett
                        if(w < w_c) {
                            if(rand(rnd) > (1 - (w/w_s))) {
                                w = w_s;
                            }
                            else{
                                alive = false;
                            }
                        }

                        // Scatter
                        double scatter_total = Es[material][energy][0] + Es[material][energy][1];
                        if(rand(rnd) < Es[material][energy][0]/scatter_total) { energy = 0; }
                        else { energy = 1; }

                    }
                }
            }
        } // End parallel portion, all particles have been transported

        // Print table headers
        if(g == 1) {
            std::cout << " Gen   Keff     Entropy\n";
            std::cout << " ----------------------\n";
        }
        else if(g == ninactive+1) {
            std::cout << " --------------------------------------------------\n";
            std::cout << " Gen   Keff     Entropy | Kavg     STD      FOM    \n";
            std::cout << " --------------------------------------------------\n";
        }

        // Calculate new eigenvalue
        //k = tally / (double)nparticles_initial;
        k = tally / (double)nparticles;

        if(g > ninactive) {
            k_sum += k;
            k_sqr_sum += (k*k);
        }

        // Calculate shannon entropy
        double H = 0.0;
        for(int i = 0; i < nHbins; i++) {
            H += -(Entropy[i]*std::log2(Entropy[i]));
        }

        // Add generations flux to flux tallies
        for(int i = 0; i < nFbins; i++) {
            Flux_sum[0][i] += Flux[0][i]/(double)(ngens - ninactive);
            Flux_sum[1][i] += Flux[1][i]/(double)(ngens - ninactive);
            Flux_sqr_sum[0][i] += (Flux[0][i]*Flux[0][i])/(double)(ngens - ninactive);
            Flux_sqr_sum[1][i] += (Flux[1][i]*Flux[1][i])/(double)(ngens - ninactive);
        }

        // Prepare for next generation
        for(int i = 0; i < nparticles_future; i++) {
            pos[i] = pos_future[i];
        }

        // Print results
        if(g > ninactive) {
            double kavg = k_sum / (g - ninactive);
            double k_sqr_avg = k_sqr_sum / ((double)g - (double)ninactive);
            double std = std::sqrt( std::abs(k_sqr_avg - (kavg*kavg)) / ((double)g - (double)ninactive - 1.0) );
            if(g > ninactive+1) {
                t = high_resolution_clock::now();
                duration<double, std::milli> time_span = t - t0;
                double T = (double)time_span.count() * 0.001;
                double FOM = 1.0/( T * (std/kavg)*(std/kavg) );
                std::cout << " " << std::setfill(' ') << std::setw(4) << g << "  " << std::fixed << std::setprecision(5) << k ;
                std::cout << "  " << std::fixed << std::setprecision(5) << H << " | " << std::fixed << std::setprecision(5);
                std::cout << kavg << "  " << std::fixed << std::setprecision(5) << std << "  " << std::setfill(' ') << std::setw(7) << (int)FOM << "\n";
            }
            else {
                std::cout << " " << std::setfill(' ') << std::setw(4) << g << "  " << std::fixed << std::setprecision(5) << k;
                std::cout << "  " << std::fixed << std::setprecision(5) << H << "\n";
            }
        }
        else {
            std::cout << " " << std::setfill(' ') << std::setw(4) << g << "  " << std::fixed << std::setprecision(5) << k;
            std::cout << "  " << std::fixed << std::setprecision(5) << H << "\n";
        }

        nparticles = nparticles_future;

        if(g == ninactive) { // Get system time for calculating FOM
            t0 = high_resolution_clock::now();
        }
    }
    double kavg = k_sum / (double)(ngens - ninactive);
    double k_sqr_avg = k_sqr_sum / ((double)(ngens - ninactive));
    double std = std::sqrt( std::abs(k_sqr_avg - (kavg*kavg)) / ( (double)(ngens - ninactive - 1) ));
    std::cout << "\n Average k = " << kavg << " +/- " << std << "\n";

    Write_Flux(Flux_sum, Flux_sqr_sum, nFbins);

    // Deallocate arrays
    delete[] Flux[0];
    delete[] Flux[1];
    delete[] Flux_sum[0];
    delete[] Flux_sum[1];
    delete[] Flux_sqr_sum[0];
    delete[] Flux_sqr_sum[1];

    delete[] Flux;
    delete[] Flux_sum;
    delete[] Flux_sqr_sum;
    delete[] Entropy;
}

int main() {
    std::cout << " No Virtual Collision Weighted Delta Tracking Two Group Reflected Slab\n\n";

    // Declare cross section variables
    double**  Et; // Total XS
    double**  Ea; // Absorption XS
    double**  Ef; // Fission XS
    double*** Es; // Scattering XS matrix
    double*** Es1; // Anisotropic scattering XS matrix
    double**  v; // Neutrons per fission

    // Load cross sections
    Cross_Sections(Et, Ea, Ef, Es, Es1, v);

    // Delta tracking majorant cross sections
    double Emaj[2];
    Emaj[0] = 0.296;//0.295;
    Emaj[1] = 2.101;//2.1;

    // Arrays for particle positions in present and future gen
    double* x;
    double* x_F;

    // Allocate position arrays and create source distribution
    Generate_Source(x, x_F);

    // Get initial eigenvalue
    double k = Initial_Eigenvalue(Emaj, Et, Ea, Ef, Es, Es1, v, x);
    std::cout << " Initial multiplication factor: " << k << "\n\n";    

    // Perform iterations for eigen value
    Eigenvalue(Emaj, Et, Ea, Ef, Es, Es1, v, x, x_F, k);
    
    // Deallocate all memory
    Deallocate(Et, Ea, Ef, Es, Es1, v, x, x_F);
    return 0;
}
