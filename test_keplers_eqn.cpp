//test_keplers_eqn.cpp 
//by Joe Hahn, jhahn@spacescience.org, 4 June 2013.
//
//This demonstrates how to use the Orbit.cpp class, and it also tests keplers_eqn() that
//solve Kepler's equation numerically. The output is used to check the numerical
//accuracy and speed of keplers_eqn(). I compile this on my Ubuntu pc via
//g++ Orbit.cpp test_keplers_eqn.cpp -O3 -o test_keplers_eqn

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <fstream>
#include "Orbit.h"

int main()
{
    using namespace std;

    //This is the maximum error allowed when keplers_eqn() is used to solve for the
    //orbit's eccentric anomaly E. Shrink this number to increase numerical accuracy.
    double max_error = 1.0e-15;

    //Solve Kepler's equation N times. 
    int N = 10000000;

    //Save a subset of the N solutions, for plotting and testing purposes.
    int Nsave = 100000;

    //The remaining quanties used below.
    double* e = new double[N];
    double* M = new double[N];
    double* E = new double[N];
    int* iterations = new int[N];
    double x, xtime_sec, rn;
    clock_t t0, t1;
    double pi = Orbit().get_pi();
    Orbit orb;

    //Solve kepler's equation N times. Each orbit is randomly generated with
    //mean anomly M uniformly distributed between -Pi<M<Pi and eccentricity
    //logarithmically distributed so that log10(x) is uniformly distributed over
    //-6<log10(x)<0 where x=1-e so that e ranges over 10^(-6)<e<0.999999.
    srand(1);
    t0 = clock();
    for (int i=0; i<N; i++){

        //Get random number between 0<rn<1, calculate x=1-e and e=1-x.
        rn = (double)rand()/(double)RAND_MAX;
        x = pow(10.0, -6.0*rn);
        e[i] = 1.0 - x;

        //Random mean anomaly M.
        rn = (double)rand()/(double)RAND_MAX;
        M[i] = pi*(2.0*rn - 1.0);

        //Construct the orbit.
        orb = Orbit(1.0, e[i], 0.0, 0.0, 0.0, M[i]);

        //Calculate eccentric anomaly and track number of iterations used by keplers_eqn().
        E[i] = orb.keplers_eqn(max_error);
        iterations[i] = orb.i_ke;
    }
    //Execution time consumed by the main loop.
    t1 = clock();
    xtime_sec = (double)(t1 - t0)/CLOCKS_PER_SEC;
    cout << "execution time (sec) = " << xtime_sec << endl;

    //Store Nsave results in file keplers_eqn.dat.
    ofstream out("keplers_eqn.dat", ios::binary);
    out.write((char *)&xtime_sec, sizeof(xtime_sec));
    out.write((char *)&Nsave, sizeof(Nsave));
    out.write((char *)e, sizeof(e)*Nsave);
    out.write((char *)M, sizeof(M)*Nsave);
    out.write((char *)E, sizeof(E)*Nsave);
    out.write((char *)iterations, sizeof(iterations)*Nsave);
    out.close();

    //Delete large arrays.
    delete[] e;
    delete[] E;
    delete[] M;
}

