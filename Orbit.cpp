//Orbit.cpp
//by Joe Hahn, jhahn@spacescience.org, 4 June 2013.
//
//These are the methods that are used by the Orbits class. This class
//describes an orbit using six standard orbit elements: a=semimajor axis, e=eccentricity,
//I=inclination, O=longitude of ascending node, w=longitude of periapse, and
//M=mean anomaly, with the angles I,O,w,M in units of radians in the inverval
//-Pi <= I,O,w,M <= Pi. Currently this class is rather limited since the only
//useful method here is keplers_eqn() which solves Kepler's equation numerically
//for an elliptic orbit having 0<=e<1. Inputs are the orbit's eccentricity e and
//mean anomaly M, and the output is the eccentric anomaly E in radians plus
//optional output i_ke=number of iterations used by keplers_eqn(). See test_keplers_eqn.cpp
//for example code that shows how to use keplers_eqn(). I intend to come back to this
//class to make it more useful by adding methods that will advance a particle along
//its orbit by some timestep dt, and methods that will convert orbit elements into
//cartesian coordinates and velocities, and vise versa.

#include "Orbit.h"
#include <iostream>
#include <cmath>

//The default constructor creates and orbit whose orbit elements are all zero.
Orbit::Orbit()
{
    a = e = I = O = w = M = 0;
}

//The primary constructor creates and orbit having the requested orbit elements.
Orbit::Orbit(double ac, double ec, double Ic, double Oc, double wc, double Mc)
{
    a = ac;
    e = ec;
    I = Ic;
    O = Oc;
    w = wc;
    M = Mc;
}

//This method returns Pi.
double Orbit::get_pi() { return atan(1.0)*4.0; }

//This method solves Kepler's equation for the eccentric anomaly E with an allowed
//error=max_error. Inputs are the eccentricity e (which has 0<=e<1 since the orbit
//is an ellipse), the mean anomaly M, and max_error=the allowed uncertainty in E.
//The output is the calculated eccentric anomaly E, keeping in mind that the true
//value in eccentric anomaly is somewhere in the interval E+-error_max. This
//algorithm is detailed in Danby's textbook, Fundamentals of Celestial Mechanics,
//2nd edition.
double Orbit::keplers_eqn(double max_error)
{
    using namespace std;
    int i;
    int i_max = 15;
    double pi, sgn, E, f, error, es, ec, df, ddf, dddf, d1, d2, d3;

    //Put M in the interval -Pi < M < Pi.
    pi = get_pi();
    M = fmod(M, 2*pi);
    if (M >  pi) M = M - 2.0*pi;
    if (M < -pi) M = M + 2.0*pi;

    //Construct the initial solultion.
    sgn = 1.0;
    if (sin(M) < 0.0) sgn = -1.0;
    E = M + sgn*(0.85)*e;

    //Solve kepler's equation iteratively to improve the solution E.
    error = 1.0;
    for(i = 0; i < i_max; i++){
      es = e*sin(E);
      ec = e*cos(E);
      f = E - es - M;
      error = fabs(f);
      if (error < max_error) break;
      df = 1.0 - ec;
      ddf = es;
      dddf = ec;
      d1 = -f/df;
      d2 = -f/(df + d1*ddf/2.0);
      d3 = -f/(df + d2*ddf/2.0 + d2*d2*dddf/6.0);
      E = E + d3;
    }

     //Warn if solution did not converge.
     if (error > max_error)
       std::cout << "***Warning*** Orbit::keplers_eqn() failed to converge***\n";

    //Update the number of iterations used, i_ke, and return the solution.
    i_ke = i;
    return E;
}

