//Orbit.h
//This header provides the prototypes for the Orbit class defined in Orbit.cpp.

//The #ifndef...#endif stuff prevents this header from being compiled more than once.
#ifndef ORBIT_H_
#define ORBIT_H_

class Orbit{

private:

public:

    //These are the usual orbit elements: semimajor axis a, eccentricity e, inclination I,
    //longitude of ascending node O, longitude of periapse w, and mean anomaly M.
    //Angles I,O,w,M are in radians.
    double a, e, I, O, w, M;

    //Default constructor.
    Orbit();

    //Primary constructor.
    Orbit(double ac, double ec, double I, double Oc, double wc, double Mc);

    //This method solves Kepler's equation for the eccentric anomaly E, which is the
    //return value, with i_ke = number of iterations used by keplers_eqn().
    double keplers_eqn(double max_error);
    int i_ke;

    //This method returns Pi.
    double get_pi();

};
#endif

