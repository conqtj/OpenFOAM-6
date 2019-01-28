#ifndef spectralMaths_H
#define spectralMaths_H

/* This is a library for the spectral model k-distribution
It generates sets of k-values for a set of weights and g-values using
Gaussian quadrature. The number of quadrature points used is normally set
as 8.
There is also an interpolation routine provided
*/

#include "absorptionEmissionModel.H"
#include "fluidThermo.H"

namespace Foam
{
//function to generate Gaussian Quadrature and weights as per Chebyshave pol.
void quadgen(int, scalarList&, scalarList&);

//Quadrature generation as per Gauss Chebyshev scheme
void gausscheb2(scalarList&, scalarList&, int);

//generating k-values based on power law scheme with given power
scalarList kPowerLaw(scalar,scalar, scalar, int);

// linear interpolation of k/g values from 64 to 8 points 
scalarList
linearInterpMono
(
    int,
    const scalarList&,
    const scalarList&,
    int,
    const scalarList&
);
}

#endif
