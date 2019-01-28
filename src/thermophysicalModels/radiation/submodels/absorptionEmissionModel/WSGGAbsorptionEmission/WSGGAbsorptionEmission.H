/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
Class
    Foam::radiation::WSGGAbsorptionEmission
Description
    Class for constructing the weighted sum of gray gases model for inhomogeneous
    medium and wall reflection
    e=sum[ai(T){1-exp(-kiPas)}]
    this coefficients satisfy a range of path lengths and temperatures values
 
    Ref: Smith, Shen and Friedman 1982
    Only the weights ai(T) are allowed to vary with temperature
    ai=sum{cij(T/Tref)^j-1}    j=0->3
    Three grey gases and one transparent gas
    Gas number zero is the transperernt gas
    The planck mean absorption coefficient is used for particles, which is 
    added to the WSGG absorption coefficient
    The planck mean coefficients scales linearly with number of particles.
SourceFiles
    WSGGAbsorptionEmission.C
\*---------------------------------------------------------------------------*/

#ifndef WSGGAbsorptionEmission_H
#define WSGGAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "fluidThermo.H"

#include <fstream>

#include "spectralMaths.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class WSGGAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class WSGGAbsorptionEmission
:
    public absorptionEmissionModel
{
private:

    // Private data
    scalar KG[4];
    scalar bij[3][4];
    scalar CAijk[3][4][4];

    //Maths Variable 
    label Nq;
    //w is the weight
    scalarList w;

    //reference state of the gas
    scalar Tref;
    
    //- Thermo package
    const fluidThermo& thermo_;

    //H2O & CO2 fields
    volScalarField H2O_, CO2_;

    //functions
        //Read Data
        void readData();

public:

    //- Runtime type information
    TypeName("WSGGAbsorptionEmission");


    // Constructors
    WSGGAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~WSGGAbsorptionEmission();


    //- Absorption coefficient for continuous phase
    tmp<volScalarField> aCont(const label bandi = 0) const;

    //- Emission coefficient for continuous phase
    tmp<volScalarField> eCont(const label bandi = 0) const;

    inline bool isGrey() const
    {
        return false;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
