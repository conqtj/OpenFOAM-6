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
    Foam::radiation::kDistAbsorptionEmission
Description
    Class for generation of k-distribution and other related parameters such as
    a-functions and scaling factors
SourceFiles
    kDistAbsorptionEmission.C
\*---------------------------------------------------------------------------*/

#ifndef kDistAbsorptionEmission_H
#define kDistAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "fluidThermo.H"

#include "spectralMaths.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class kDistAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class kDistAbsorptionEmission
:
    public absorptionEmissionModel
{
private:

    // Private data
    label nop;
    scalar dlmn[4][4][4];
    scalar almn[4][4][4];
    scalar blmn[3][3][2];

    scalarList k1;
    scalar kmin, kmax;

    //Maths Variable 
    label Nq;
    // g is the quadrature point and w is the weight
    mutable  scalarList w, g;
    
    //reference state of the gas
    mutable   scalar Tp,XCO2p,XH2Op;
    
    //- Thermo package
    const fluidThermo& thermo_;

    //H2O & CO2 fields
    mutable volScalarField H2O_, CO2_;
    
    //functions
        //Read Data
        void readData();
        //Set Reference
        void setRefState() const;
        //kDist Functions
        scalarList aFunction(const scalarList&, const scalarList&, label) const;
        scalarList fskdistmix(scalar,scalar,scalar,scalar,const scalarList&,label) const;
        scalar fskDco2(scalar, scalar, scalar) const;
        scalar fskDh2o(scalar, scalar, scalar, scalar) const;

public:

    //- Runtime type information
    TypeName("kDistAbsorptionEmission");


    // Constructors
    kDistAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~kDistAbsorptionEmission();


    //- Absorption coefficient for continuous phase
    tmp<volScalarField> aCont(const label bandi = 0) const;

    //- Emission coefficient for continuous phase
    tmp<volScalarField> eCont(const label bandi = 0) const;

    inline bool isGrey() const
    {
        return false;
    }



    friend  void quadgen(int, scalarList, scalarList);

    //Quadrature generation as per Gauss Chebyshev scheme
    friend void gausscheb2(scalarList&,scalarList&, int);

    //generating k-values based on power law scheme with given power
    friend scalarList kPowerLaw(scalar,scalar, scalar, int);

    // linear interpolation of k/g values from 64 to 8 points 
    friend scalarList linearInterpMono
    (
        int,
        const scalarList&,
        const scalarList&,
        int,
        const scalarList &
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
