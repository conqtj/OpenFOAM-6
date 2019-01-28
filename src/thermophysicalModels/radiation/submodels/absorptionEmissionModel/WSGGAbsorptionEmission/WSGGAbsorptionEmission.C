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
\*---------------------------------------------------------------------------*/

#include "WSGGAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(WSGGAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            WSGGAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::WSGGAbsorptionEmission::WSGGAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    Nq(8),
    w(Nq, 1.0),
    Tref(1200),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    H2O_
    (
        IOobject
        (
            "H2O",
            mesh_.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    CO2_
    (
        IOobject
        (
            "CO2",
            mesh_.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    )
{
    readData();
}

// * * * * * * * * * * * * * Public Functions * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::WSGGAbsorptionEmission::aCont(const label bandi) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0)
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();

    scalarList kg(Nq);

    forAll(a, celli)
    {
        scalar XH2O = std::max(1e-06, H2O_[celli]);
        scalar XCO2 = std::max(1e-06, CO2_[celli]);

        if(!bandi)kg[0]=1e-12;
        else kg[bandi]=KG[bandi]*(XCO2+XH2O);

        a[celli] = kg[bandi];
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::WSGGAbsorptionEmission::eCont(const label bandi) const
{
    const volScalarField& T = thermo_.T();

    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, 0.5)
        )
    );

   return te;
    scalarField& e = te.ref().primitiveFieldRef();

    scalarList kg(Nq);
    scalarList aFg(Nq);

    forAll(e, celli)
    {
        scalar XH2O = std::max(1e-06, H2O_[celli]);
        scalar XCO2 = std::max(1e-06, CO2_[celli]);

        scalar Tg = T[celli];



        kg[0]=1e-12;
	kg[bandi]=KG[bandi]*(XCO2+XH2O);
      
        for(int i=0;i<3;i++)
        {
            aFg[i+1]=0.0;
            for(int j=0;j<4;j++)
            aFg[i+1]+=bij[i][j]*pow(Tg,j);
        }
        
        //clear gas
        aFg[0]=1.0-aFg[1]-aFg[2]-aFg[3];

        e[celli] = kg[bandi]*aFg[bandi];
    }

    return te;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::WSGGAbsorptionEmission::~WSGGAbsorptionEmission()
{}

// * * * * * * * * * * * * * * Read Data * * * * * * * * * * * * * * * * * //

void Foam::radiation::WSGGAbsorptionEmission::readData()
{
    //read coefficients for the WSGG Model 
    string s;

    std::ifstream dataIn("wsgg.xml");

    if(dataIn)
    {
        getline(dataIn,s);
        getline(dataIn,s);

        for(int i=0;i<3;i++)
            dataIn>>KG[i]>>bij[i][0]>>bij[i][1]>>bij[i][2]>>bij[i][3];


    for(int j=0;j<4;j++)
    {
        for(int i=0;i<3;i++)
            dataIn>>CAijk[i][j][0]>>CAijk[i][j][1]>>CAijk[i][j][2]>>CAijk[i][j][3];
    }

    }

    dataIn.close();

    return;
}
