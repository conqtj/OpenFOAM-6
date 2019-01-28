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
#include <fstream>
#include "kDistAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(kDistAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            kDistAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::kDistAbsorptionEmission::kDistAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    nop(64),
    kmin(1.e-6),
    kmax(10.0),
    Nq(8),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    H2O_
    (
        IOobject
        (
            "H2O",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    CO2_
    (
        IOobject
        (
            "CO2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    readData();
    k1.resize(nop);
    k1=kPowerLaw(kmin, kmax, 0.1,nop);
    w.resize(Nq);
    g.resize(Nq);
    quadgen(Nq, w, g);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::kDistAbsorptionEmission::~kDistAbsorptionEmission()
{}

// * * * * * * * * * * * * * Public Functions * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::kDistAbsorptionEmission::aCont(const label bandi) const
{
    const volScalarField& T = thermo_.T();

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("a", dimless/dimLength, 0)
        )
    );

 
    setRefState();

    scalarField& a = ta.ref().primitiveFieldRef();

    scalarList g1(nop);
    scalarList kg(Nq);

    forAll(a, celli)
    {
        scalar XH2O = std::max(1e-06, H2O_[celli]);
        scalar XCO2 = std::max(1e-06, CO2_[celli]);

        scalar Tg = T[celli];

        //local absorption coefficient 
        g1=fskdistmix(XCO2,XH2O,Tg,Tp,k1,2);

        //interpolate from 64 points to 8 point standard Gaussian Quadrature scheme
        kg=linearInterpMono(nop,g1,k1,Nq,g);

        scalar kgw = 0;
        forAll(kg, point)
        {
            kgw += kg[point]*w[point];
        }

        a[celli] = kgw;
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::kDistAbsorptionEmission::eCont(const label bandi ) const
{
    const volScalarField& T = thermo_.T();

    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("e", dimless/dimLength, 0)
        )
    );

    setRefState();

    scalarField& e = te.ref().primitiveFieldRef();

    scalarList g1(nop);
    scalarList g2(nop);
    scalarList kg(Nq);
    scalarList aFg(Nq);
    
    g2 = fskdistmix(XCO2p,XH2Op,Tp,Tp,k1,2);

    forAll(e, celli)
    {
        scalar XH2O = std::max(1e-06, H2O_[celli]);
        scalar XCO2 = std::max(1e-06, CO2_[celli]);

        scalar Tg = T[celli];

        // generate k* and a-function for each zone 
        g1=fskdistmix(XCO2p,XH2Op,Tp,Tg,k1,2);

        //weight function for the k-distribution  
        scalarList aF2=aFunction(g1,g2,nop);

        //local absorption coefficient 
        g1=fskdistmix(XCO2,XH2O,Tg,Tp,k1,2);

        //interpolate from 64 points to 8 point standard Gaussian Quadrature scheme
        kg=linearInterpMono(nop,g1,k1,Nq,g);
        aFg=linearInterpMono(nop,g2,aF2,Nq,g);

        scalar kgw = 0;
        forAll(kg, point)
        {
            kgw += kg[point]*w[point];
        }
        scalar aFgw = 0;
        forAll(aFg, point)
        {
            aFgw += aFg[point]*w[point];
        }

        e[celli] = kgw*aFgw;
    }
    return te;
}

// * * * * * * * * * * * * * * Read Data * * * * * * * * * * * * * * * * * //

void Foam::radiation::kDistAbsorptionEmission::readData()
{
  string s;
  //read data for CO2 k-Distribution Correlation
  std::ifstream dataIn;
  dataIn.open("fskco2.xml");
  dataIn>>s;

  for(label l=0;l<4;l++)
  {
    for(label m=0;m<4;m++)
    {
      for(label n=0;n<4;n++)
        dataIn>>dlmn[l][m][n]>>s;
    }
  } 
  dataIn.close();

  // readData for H20

  dataIn.open("fskh2o.xml");

  dataIn>>s;

  for(label l=0;l<4;l++)
  {
    for(label m=0;m<4;m++)
    {
      for(label n=0;n<4;n++)
        dataIn>>almn[l][m][n]>>s;
    }
  }


  dataIn>>s;

  for(label l=0;l<3;l++)
  {
    for(label m=0;m<3;m++)
    {
      for(label n=0;n<2;n++)
        dataIn>>blmn[l][m][n]>>s;
    }
  }

  dataIn.close();

  return ;
}

// * * * * * * * * * * * * * * Set Reference * * * * * * * * * * * * * * * //

void Foam::radiation::kDistAbsorptionEmission::setRefState() const
{
    bool refStateFlag=false;

    if(refStateFlag)
    { 
        Tp=1500.;
        XCO2p=0.1;
        XH2Op=0.1;
    }
   else
    {
        // an average of the entire volume field is to be taken
    	scalar totalV=0.;
        Tp=0.0;
        XCO2p=1e-06;
        XH2Op=1e-06;
    
	    const volScalarField& T = thermo_.T();

		forAll(mesh_.cells(), celli)
        {
            Tp += T[celli]*mesh_.V()[celli];
            XCO2p += CO2_[celli]*mesh_.V()[celli];
            XH2Op += H2O_[celli]*mesh_.V()[celli];
            totalV += mesh_.V()[celli];   
        }

        Tp/=totalV;
        XCO2p/=totalV;
        XH2Op/=totalV;
        XCO2p=std::max(1e-06,XCO2p);
        XH2Op=std::max(1e-06,XH2Op);    
    } 

    return ;
}

// * * * * * * * * * * * * * * kDistFunctions * * * * * * * * * * * * * * * //

Foam::scalarList
Foam::radiation::kDistAbsorptionEmission::aFunction(const scalarList& g1,const scalarList & g2,label n) const
{
  scalarList a(n);

  for (label i=1;i<n-1;i++)
  {
    a[i]=(g1[i]-g1[i+1]+1.0e-12)/(g2[i]-g2[i+1]+1.0e-12);
  }

  a[n-1]=(g1[n-2]-g1[n-1]+1.0e-12)/(g2[n-2]-g2[n-1]+1.0e-12);
  a[0]=(g1[0]-g1[1]+1.0e-12)/(g2[0]-g2[1]+1.0e-12);
       
  return a;
}

Foam::scalarList
Foam::radiation::kDistAbsorptionEmission::fskdistmix
(scalar xCO2,scalar xH2O,scalar Tg,scalar Tp,const scalarList & k,label mix) const
{
    /* Adapted from M. F. Modest 
  This function calculates full-spectrum k-distribution of a CO2 and H2O
  mixture   using full-spectrum k-distribution of each individule species.
  The total pressure of the mixture is assumed to be 1 bar.
  Input:
       xCO2: CO2 mole fraction
       xH2O: H2O mole fraction
       Tg: (local) gas temperature 
       Tp: Planck function temperature
        m: integer to specify mixing model
          mix = 1: superposition
          mix = 2: multiplication
          mix = 3: uncorrelated mixture (Modest and Riazzi)
       k: units are in m-1   
    */  

  scalar kp;
  scalarList g(nop);
  scalarList gCO2(nop),gH2O(nop);

  //k is multiplied by 0.01 to convert to in cm-1; unit used in correlations
  for(label ik=0;ik<nop;ik++)
  {
    kp=k[ik]*1.0e-02/xCO2;
    gCO2[ik]=fskDco2(Tg, Tp, kp);
  }

  for(label ik=0;ik< nop;ik++)
  {
    kp=k[ik]*1.0e-02/xH2O;
    gH2O[ik]=fskDh2o(Tg, Tp, kp, xH2O);
  }


  switch(mix)
  {
    case 1: //superposition
      for(label i=0;i<nop;i++)
        g[i]=gH2O[i]+gCO2[i]-1.0;
      break;
    case 2: // multiplication 
      for(label i=0;i<nop;i++)
        g[i]=gH2O[i]*gCO2[i];
      break;
  }

  return g;

} 

Foam::scalar Foam::radiation::kDistAbsorptionEmission::fskDco2(scalar Tg,scalar Tb,scalar absco) const
{
  /*
     Modest & Mehta correlation for CO2
      Input :
      	Tg : local gas temperature in Kelvin
      	Tb : Planck function temperature in Kelvin
      	absco : pressure based absorption coefficient 
    Output variables:
   	gcal : corresponding g value
  */
  label l,m,n;
  scalar bigp,ksai,gcal;

  ksai=log10(absco);
  bigp=0.;
  for(l=0;l<=3;l++)
  {
    for(m=0;m<=3;m++)
    {
      for(n=0;n<=3;n++)
        bigp=bigp+dlmn[l][m][n]*pow(Tg/1000.,n)*pow(Tb/1000.,m)*pow(ksai,l);
    } 
  } 
  gcal=0.5*tanh(bigp)+0.5;
  return gcal;
}

Foam::scalar Foam::radiation::kDistAbsorptionEmission::fskDh2o(scalar Tg,scalar Tb,scalar k,scalar x) const
{
  /*
      Modest & Singh correlation for H2O
      Input :
     	Tg : gas temperature
     	Tb : Planck function temperature
     	k  : pressure based absorption coefficient 
     	x  : mole fraction of water
      Output variables:
     	gcal: correspond g value
 */ 

  label l,m,n;
  scalar bigp,xi,xis,gcal;

  xi=log10(k);//  ! log 10
  xis=0.;  //! if x=0, xis=0 because of x**(m+1)
  for(l=0;l<=2;l++)
  {
    for(m=0;m<=2;m++)
    {
      for(n=0;n<=1;n++) 
        xis=xis+blmn[l][m][n]*pow(Tg/1000.0,l)*pow(xi,n)*pow(x,m+1);       
    }
  }

  bigp=0.;
  for(l=0;l<=3;l++)
  { 
    for(m=0;m<=3;m++)
    {
      for(n=0;n<=3;n++)
      { 
        bigp=bigp+almn[l][m][n]*pow(Tg/1000.0,n)*pow(Tb/1000.,m)*pow(xi+xis,l);
      }
    }
  }
  gcal=0.50*tanh(bigp)+0.50;
   
  return gcal;  
}
