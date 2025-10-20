/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "I2GGU.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(I2GGU, 0);
    addToRunTimeSelectionTable(functionObject, I2GGU, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::I2GGU::calc()
{
    if (foundObject<volVectorField>(fieldName_))
    {
        tmp<volVectorField> tU = lookupObject<volVectorField>(fieldName_);
        const volVectorField& U = tU();
        const volVectorField GU0 = fvc::grad(U.component(0));
        const volTensorField GGU0 = fvc::grad(GU0);
        const volVectorField GU1 = fvc::grad(U.component(1));
        const volTensorField GGU1 = fvc::grad(GU1);
        const volVectorField GU2 = fvc::grad(U.component(2));
        const volTensorField GGU2 = fvc::grad(GU2);
        const volScalarField GGU000 = GGU0.component(0*3+0);
        const volScalarField GGU001 = GGU0.component(0*3+1);
        const volScalarField GGU002 = GGU0.component(0*3+2);
        const volScalarField GGU010 = GGU0.component(1*3+0);
        const volScalarField GGU011 = GGU0.component(1*3+1);
        const volScalarField GGU012 = GGU0.component(1*3+2);
        const volScalarField GGU020 = GGU0.component(2*3+0);
        const volScalarField GGU021 = GGU0.component(2*3+1);
        const volScalarField GGU022 = GGU0.component(2*3+2);
        const volScalarField GGU100 = GGU1.component(0*3+0);
        const volScalarField GGU101 = GGU1.component(0*3+1);
        const volScalarField GGU102 = GGU1.component(0*3+2);
        const volScalarField GGU110 = GGU1.component(1*3+0);
        const volScalarField GGU111 = GGU1.component(1*3+1);
        const volScalarField GGU112 = GGU1.component(1*3+2);
        const volScalarField GGU120 = GGU1.component(2*3+0);
        const volScalarField GGU121 = GGU1.component(2*3+1);
        const volScalarField GGU122 = GGU1.component(2*3+2);
        const volScalarField GGU200 = GGU2.component(0*3+0);
        const volScalarField GGU201 = GGU2.component(0*3+1);
        const volScalarField GGU202 = GGU2.component(0*3+2);
        const volScalarField GGU210 = GGU2.component(1*3+0);
        const volScalarField GGU211 = GGU2.component(1*3+1);
        const volScalarField GGU212 = GGU2.component(1*3+2);
        const volScalarField GGU220 = GGU2.component(2*3+0);
        const volScalarField GGU221 = GGU2.component(2*3+1);
        const volScalarField GGU222 = GGU2.component(2*3+2);
        return store
        (
            resultName_,
            sqr(GGU000+GGU011+GGU022)+sqr(GGU100+GGU111+GGU122)+sqr(GGU200+GGU211+GGU222)
        );
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::I2GGU::I2GGU
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U")
{
    setResultName(typeName, fieldName_);
}


// ************************************************************************* //
