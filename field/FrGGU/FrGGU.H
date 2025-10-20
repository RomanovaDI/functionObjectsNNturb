/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::FrGGU

Group
    grpFieldFunctionObjects

Description
    Computes the quadratic invariant of gradient of velocity gradient 3rd rank tensor.

    \f[
        FrGGU = \sum\limits_{i, j, k} (\frac{\partial^2 U_i}{\partial x_j \partial x_k})^2
    \f]

    where
    \vartable
        FrGGU  | Frobenius norm of gradient of velocity gradient 3rd rank tensor    [1/m s]
        U      | Velocity     [m/s]
    \endvartable

    Operands:
    \table
      Operand        | Type           | Location
      input          | volVectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -              | -
      output field   | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    FrGGU
    {
        // Mandatory entries (unmodifiable)
        type        FrGGU;
        libs        (fieldNNturbFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                              | Type | Req'd | Dflt
      type         | Type name: FrGGU                         | word |  yes  | -
      libs         | Library name: fieldFunctionObjectsNNturb | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func FrGGU
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression

SourceFiles
    FrGGU.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_FrGGU_H
#define functionObjects_FrGGU_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class FrGGU Declaration
\*---------------------------------------------------------------------------*/

class FrGGU
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the FrGGU field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("FrGGU");


    // Constructors

        //- Construct from Time and dictionary
        FrGGU
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        FrGGU(const FrGGU&) = delete;

        //- No copy assignment
        void operator=(const FrGGU&) = delete;


    //- Destructor
    virtual ~FrGGU() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
