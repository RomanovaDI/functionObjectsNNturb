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
    Foam::functionObjects::I0GGU

Group
    grpFieldFunctionObjects

Description
    Computes the quadratic invariant of gradient of velocity gradient 3rd rank tensor.

    \f[
        I0GGU = \sum\limits_{i, j, k} \frac{\partial^2 U_i}{\partial x_i \partial x_j} \frac{\partial^2 U_k}{\partial x_k \partial x_j}
    \f]

    where
    \vartable
        I0GGUU | Quadratic invariant of gradient of velocity gradient 3rd rank tensor    [1/m^2 s^2]
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
    I0GGU
    {
        // Mandatory entries (unmodifiable)
        type        I0GGU;
        libs        (fieldNNturbFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                              | Type | Req'd | Dflt
      type         | Type name: I0GGU                         | word |  yes  | -
      libs         | Library name: fieldFunctionObjectsNNturb | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func I0GGU
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression

SourceFiles
    I0GGU.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_I2gradU_H
#define functionObjects_I2gradU_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class I0GGU Declaration
\*---------------------------------------------------------------------------*/

class I0GGU
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the I0GGU field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("I0GGU");


    // Constructors

        //- Construct from Time and dictionary
        I0GGU
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        I0GGU(const I0GGU&) = delete;

        //- No copy assignment
        void operator=(const I0GGU&) = delete;


    //- Destructor
    virtual ~I0GGU() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
