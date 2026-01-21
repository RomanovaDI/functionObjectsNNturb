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
    Foam::functionObjects::I4

Group
    grpFieldFunctionObjects

Description
    Computes the fourth Pope's invariant.

    \f[
        I_4 = tr(r & r & s)
    \f]

    where
    \vartable
        I_4             | Fourth Pope's invariant    [1/s^3]
        r               | Rotation rate tensor    [1/s]
        s               | Strain rate tensor    [1/s]
        U               | Velocity     [m/s]
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
    I4
    {
        // Mandatory entries (unmodifiable)
        type        I4;
        libs        (fieldNNturbFunctionObjects);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                              | Type | Req'd | Dflt
      type         | Type name: I4                            | word |  yes  | -
      libs         | Library name: fieldFunctionObjectsNNturb | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func I4
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression

SourceFiles
    I4.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_I4_H
#define functionObjects_I4_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class I4 Declaration
\*---------------------------------------------------------------------------*/

class I4
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the I4 field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("I4");


    // Constructors

        //- Construct from Time and dictionary
        I4
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        I4(const I4&) = delete;

        //- No copy assignment
        void operator=(const I4&) = delete;


    //- Destructor
    virtual ~I4() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
