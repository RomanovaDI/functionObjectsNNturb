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
    Foam::functionObjects::T6

Group
    grpFieldFunctionObjects

Description
    Computes the T6 Pope's tensor.

    \f[
        T6 = r & r & s + s & r & r - I * 2 * tr(s & r & r) / 3
    \f]

    where
    \vartable
        T6       | Sixth Pope's tensor [1/s]
        s        | Strain tensor       [1/s]
        r        | Rotation tensor     [1/s]
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
    strainRate1
    {
        // Mandatory entries (unmodifiable)
        type        T6;
        libs        (fieldFunctionObjectsNNturb);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: T6                      | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func T6
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::T6

SourceFiles
    T6.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_T6_H
#define functionObjects_T6_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class T6 Declaration
\*---------------------------------------------------------------------------*/

class T6
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the T6 field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("T6");


    // Constructors

        //- Construct from Time and dictionary
        T6
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        T6(const T6&) = delete;

        //- No copy assignment
        void operator=(const T6&) = delete;


    //- Destructor
    virtual ~T6() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
