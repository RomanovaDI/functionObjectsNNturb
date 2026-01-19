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
    Foam::functionObjects::T5

Group
    grpFieldFunctionObjects

Description
    Computes the T5 Pope's tensor.

    \f[
        T5 = r & s & s - s & s & r
    \f]

    where
    \vartable
        T5       | Fifth Pope's tensor [1/s]
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
        type        T5;
        libs        (fieldFunctionObjectsNNturb);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: T5                      | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func T5
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::T5

SourceFiles
    T5.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_T5_H
#define functionObjects_T5_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class T5 Declaration
\*---------------------------------------------------------------------------*/

class T5
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the T5 field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("T5");


    // Constructors

        //- Construct from Time and dictionary
        T5
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        T5(const T5&) = delete;

        //- No copy assignment
        void operator=(const T5&) = delete;


    //- Destructor
    virtual ~T5() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
