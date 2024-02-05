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
    Foam::functionObjects::strainTensor

Group
    grpFieldFunctionObjects

Description
    Computes the strain tensor.

    \f[
        s = \sqrt 0.5 (\nabla U + (\nabla U)^T)
    \f]

    where
    \vartable
        U        | Velocity    [m/s]
    \endvartable

    Operands:
    \table
      Operand        | Type           | Location
      input          | volVectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file    | -              | -
      output field   | volTensorField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    strainTensor1
    {
        // Mandatory entries (unmodifiable)
        type        strainTensor;
        libs        (fieldFunctionObjectsNNturb);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                              | Type | Req'd | Dflt
      type         | Type name: strainTensor                  | word |  yes  | -
      libs         | Library name: fieldFunctionObjectsCustom | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func strainTensor
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression

SourceFiles
    strainTensor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_strainTensor_H
#define functionObjects_strainTensor_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class strainTensor Declaration
\*---------------------------------------------------------------------------*/

class strainTensor
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the strainTensor field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("strainTensor");


    // Constructors

        //- Construct from Time and dictionary
        strainTensor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        strainTensor(const strainTensor&) = delete;

        //- No copy assignment
        void operator=(const strainTensor&) = delete;


    //- Destructor
    virtual ~strainTensor() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
