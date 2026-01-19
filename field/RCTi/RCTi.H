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
    Foam::functionObjects::RCTi

Group
    grpFieldFunctionObjects

Description
    Write Reynolds stress calculated via CTi and Ti.

    \f[
        RCTi_j = symm(\sum\limits_{i=1}^{10} C_i T_{i, j}).
    \f]

    where
    \vartable
        RCTi        | Reynolds stress tensor (6 components) [1/s]
        T_i         | Pope's tensors [1/s]
        CT_i        | Coefficients for Pope's tensors
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
        type        RCTi;
        libs        (fieldFunctionObjectsNNturb);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: RCTi                      | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func RCTi
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::RCTi

SourceFiles
    RCTi.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_RCTi_H
#define functionObjects_RCTi_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

class RCTi
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Result Reynolds stress tensor field (non-symmetric)
        volSymmTensorField RCTiField_;


    // Private Member Functions

        //- Calculate the Reynolds stress tensor
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("RCTi");


    // Constructors

        //- Construct from Time and dictionary
        RCTi
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        RCTi(const RCTi&) = delete;

        //- No copy assignment
        void operator=(const RCTi&) = delete;


    //- Destructor
    virtual ~RCTi() = default;


    // Member Functions

        //- Read the RCTi data
        virtual bool read(const dictionary&);

        //- Execute, calculate the field
        virtual bool execute();

        //- Write the RCTi field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
