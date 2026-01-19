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
    Foam::functionObjects::SVDi

Group
    grpFieldFunctionObjects

Description
    Computes the c vector of coefficients of C_i (cel centers coordinates) as a pseudosolution of a linear system mag(U) = \sum\limits_{i=1}^{3} c_i C_i.

    \f[
        c_i = C_i^+ * mag(U)
    \f]

    where
    \vartable
        c_i     | Coefficients of C_i
        C_i     | Cell center coordinate [m]
        mag(U)  | Velocity magnitude [m/s]
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
    SVDi
    {
        // Mandatory entries (unmodifiable)
        type        SVDi;
        libs        (fieldNNturbFunctionObjects);

        // Mandatory (inherited) entry (runtime modifiable)
        fields      (<field1> <field2> ... <fieldN>);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                              | Type | Req'd | Dflt
      type         | Type name: SVDi                           | word |  yes  | -
      libs         | Library name: fieldFunctionObjectsNNturb | word |  yes  | -
      fields       | Names of the operand fields        | wordList | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func SVDi
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression

SourceFiles
    SVDi.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_SVDi_H
#define functionObjects_SVDi_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class SVDi Declaration
\*---------------------------------------------------------------------------*/

class SVDi
:
    public fvMeshFunctionObject
{
    // Private Data

        //- List of result fields
        PtrList<volScalarField> resultFields_;


    // Private Member Functions

        //- Calculate the SVDi fields
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("SVDi");


    // Constructors

        //- Construct from Time and dictionary
        SVDi
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        SVDi(const SVDi&) = delete;

        //- No copy assignment
        void operator=(const SVDi&) = delete;


    //- Destructor
    virtual ~SVDi() = default;


    // Member Functions

        //- Read the SVDi data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the SVDi fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
