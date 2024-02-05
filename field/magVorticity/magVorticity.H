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
    Foam::functionObjects::magVorticity

Group
    grpFieldFunctionObjects

Description
    Computes the magVorticity, the magnitude of the curl of the velocity.

    \f[
        |\vec \omega| = |\nabla \times \vec U|
    \f]

    where
    \vartable
        |\vec \omega|   | Magnitude of vorticity vector    [1/s]
        \vec U          | Velocity vector     [m/s]
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
    magVorticity1
    {
        // Mandatory entries (unmodifiable)
        type        magVorticity;
        libs        (fieldFunctionObjectsNNturb);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: magVorticity               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func magVorticity
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::magVorticity

SourceFiles
    magVorticity.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_magVorticity_H
#define functionObjects_magVorticity_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class magVorticity Declaration
\*---------------------------------------------------------------------------*/

class magVorticity
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the magVorticity field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("magVorticity");


    // Constructors

        //- Construct from Time and dictionary
        magVorticity
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        magVorticity(const magVorticity&) = delete;

        //- No copy assignment
        void operator=(const magVorticity&) = delete;


    //- Destructor
    virtual ~magVorticity() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
