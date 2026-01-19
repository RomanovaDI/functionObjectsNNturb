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

#include "CTi.H"
#include "fvcGrad.H"
#include "SVD.H"
#include "RectangularMatrix.H"
#include "scalar.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CTi, 0);
    addToRunTimeSelectionTable(functionObject, CTi, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::CTi::CTi
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    resultFields_()
{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CTi::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    // Create result fields CT1-CT10
    resultFields_.clear();
    resultFields_.setSize(10);

    for (label i = 0; i < 10; ++i)
    {
        word fieldName = "CT" + Foam::name(i+1);
        
        resultFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    return true;
}


bool Foam::functionObjects::CTi::execute()
{
    return calc();
}


bool Foam::functionObjects::CTi::write()
{
    Log << type() << " " << name() << " write:" << endl;
    
    // Write all result fields
    forAll(resultFields_, i)
    {
        resultFields_[i].write();
    }
    
    return true;
}


bool Foam::functionObjects::CTi::calc()
{
    Log << type() << ' ' << name() << " execute:" << nl;

    // Check if required tensor fields T1-T10 exist
/*    for (label i = 1; i <= 10; ++i)
    {
        word fieldName = "T" + Foam::name(i);
        if (!mesh_.foundObject<volTensorField>(fieldName))
        {
            WarningInFunction
                << "Required tensor field " << fieldName << " not found." << endl;
            return false;
        }
    }

    // Check if reference tensor field R exists
    if (!mesh_.foundObject<volTensorField>("RT"))
    {
        WarningInFunction
            << "Reference tensor field RT not found." << endl;
        return false;
    }*/

    // Get references to the input tensor fields
    const volTensorField& T1 = mesh_.lookupObject<volTensorField>("T1");
    const volTensorField& T2 = mesh_.lookupObject<volTensorField>("T2");
    const volTensorField& T3 = mesh_.lookupObject<volTensorField>("T3");
    const volTensorField& T4 = mesh_.lookupObject<volTensorField>("T4");
    const volTensorField& T5 = mesh_.lookupObject<volTensorField>("T5");
    const volTensorField& T6 = mesh_.lookupObject<volTensorField>("T6");
    const volTensorField& T7 = mesh_.lookupObject<volTensorField>("T7");
    const volTensorField& T8 = mesh_.lookupObject<volTensorField>("T8");
    const volTensorField& T9 = mesh_.lookupObject<volTensorField>("T9");
    const volTensorField& T10 = mesh_.lookupObject<volTensorField>("T10");

    // Get reference to the reference tensor field R (predefined name)
    const volSymmTensorField& R = mesh_.lookupObject<volSymmTensorField>("R");

    volScalarField magR(mag(R));
    scalar RmaxMag = gMax(magR.internalField());
    volScalarField magT1(mag(T1));
    scalar T1maxMag = gMax(magT1.internalField());
    volScalarField magT2(mag(T2));
    scalar T2maxMag = gMax(magT2.internalField());
    volScalarField magT3(mag(T3));
    scalar T3maxMag = gMax(magT3.internalField());
    volScalarField magT4(mag(T4));
    scalar T4maxMag = gMax(magT4.internalField());
    volScalarField magT5(mag(T5));
    scalar T5maxMag = gMax(magT5.internalField());
    volScalarField magT6(mag(T6));
    scalar T6maxMag = gMax(magT6.internalField());
    volScalarField magT7(mag(T7));
    scalar T7maxMag = gMax(magT7.internalField());
    volScalarField magT8(mag(T8));
    scalar T8maxMag = gMax(magT8.internalField());
    volScalarField magT9(mag(T9));
    scalar T9maxMag = gMax(magT9.internalField());
    volScalarField magT10(mag(T10));
    scalar T10maxMag = gMax(magT10.internalField());
    scalarList maxMags = {T1maxMag, T2maxMag, T3maxMag, T4maxMag, T5maxMag, T6maxMag, T7maxMag, T8maxMag, T9maxMag, T10maxMag, RmaxMag};
    scalar maxMag = gMax(maxMags);

    // Initialize fields to zero
    forAll(resultFields_, i)
    {
        resultFields_[i].primitiveFieldRef() = 0.0;
    }

    // Matrix dimensions: 9 tensor components × 10 tensor fields
    const label numComponents = 9;
    const label numFields = 10;

    // Loop over all cells
    forAll(R.internalField(), cellI)
    {
        // Create the Ti matrix (9×10)
        scalarRectangularMatrix TiMatrix(numComponents, numFields);

        // Fill the Ti matrix with tensor components from T1-T10
        // Column 0: T1 components, Column 1: T2 components, etc.
        //for (label comp = 0; comp < numComponents; ++comp)
        for (label comp = 0; comp < numComponents; comp++)
        {
            TiMatrix(comp, 0) = T1.primitiveField()[cellI].component(comp)/maxMag;
            TiMatrix(comp, 1) = T2.primitiveField()[cellI].component(comp)/maxMag;
            TiMatrix(comp, 2) = T3.primitiveField()[cellI].component(comp)/maxMag;
            TiMatrix(comp, 3) = T4.primitiveField()[cellI].component(comp)/maxMag;
            TiMatrix(comp, 4) = T5.primitiveField()[cellI].component(comp)/maxMag;
            TiMatrix(comp, 5) = T6.primitiveField()[cellI].component(comp)/maxMag;
            TiMatrix(comp, 6) = T7.primitiveField()[cellI].component(comp)/maxMag;
            TiMatrix(comp, 7) = T8.primitiveField()[cellI].component(comp)/maxMag;
            TiMatrix(comp, 8) = T9.primitiveField()[cellI].component(comp)/maxMag;
            TiMatrix(comp, 9) = T10.primitiveField()[cellI].component(comp)/maxMag;
        }

        // SVD decomposition: TiMatrix = U * D * V^T
        SVD svdA(TiMatrix);
        scalarRectangularMatrix TiMatrixPInv = svdA.VSinvUt();
        
        // Create R vector as column matrix (9×1) from tensor R components
        scalarRectangularMatrix Rm(numComponents, 1);
        //for (label comp = 0; comp < numComponents; ++comp)
        /*for (label comp = 0; comp < numComponents; comp++)
        {
            Rm(comp, 0) = R.primitiveField()[cellI].component(comp);
        }*/
        Rm(0, 0) = R.primitiveField()[cellI].component(0)/maxMag;
        Rm(1, 0) = R.primitiveField()[cellI].component(1)/maxMag;
        Rm(2, 0) = R.primitiveField()[cellI].component(2)/maxMag;
        Rm(3, 0) = R.primitiveField()[cellI].component(1)/maxMag;
        Rm(4, 0) = R.primitiveField()[cellI].component(3)/maxMag;
        Rm(5, 0) = R.primitiveField()[cellI].component(4)/maxMag;
        Rm(6, 0) = R.primitiveField()[cellI].component(2)/maxMag;
        Rm(7, 0) = R.primitiveField()[cellI].component(4)/maxMag;
        Rm(8, 0) = R.primitiveField()[cellI].component(5)/maxMag;

        // Solve: CTi = pinv(TiMatrix) * R
        scalarRectangularMatrix CTi(numFields, 1);
        multiply(CTi, TiMatrixPInv, Rm);

        // Assign results to output fields
        for (label field = 0; field < numComponents; field++)
            resultFields_[field].primitiveFieldRef()[cellI] = CTi(field, 0);  // CT_field
        /*resultFields_[1].primitiveFieldRef()[cellI] = CTi(1, 0);  // CT2
        resultFields_[2].primitiveFieldRef()[cellI] = CTi(2, 0);  // CT3
        resultFields_[3].primitiveFieldRef()[cellI] = CTi(3, 0);  // CT4
        resultFields_[4].primitiveFieldRef()[cellI] = CTi(4, 0);  // CT5
        resultFields_[5].primitiveFieldRef()[cellI] = CTi(5, 0);  // CT6
        resultFields_[6].primitiveFieldRef()[cellI] = CTi(6, 0);  // CT7
        resultFields_[7].primitiveFieldRef()[cellI] = CTi(7, 0);  // CT8
        resultFields_[8].primitiveFieldRef()[cellI] = CTi(8, 0);  // CT9
        resultFields_[9].primitiveFieldRef()[cellI] = CTi(9, 0);  // CT10*/
    }

    Log << "    Calculated CT1-CT10 fields" << endl;

    return true;
}


// ************************************************************************* //

