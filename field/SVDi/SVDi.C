#include "SVDi.H"
#include "fvcGrad.H"
#include "SVD.H"
#include "RectangularMatrix.H"
#include "scalar.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"
#include "symmTransformField.H"  // Needed for symm tensor operations

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(SVDi, 0);
    addToRunTimeSelectionTable(functionObject, SVDi, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::SVDi::SVDi
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    resultFields_(5)  // Initialize with size 5
{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::SVDi::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    // Clear and resize to 5 fields
    resultFields_.clear();
    resultFields_.resize(5);

    // Create scalar fields for coefficients
    resultFields_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                "c1",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    
    resultFields_.set
    (
        1,
        new volScalarField
        (
            IOobject
            (
                "c2",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    
    resultFields_.set
    (
        2,
        new volScalarField
        (
            IOobject
            (
                "c3",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    
    // Create symmTensor fields for R and R_SVD
    resultFields_.set
    (
        3,
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor("zero", dimVelocity*dimVelocity, Zero)
        )
    );
    
    resultFields_.set
    (
        4,
        new volSymmTensorField
        (
            IOobject
            (
                "R_SVD",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor("zero", dimVelocity*dimVelocity, Zero)
        )
    );

    return true;
}


bool Foam::functionObjects::SVDi::execute()
{
    return calc();
}


bool Foam::functionObjects::SVDi::write()
{
    Log << type() << " " << name() << " write:" << endl;
    
    // Write each field
    forAll(resultFields_, i)
    {
        // Cast to appropriate type and write
        if (i < 3)  // scalar fields
        {
            const volScalarField& fld = 
                dynamicCast<const volScalarField&>(*resultFields_[i]);
            fld.write();
        }
        else  // symmTensor fields
        {
            const volSymmTensorField& fld = 
                dynamicCast<const volSymmTensorField&>(*resultFields_[i]);
            fld.write();
        }
    }
    
    return true;
}


bool Foam::functionObjects::SVDi::calc()
{
    Log << type() << ' ' << name() << " execute:" << nl;

    // Get required fields
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const volScalarField& nut = mesh_.lookupObject<volScalarField>("nut");

    // Calculate strain rate tensor and its powers
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU();
    
    volSymmTensorField s(symm(gradU));
    volSymmTensorField ss(s & s);
    volSymmTensorField sss(s & s & s);
    
    // Calculate R = 2 * nut * s
    volSymmTensorField R(2 * nut * s);
    
    // Get references to result fields with proper casting
    volScalarField& c1 = dynamicCast<volScalarField&>(*resultFields_[0]);
    volScalarField& c2 = dynamicCast<volScalarField&>(*resultFields_[1]);
    volScalarField& c3 = dynamicCast<volScalarField&>(*resultFields_[2]);
    volSymmTensorField& RField = dynamicCast<volSymmTensorField&>(*resultFields_[3]);
    volSymmTensorField& R_SVD = dynamicCast<volSymmTensorField&>(*resultFields_[4]);
    
    // Store original R
    RField = R;

    // Matrix dimensions: 6 symm tensor components × 3 basis tensors
    const label numComponents = 6;  // Symmetric tensor has 6 independent components
    const label numFields = 3;

    // Loop over all cells
    forAll(R.internalField(), cellI)
    {
        // Create the Ti matrix (6×3)
        scalarRectangularMatrix TiMatrix(numComponents, numFields, 0.0);

        // Fill the Ti matrix with tensor components from s, ss, sss
        // Symmetric tensor components: xx, xy, xz, yy, yz, zz
        for (label comp = 0; comp < numComponents; comp++)
        {
            TiMatrix(comp, 0) = s.internalField()[cellI].component(comp);
            TiMatrix(comp, 1) = ss.internalField()[cellI].component(comp);
            TiMatrix(comp, 2) = sss.internalField()[cellI].component(comp);
        }

        // SVD decomposition
        SVD svdA(TiMatrix);
        scalarRectangularMatrix TiMatrixPInv = svdA.VSinvUt();
        
        // Create R vector as column matrix (6×1)
        scalarRectangularMatrix Rm(numComponents, 1, 0.0);
        for (label comp = 0; comp < numComponents; comp++)
        {
            Rm(comp, 0) = R.internalField()[cellI].component(comp);
        }

        // Solve: SVDi = pinv(TiMatrix) * R
        scalarRectangularMatrix SVDi(numFields, 1, 0.0);
        
        // Matrix multiplication: SVDi = TiMatrixPInv * Rm
        for (label i = 0; i < numFields; i++)
        {
            for (label j = 0; j < numComponents; j++)
            {
                SVDi(i, 0) += TiMatrixPInv(i, j) * Rm(j, 0);
            }
        }

        // Store coefficients
        c1.internalField()[cellI] = SVDi(0, 0);
        c2.internalField()[cellI] = SVDi(1, 0);
        c3.internalField()[cellI] = SVDi(2, 0);
    }

    // Reconstruct R_SVD from coefficients and basis tensors
    forAll(R_SVD.internalField(), cellI)
    {
        symmTensor sumTensor = 
            c1.internalField()[cellI] * s.internalField()[cellI] +
            c2.internalField()[cellI] * ss.internalField()[cellI] +
            c3.internalField()[cellI] * sss.internalField()[cellI];
            
        R_SVD.internalField()[cellI] = sumTensor;
    }

    // Process boundary fields
    forAll(R_SVD.boundaryField(), patchI)
    {
        fvPatchSymmTensorField& r_svdPatch = R_SVD.boundaryFieldRef()[patchI];
        
        const fvPatchScalarField& c0Patch = c1.boundaryField()[patchI];
        const fvPatchScalarField& c1Patch = c2.boundaryField()[patchI];
        const fvPatchScalarField& c2Patch = c3.boundaryField()[patchI];
        
        const fvPatchSymmTensorField& sPatch = s.boundaryField()[patchI];
        const fvPatchSymmTensorField& ssPatch = ss.boundaryField()[patchI];
        const fvPatchSymmTensorField& sssPatch = sss.boundaryField()[patchI];
        
        forAll(r_svdPatch, faceI)
        {
            symmTensor sumTensor = 
                c0Patch[faceI] * sPatch[faceI] +
                c1Patch[faceI] * ssPatch[faceI] +
                c2Patch[faceI] * sssPatch[faceI];
            
            r_svdPatch[faceI] = sumTensor;
        }
    }

    // Update boundary conditions
    c1.correctBoundaryConditions();
    c2.correctBoundaryConditions();
    c3.correctBoundaryConditions();
    R_SVD.correctBoundaryConditions();

    Log << "    Calculated C1-C3 fields" << endl;

    return true;
}
