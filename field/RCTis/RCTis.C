#include "RCTis.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(RCTis, 0);
    addToRunTimeSelectionTable(functionObject, RCTis, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::RCTis::RCTis
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    RCTiField_
    (
        IOobject
        (
            "RCTis",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("zero", dimVelocity*dimVelocity, Zero)
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::RCTis::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    return true;
}


bool Foam::functionObjects::RCTis::execute()
{
    return calc();
}


bool Foam::functionObjects::RCTis::write()
{
    Log << type() << " " << name() << " write:" << endl;
    
    // Write the Reynolds stress tensor field
    RCTiField_.write();
    
    return true;
}


bool Foam::functionObjects::RCTis::calc()
{
    Log << type() << ' ' << name() << " execute:" << nl;

    // Check if all required CTi scalar fields exist (CT1-CT10)
    for (label i = 1; i <= 10; ++i)
    {
        word fieldName = "CT" + Foam::name(i) + "s";
        if (!mesh_.foundObject<volScalarField>(fieldName))
        {
            WarningInFunction
                << "Required scalar field " << fieldName << " not found." << endl;
            return false;
        }
    }

    // Check if all required Ti tensor fields exist (T1-T10)
    for (label i = 1; i <= 10; ++i)
    {
        word fieldName = "T" + Foam::name(i) + "s";
        if (!mesh_.foundObject<volSymmTensorField>(fieldName))
        {
            WarningInFunction
                << "Required tensor field " << fieldName << " not found." << endl;
            return false;
        }
    }

    // Get references to all CTi scalar fields
    const volScalarField& CT1 = mesh_.lookupObject<volScalarField>("CT1s");
    const volScalarField& CT2 = mesh_.lookupObject<volScalarField>("CT2s");
    const volScalarField& CT3 = mesh_.lookupObject<volScalarField>("CT3s");
    const volScalarField& CT4 = mesh_.lookupObject<volScalarField>("CT4s");
    const volScalarField& CT5 = mesh_.lookupObject<volScalarField>("CT5s");
    const volScalarField& CT6 = mesh_.lookupObject<volScalarField>("CT6s");
    const volScalarField& CT7 = mesh_.lookupObject<volScalarField>("CT7s");
    const volScalarField& CT8 = mesh_.lookupObject<volScalarField>("CT8s");
    const volScalarField& CT9 = mesh_.lookupObject<volScalarField>("CT9s");
    const volScalarField& CT10 = mesh_.lookupObject<volScalarField>("CT10s");

    // Get references to all Ti tensor fields
    const volSymmTensorField& T1 = mesh_.lookupObject<volSymmTensorField>("T1s");
    const volSymmTensorField& T2 = mesh_.lookupObject<volSymmTensorField>("T2s");
    const volSymmTensorField& T3 = mesh_.lookupObject<volSymmTensorField>("T3s");
    const volSymmTensorField& T4 = mesh_.lookupObject<volSymmTensorField>("T4s");
    const volSymmTensorField& T5 = mesh_.lookupObject<volSymmTensorField>("T5s");
    const volSymmTensorField& T6 = mesh_.lookupObject<volSymmTensorField>("T6s");
    const volSymmTensorField& T7 = mesh_.lookupObject<volSymmTensorField>("T7s");
    const volSymmTensorField& T8 = mesh_.lookupObject<volSymmTensorField>("T8s");
    const volSymmTensorField& T9 = mesh_.lookupObject<volSymmTensorField>("T9s");
    const volSymmTensorField& T10 = mesh_.lookupObject<volSymmTensorField>("T10s");

    // Reset RCTis field to zero
    RCTiField_.primitiveFieldRef() = symmTensor::zero;

    // Calculate Reynolds stress tensor: RCTis = sum(CTi * Ti)
    // Using primitive field access for efficiency
    forAll(RCTiField_.primitiveField(), cellI)
    {
        // Initialize with first term
        symmTensor sumTensor = CT1.primitiveField()[cellI] * T1.primitiveField()[cellI];
        
        // Add the remaining terms
        sumTensor += CT2.primitiveField()[cellI] * T2.primitiveField()[cellI];
        sumTensor += CT3.primitiveField()[cellI] * T3.primitiveField()[cellI];
        sumTensor += CT4.primitiveField()[cellI] * T4.primitiveField()[cellI];
        sumTensor += CT5.primitiveField()[cellI] * T5.primitiveField()[cellI];
        sumTensor += CT6.primitiveField()[cellI] * T6.primitiveField()[cellI];
        sumTensor += CT7.primitiveField()[cellI] * T7.primitiveField()[cellI];
        sumTensor += CT8.primitiveField()[cellI] * T8.primitiveField()[cellI];
        sumTensor += CT9.primitiveField()[cellI] * T9.primitiveField()[cellI];
        sumTensor += CT10.primitiveField()[cellI] * T10.primitiveField()[cellI];
        
        // Store as full tensor (not symmetric)
        RCTiField_.primitiveFieldRef()[cellI] = sumTensor;
    }

    // Handle boundary fields
    forAll(RCTiField_.boundaryField(), patchI)
    {
        // Get references to patch fields
        fvPatchSymmTensorField& rctiPatch = RCTiField_.boundaryFieldRef()[patchI];
        
        const fvPatchScalarField& ct1Patch = CT1.boundaryField()[patchI];
        const fvPatchScalarField& ct2Patch = CT2.boundaryField()[patchI];
        const fvPatchScalarField& ct3Patch = CT3.boundaryField()[patchI];
        const fvPatchScalarField& ct4Patch = CT4.boundaryField()[patchI];
        const fvPatchScalarField& ct5Patch = CT5.boundaryField()[patchI];
        const fvPatchScalarField& ct6Patch = CT6.boundaryField()[patchI];
        const fvPatchScalarField& ct7Patch = CT7.boundaryField()[patchI];
        const fvPatchScalarField& ct8Patch = CT8.boundaryField()[patchI];
        const fvPatchScalarField& ct9Patch = CT9.boundaryField()[patchI];
        const fvPatchScalarField& ct10Patch = CT10.boundaryField()[patchI];
        
        const fvPatchSymmTensorField& t1Patch = T1.boundaryField()[patchI];
        const fvPatchSymmTensorField& t2Patch = T2.boundaryField()[patchI];
        const fvPatchSymmTensorField& t3Patch = T3.boundaryField()[patchI];
        const fvPatchSymmTensorField& t4Patch = T4.boundaryField()[patchI];
        const fvPatchSymmTensorField& t5Patch = T5.boundaryField()[patchI];
        const fvPatchSymmTensorField& t6Patch = T6.boundaryField()[patchI];
        const fvPatchSymmTensorField& t7Patch = T7.boundaryField()[patchI];
        const fvPatchSymmTensorField& t8Patch = T8.boundaryField()[patchI];
        const fvPatchSymmTensorField& t9Patch = T9.boundaryField()[patchI];
        const fvPatchSymmTensorField& t10Patch = T10.boundaryField()[patchI];
        
        // Calculate for each face in the patch
        forAll(rctiPatch, faceI)
        {
            symmTensor sumTensor = ct1Patch[faceI] * t1Patch[faceI];
            sumTensor += ct2Patch[faceI] * t2Patch[faceI];
            sumTensor += ct3Patch[faceI] * t3Patch[faceI];
            sumTensor += ct4Patch[faceI] * t4Patch[faceI];
            sumTensor += ct5Patch[faceI] * t5Patch[faceI];
            sumTensor += ct6Patch[faceI] * t6Patch[faceI];
            sumTensor += ct7Patch[faceI] * t7Patch[faceI];
            sumTensor += ct8Patch[faceI] * t8Patch[faceI];
            sumTensor += ct9Patch[faceI] * t9Patch[faceI];
            sumTensor += ct10Patch[faceI] * t10Patch[faceI];
            
            rctiPatch[faceI] = sumTensor;
        }
    }

    // Update boundary conditions if needed
    RCTiField_.correctBoundaryConditions();

    Log << "    Calculated Reynolds stress tensor RCTis (symmetric)" << endl;

    return true;
}


// ************************************************************************* //
