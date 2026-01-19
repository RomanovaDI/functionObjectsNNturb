/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

Description
    See reynoldsStressSVD.H

\*---------------------------------------------------------------------------*/

#include "reynoldsStressSVD.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "SVD.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reynoldsStressSVD, 0);
    
    addToRunTimeSelectionTable
    (
        functionObject,
        reynoldsStressSVD,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::reynoldsStressSVD::calcStrainRate()
{
    // Get velocity field
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    
    // Calculate velocity gradient
    volTensorField gradU = fvc::grad(U);
    
    // Calculate strain rate tensor: s = 0.5*(grad(U) + grad(U).T())
    s_ = 0.5 * (gradU + gradU.T());
    
    // Calculate square: ss = s & s
    ss_ = s_ & s_;
    
    // Calculate cube: sss = s & s & s
    sss_ = s_ & s_ & s_;
    
    if (writeTensors_ && log)
    {
        Info<< type() << " " << name() << ": "
            << "Calculated strain rate tensor and its powers" << endl;
        
        // Print some statistics
        dimensionedTensor sMin = min(s_);
        dimensionedTensor sMax = max(s_);
        Info<< "  s range (xx component): " << sMin.value().xx() 
            << " to " << sMax.value().xx() << endl;
    }
}


void Foam::functionObjects::reynoldsStressSVD::calcReynoldsStress()
{
    // Get turbulent viscosity
    const volScalarField& nut = mesh_.lookupObject<volScalarField>(nutName_);
    
    // Get turbulence kinetic energy if available
    const volScalarField* kPtr = mesh_.findObject<volScalarField>(kName_);
    
    if (writeTensors_ && log)
    {
        Info<< type() << " " << name() << ": "
            << "Calculating Reynolds stress tensor" << endl;
    }
    
    // Calculate Reynolds stress using Boussinesq hypothesis:
    // R = 2/3*k*I - 2*nu_t*s
    // where s is strain rate tensor
    
    forAll(mesh_.C(), cellI)
    {
        // Get turbulent viscosity and strain rate tensor for this cell
        const scalar nutCell = nut[cellI];
        const tensor& sCell = s_[cellI];
        
        // Calculate -2*nu_t*s
        tensor RCell = -2.0 * nutCell * sCell;
        
        // Add turbulent kinetic energy contribution if available
        if (kPtr)
        {
            const scalar kCell = (*kPtr)[cellI];
            // Add 2/3*k*I using tensor::I
            RCell += (2.0/3.0) * kCell * tensor::I;
        }
        
        R_[cellI] = RCell;
    }
    
    if (log)
    {
        dimensionedTensor Rmin = min(R_);
        dimensionedTensor Rmax = max(R_);
        Info<< type() << " " << name() << ": "
            << "Reynolds stress tensor range (xx component): "
            << Rmin.value().xx() << " to " 
            << Rmax.value().xx() << endl;
    }
}


void Foam::functionObjects::reynoldsStressSVD::performSVD()
{
    if (writeTensors_ && log)
    {
        Info<< type() << " " << name() << ": "
            << "Performing SVD decomposition in each cell" << endl;
    }
    
    // Counters for statistics
    label nCells = mesh_.nCells();
    label nFailedSVD = 0;
    scalar totalError = 0.0;
    
    // Loop over all cells
    forAll(mesh_.C(), cellI)
    {
        // Get tensor components for this cell
        const tensor& s = s_[cellI];
        const tensor& ss = ss_[cellI];
        const tensor& sss = sss_[cellI];
        const tensor& R = R_[cellI];
        
        // Create matrix A (9x3) and vector b (9) for linear system
        // R_ij = c_s * s_ij + c_ss * ss_ij + c_sss * sss_ij
        // We have 9 equations for 3 unknowns
        
        // Create RectangularMatrix for A (9 rows, 3 columns)
        RectangularMatrix<scalar> A(9, 3);
        // Create vector for b (9 elements)
        List<scalar> b(9);
        
        // Fill matrix A and vector b using tensor components
        // Order: xx, xy, xz, yx, yy, yz, zx, zy, zz
        label row = 0;
        for (label i = 0; i < 3; ++i)
        {
            for (label j = 0; j < 3; ++j)
            {
                A(row, 0) = s(i, j);
                A(row, 1) = ss(i, j);
                A(row, 2) = sss(i, j);
                b[row] = R(i, j);
                row++;
            }
        }
        
        try
        {
            // Perform SVD decomposition using OpenFOAM's SVD class
            SVD svd(A, svdTolerance_);
            
            // Check if SVD was successful by checking singular values
            // SVD returns a scalarDiagonalMatrix for singular values
            const scalarDiagonalMatrix& S = svd.S();
            
            // Count non-zero singular values
            label nNonZero = 0;
            for (label i = 0; i < S.size(); ++i)
            {
                if (S[i] > svdTolerance_)
                {
                    nNonZero++;
                }
            }
            
            if (nNonZero < 3)
            {
                nFailedSVD++;
                c_s_[cellI] = 0.0;
                c_ss_[cellI] = 0.0;
                c_sss_[cellI] = 0.0;
                continue;
            }
            
            // Solve using the SVD decomposition
            // Get U, V matrices
            const RectangularMatrix<scalar>& U = svd.U();
            const RectangularMatrix<scalar>& V = svd.V();
            
            // Compute solution using SVD: x = V * inv(S) * U.T() * b
            // 1. Compute U.T() * b
            List<scalar> UTb(3, 0.0);
            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 9; ++j)
                {
                    UTb[i] += U(j, i) * b[j];
                }
            }
            
            // 2. Compute inv(S) * (U.T() * b)
            List<scalar> SinvUTb(3, 0.0);
            for (label i = 0; i < 3; ++i)
            {
                if (S[i] > svdTolerance_)
                {
                    SinvUTb[i] = UTb[i] / S[i];
                }
                else
                {
                    // Singular value is too small, set to zero
                    SinvUTb[i] = 0.0;
                }
            }
            
            // 3. Compute V * (inv(S) * U.T() * b)
            List<scalar> x(3, 0.0);
            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    x[i] += V(i, j) * SinvUTb[j];
                }
            }
            
            // Store coefficients
            c_s_[cellI] = x[0];
            c_ss_[cellI] = x[1];
            c_sss_[cellI] = x[2];
            
            // Calculate reconstruction error for statistics
            tensor R_reconstructed = 
                c_s_[cellI] * s + 
                c_ss_[cellI] * ss + 
                c_sss_[cellI] * sss;
            
            scalar error = mag(R_reconstructed - R);
            totalError += error;
            
        }
        catch (const Foam::error& err)
        {
            nFailedSVD++;
            c_s_[cellI] = 0.0;
            c_ss_[cellI] = 0.0;
            c_sss_[cellI] = 0.0;
        }
    }
    
    // Write coefficient statistics
    if (log)
    {
        scalar avgError = 0.0;
        if (nCells - nFailedSVD > 0)
        {
            avgError = totalError / (nCells - nFailedSVD);
        }
        
        Info<< type() << " " << name() << ": "
            << "SVD Statistics:" << nl
            << "  Total cells: " << nCells << nl
            << "  Failed SVD: " << nFailedSVD 
            << " (" << (100.0 * nFailedSVD / max(nCells, 1)) << "%)" << nl
            << "  Average reconstruction error: " << avgError << nl
            << "  Coefficient statistics:" << nl
            << "    c_s:   min = " << min(c_s_).value() 
            << ", max = " << max(c_s_).value()
            << ", avg = " << average(c_s_).value() << nl
            << "    c_ss:  min = " << min(c_ss_).value() 
            << ", max = " << max(c_ss_).value()
            << ", avg = " << average(c_ss_).value() << nl
            << "    c_sss: min = " << min(c_sss_).value() 
            << ", max = " << max(c_sss_).value()
            << ", avg = " << average(c_sss_).value() << endl;
    }
}


void Foam::functionObjects::reynoldsStressSVD::reconstructRSVD()
{
    if (writeTensors_ && log)
    {
        Info<< type() << " " << name() << ": "
            << "Reconstructing R_SVD from coefficients" << endl;
    }
    
    // Reconstruct R_SVD = c_s * s + c_ss * ss + c_sss * sss
    forAll(mesh_.C(), cellI)
    {
        const scalar cs = c_s_[cellI];
        const scalar css = c_ss_[cellI];
        const scalar csss = c_sss_[cellI];
        
        const tensor& s = s_[cellI];
        const tensor& ss = ss_[cellI];
        const tensor& sss = sss_[cellI];
        
        R_SVD_[cellI] = cs * s + css * ss + csss * sss;
    }
    
    // Calculate and display reconstruction error
    if (log)
    {
        scalarField errorField(mag(R_SVD_ - R_));
        scalar maxError = gMax(errorField);
        scalar avgError = gAverage(errorField);
        
        Info<< type() << " " << name() << ": "
            << "Final reconstruction error: max = " << maxError
            << ", avg = " << avgError << endl;
        
        // Print sample comparison for first cell
        if (mesh_.nCells() > 0)
        {
            Info<< "Sample comparison (cell 0):" << nl
                << "  Original R: " << R_[0] << nl
                << "  Reconstructed R_SVD: " << R_SVD_[0] << nl
                << "  Difference: " << R_SVD_[0] - R_[0] << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reynoldsStressSVD::reynoldsStressSVD
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    UName_("U"),
    nutName_("nut"),
    kName_("k"),
    writeTensors_(true),
    svdTolerance_(1e-6),
    s_
    (
        IOobject
        (
            "s",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimVelocity/dimLength, Zero)
    ),
    ss_
    (
        IOobject
        (
            "ss",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(sqr(dimVelocity/dimLength), Zero)
    ),
    sss_
    (
        IOobject
        (
            "sss",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(pow(dimVelocity/dimLength, 3), Zero)
    ),
    R_
    (
        IOobject
        (
            "R",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(sqr(dimVelocity), Zero)
    ),
    c_s_
    (
        IOobject
        (
            "c_s",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength/dimTime, Zero)
    ),
    c_ss_
    (
        IOobject
        (
            "c_ss",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(sqr(dimLength/dimTime), Zero)
    ),
    c_sss_
    (
        IOobject
        (
            "c_sss",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(pow(dimLength/dimTime, 3), Zero)
    ),
    R_SVD_
    (
        IOobject
        (
            "R_SVD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(sqr(dimVelocity), Zero)
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::reynoldsStressSVD::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    
    // Read field names
    dict.readIfPresent("U", UName_);
    dict.readIfPresent("nut", nutName_);
    dict.readIfPresent("k", kName_);
    dict.readIfPresent("writeTensors", writeTensors_);
    dict.readIfPresent("svdTolerance", svdTolerance_);
    
    // Check if fields exist
    if (!mesh_.foundObject<volVectorField>(UName_))
    {
        FatalErrorInFunction
            << "Velocity field " << UName_ << " not found in mesh"
            << exit(FatalError);
    }
    
    if (!mesh_.foundObject<volScalarField>(nutName_))
    {
        FatalErrorInFunction
            << "Turbulent viscosity field " << nutName_ << " not found in mesh"
            << exit(FatalError);
    }
    
    // Check if k field exists (optional)
    if (!mesh_.foundObject<volScalarField>(kName_))
    {
        if (log)
        {
            Info<< type() << " " << name() << ": "
                << "Turbulence kinetic energy field " << kName_ 
                << " not found. Reynolds stress will be calculated without k term."
                << endl;
        }
    }
    
    if (log)
    {
        Info<< type() << " " << name() << ":" << nl
            << "    U: " << UName_ << nl
            << "    nut: " << nutName_ << nl
            << "    k: " << kName_ << nl
            << "    writeTensors: " << writeTensors_ << nl
            << "    SVD tolerance: " << svdTolerance_ << endl;
    }
    
    return true;
}


bool Foam::functionObjects::reynoldsStressSVD::execute()
{
    if (log)
    {
        Info<< type() << " " << name() << " executing..." << endl;
    }
    
    // Calculate strain rate tensor and its powers
    calcStrainRate();
    
    // Calculate Reynolds stress tensor
    calcReynoldsStress();
    
    // Perform SVD decomposition to find coefficients
    performSVD();
    
    // Reconstruct R_SVD from coefficients
    reconstructRSVD();
    
    return true;
}


bool Foam::functionObjects::reynoldsStressSVD::write()
{
    if (log)
    {
        Info<< type() << " " << name() << " writing fields:" << endl;
    }
    
    // Always write coefficients and reconstructed field
    c_s_.write();
    c_ss_.write();
    c_sss_.write();
    R_SVD_.write();
    
    // Write intermediate tensors if requested
    if (writeTensors_)
    {
        s_.write();
        ss_.write();
        sss_.write();
        R_.write();
    }
    
    return true;
}


// ************************************************************************* //
