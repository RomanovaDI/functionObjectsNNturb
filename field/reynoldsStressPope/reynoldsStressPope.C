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
    See reynoldsStressPope.H

\*---------------------------------------------------------------------------*/

#include "reynoldsStressPope.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reynoldsStressPope, 0);
    
    addToRunTimeSelectionTable
    (
        functionObject,
        reynoldsStressPope,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::reynoldsStressPope::calcStrainAndRotation()
{
    // Get velocity field
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    
    // Calculate velocity gradient
    volTensorField gradU = fvc::grad(U);
    
    // Calculate strain rate tensor: s = 0.5*(grad(U) + grad(U).T())
    s_ = 0.5 * (gradU + gradU.T());
    
    // Calculate rotation rate tensor: r = 0.5*(grad(U) - grad(U).T())
    r_ = 0.5 * (gradU - gradU.T());
    
    if (writeTensors_ && log)
    {
        Info<< type() << " " << name() << ": "
            << "Calculated strain rate and rotation rate tensors" << endl;
    }
}


void Foam::functionObjects::reynoldsStressPope::calcPopeTensors()
{
    // Create identity tensor
    tensor I = tensor::I;
    
    // Define dimensions for tensor operations
    dimensionSet strainDim = dimVelocity/dimLength;  // [0 0 -1 0 0 0 0]
    dimensionSet strainSqDim = sqr(strainDim);       // [0 0 -2 0 0 0 0]
    dimensionSet strainCubeDim = strainDim * strainSqDim;  // [0 0 -3 0 0 0 0]
    dimensionSet strainQuartDim = strainSqDim * strainSqDim; // [0 0 -4 0 0 0 0]
    dimensionSet strainQuintDim = strainCubeDim * strainSqDim; // [0 0 -5 0 0 0 0]
    
    // Temporary fields
    volTensorField ss(s_ & s_);
    volTensorField rr(r_ & r_);
    volTensorField sr(s_ & r_);
    volTensorField rs(r_ & s_);
    volTensorField srr(s_ & rr);
    volTensorField ssrr(ss & rr);
    
    // T1 = s
    T1_ = s_;
    
    // T2 = s & r - r & s
    T2_ = sr - rs;
    T2_.dimensions().reset(strainSqDim);
    
    // T3 = s & s - I * tr(s & s) / 3
    T3_ = ss - (tr(ss) / 3 * I);
    T3_.dimensions().reset(strainSqDim);
    
    // T4 = r & r - I * tr(r & r) / 3
    T4_ = rr - (I * tr(rr) / 3);
    T4_.dimensions().reset(strainSqDim);
    
    // T5 = r & s & s - s & s & r
    T5_ = (r_ & ss) - (ss & r_);
    T5_.dimensions().reset(strainCubeDim);
    
    // T6 = r & r & s + s & r & r - I * 2 * tr(s & r & r) / 3
    T6_ = (rr & s_) + srr - (I * 2 * tr(srr) / 3);
    T6_.dimensions().reset(strainCubeDim);
    
    // T7 = r & s & r & r - r & r & s & r
    T7_ = (rs & rr) - (rr & sr);
    T7_.dimensions().reset(strainQuartDim);
    
    // T8 = s & r & s & s - s & s & r & s
    T8_ = (sr & ss) - (ss & rs);
    T8_.dimensions().reset(strainQuartDim);
    
    // T9 = r & r & s & s + s & s & r & r - I * 2 * tr(s & s & r & r) / 3
    T9_ = (rr & ss) + ssrr - (I * 2 * tr(ssrr) / 3);
    T9_.dimensions().reset(strainQuartDim);
    
    // T10 = r & s & s & r & r - r & r & s & s & r
    T10_ = (r_ & ss & rr) - (rr & ss & r_);
    T10_.dimensions().reset(strainQuintDim);

    if (writeTensors_ && log)
    {
        Info<< type() << " " << name() << ": "
            << "Calculated Pope's 10 tensors T1-T10" << endl;
    }
}


void Foam::functionObjects::reynoldsStressPope::performSVD()
{
    if (!Rptr_)
    {
        FatalErrorInFunction
            << "Reynolds stress tensor field " << RName_ << " not found in mesh"
            << exit(FatalError);
    }
    
    const volSymmTensorField& R = *Rptr_;
    
    if (writeTensors_ && log)
    {
        Info<< type() << " " << name() << ": "
            << "Performing SVD decomposition in each cell to find Pope coefficients" << endl;
        
        if (useSymmetricEqns_)
        {
            Info<< "  Using 6 independent equations for symmetric tensor" << endl;
        }
        else
        {
            Info<< "  Using 9 equations (with symmetric components repeated)" << endl;
        }
    }
    
    // Counters for statistics
    label nCells = mesh_.nCells();
    label nRankDeficient = 0;
    scalar totalError = 0.0;
    
    // Get references to all Pope tensor fields in an array for easier access
    const volTensorField* T_fields[10] = {&T1_, &T2_, &T3_, &T4_, &T5_, 
                                           &T6_, &T7_, &T8_, &T9_, &T10_};
    
    // Determine matrix dimensions based on symmetry
    label nRows = useSymmetricEqns_ ? 6 : 9;  // 6 independent components for symmTensor
    label nCols = 10;  // 10 Pope coefficients
    
    // Loop over all cells
    forAll(mesh_.C(), cellI)
    {
        // Get Reynolds stress tensor for this cell
        const symmTensor& Rcell = R[cellI];
        
        // Create matrix A and vector b for linear system
        // R_ij = sum_{k=1}^{10} cTk * Tk_ij
        
        // Create RectangularMatrix for A
        RectangularMatrix<scalar> A(nRows, nCols);
        // Create vector for b
        List<scalar> b(nRows);
        
        // Fill matrix A and vector b
        if (useSymmetricEqns_)
        {
            // Use only 6 independent equations for symmetric tensor
            // Order: xx, xy, xz, yy, yz, zz
            
            for (label eqn = 0; eqn < 6; ++eqn)
            {
                // Get tensor indices for this equation
                label i, j;
                switch(eqn)
                {
                    case 0: i = 0; j = 0; break; // xx
                    case 1: i = 0; j = 1; break; // xy
                    case 2: i = 0; j = 2; break; // xz
                    case 3: i = 1; j = 1; break; // yy
                    case 4: i = 1; j = 2; break; // yz
                    case 5: i = 2; j = 2; break; // zz
                    default: i = 0; j = 0; break;
                }
                
                // Fill Pope tensor components for this equation
                for (label t = 0; t < 10; ++t)
                {
                    A(eqn, t) = (*T_fields[t])[cellI](i, j);
                }
                
                // Fill right-hand side with Reynolds stress component
                b[eqn] = Rcell.component(eqn);
            }
        }
        else
        {
            // Use 9 equations (with symmetric components repeated)
            // Order: xx, xy, xz, yx, yy, yz, zx, zy, zz
            
            label row = 0;
            for (label i = 0; i < 3; ++i)
            {
                for (label j = 0; j < 3; ++j)
                {
                    // Fill Pope tensor components for this equation
                    for (label t = 0; t < 10; ++t)
                    {
                        A(row, t) = (*T_fields[t])[cellI](i, j);
                    }
                    
                    // Fill right-hand side with Reynolds stress component
                    scalar Rcomp = 0.0;
                    if (i == 0 && j == 0) Rcomp = Rcell.xx();
                    else if (i == 0 && j == 1) Rcomp = Rcell.xy();
                    else if (i == 0 && j == 2) Rcomp = Rcell.xz();
                    else if (i == 1 && j == 0) Rcomp = Rcell.xy();  // yx = xy
                    else if (i == 1 && j == 1) Rcomp = Rcell.yy();
                    else if (i == 1 && j == 2) Rcomp = Rcell.yz();
                    else if (i == 2 && j == 0) Rcomp = Rcell.xz();  // zx = xz
                    else if (i == 2 && j == 1) Rcomp = Rcell.yz();  // zy = yz
                    else if (i == 2 && j == 2) Rcomp = Rcell.zz();
                    
                    b[row] = Rcomp;
                    row++;
                }
            }
        }
        
        try
        {
            // Perform SVD decomposition using OpenFOAM's SVD class
            SVD svd(A, svdTolerance_);
            
            // Get singular values
            const scalarDiagonalMatrix& S = svd.S();
            
            // Count non-zero singular values for statistics
            label nNonZero = 0;
            for (label i = 0; i < S.size(); ++i)
            {
                if (S[i] > svdTolerance_)
                {
                    nNonZero++;
                }
            }
            
            // Track rank-deficient cases
            if (nNonZero < min(nRows, nCols))
            {
                nRankDeficient++;
            }
            
            // Solve using the SVD decomposition
            const RectangularMatrix<scalar>& U = svd.U();
            const RectangularMatrix<scalar>& V = svd.V();
            
            // Compute solution using SVD: x = V * inv(S) * U.T() * b
            // 1. Compute U.T() * b
            List<scalar> UTb(nCols, 0.0);
            for (label i = 0; i < nCols; ++i)
            {
                for (label j = 0; j < nRows; ++j)
                {
                    UTb[i] += U(j, i) * b[j];
                }
            }
            
            // 2. Compute inv(S) * (U.T() * b)
            List<scalar> SinvUTb(nCols, 0.0);
            for (label i = 0; i < nCols; ++i)
            {
                if (S[i] > svdTolerance_)
                {
                    SinvUTb[i] = UTb[i] / S[i];
                }
                // else: Singular value is too small, contribution is zero
            }
            
            // 3. Compute V * (inv(S) * U.T() * b)
            List<scalar> x(nCols, 0.0);
            for (label i = 0; i < nCols; ++i)
            {
                for (label j = 0; j < nCols; ++j)
                {
                    x[i] += V(i, j) * SinvUTb[j];
                }
            }
            
            // Store coefficients
            cT1_[cellI] = x[0];
            cT2_[cellI] = x[1];
            cT3_[cellI] = x[2];
            cT4_[cellI] = x[3];
            cT5_[cellI] = x[4];
            cT6_[cellI] = x[5];
            cT7_[cellI] = x[6];
            cT8_[cellI] = x[7];
            cT9_[cellI] = x[8];
            cT10_[cellI] = x[9];
            
            // Calculate reconstruction error for statistics
            tensor R_reconstructed = 
                cT1_[cellI] * T1_[cellI] +
                cT2_[cellI] * T2_[cellI] +
                cT3_[cellI] * T3_[cellI] +
                cT4_[cellI] * T4_[cellI] +
                cT5_[cellI] * T5_[cellI] +
                cT6_[cellI] * T6_[cellI] +
                cT7_[cellI] * T7_[cellI] +
                cT8_[cellI] * T8_[cellI] +
                cT9_[cellI] * T9_[cellI] +
                cT10_[cellI] * T10_[cellI];
            
            // Convert reconstructed tensor to symmTensor for comparison
            symmTensor R_recon_symm(
                R_reconstructed.xx(), R_reconstructed.xy(), R_reconstructed.xz(),
                R_reconstructed.yy(), R_reconstructed.yz(),
                R_reconstructed.zz()
            );
            
            scalar error = mag(R_recon_symm - Rcell);
            totalError += error;
            
        }
        catch (const Foam::error& err)
        {
            nRankDeficient++;
            cT1_[cellI] = 0.0;
            cT2_[cellI] = 0.0;
            cT3_[cellI] = 0.0;
            cT4_[cellI] = 0.0;
            cT5_[cellI] = 0.0;
            cT6_[cellI] = 0.0;
            cT7_[cellI] = 0.0;
            cT8_[cellI] = 0.0;
            cT9_[cellI] = 0.0;
            cT10_[cellI] = 0.0;
        }
    }
    
    // Write coefficient statistics
    if (log)
    {
        scalar avgError = 0.0;
        if (nCells > 0)
        {
            avgError = totalError / nCells;
        }
        
        Info<< type() << " " << name() << ": "
            << "Pope Tensor SVD Statistics:" << nl
            << "  Total cells: " << nCells << nl
            << "  Rank-deficient cases: " << nRankDeficient 
            << " (" << (100.0 * nRankDeficient / max(nCells, 1)) << "%)" << nl
            << "  Average reconstruction error: " << avgError << nl
            << "  Coefficient statistics:" << nl
            << "    cT1:  min = " << min(cT1_).value() 
            << ", max = " << max(cT1_).value()
            << ", avg = " << average(cT1_).value() << nl
            << "    cT2:  min = " << min(cT2_).value() 
            << ", max = " << max(cT2_).value()
            << ", avg = " << average(cT2_).value() << nl
            << "    cT3:  min = " << min(cT3_).value() 
            << ", max = " << max(cT3_).value()
            << ", avg = " << average(cT3_).value() << nl
            << "    cT4:  min = " << min(cT4_).value() 
            << ", max = " << max(cT4_).value()
            << ", avg = " << average(cT4_).value() << nl
            << "    cT5:  min = " << min(cT5_).value() 
            << ", max = " << max(cT5_).value()
            << ", avg = " << average(cT5_).value() << nl
            << "    cT6:  min = " << min(cT6_).value() 
            << ", max = " << max(cT6_).value()
            << ", avg = " << average(cT6_).value() << nl
            << "    cT7:  min = " << min(cT7_).value() 
            << ", max = " << max(cT7_).value()
            << ", avg = " << average(cT7_).value() << nl
            << "    cT8:  min = " << min(cT8_).value() 
            << ", max = " << max(cT8_).value()
            << ", avg = " << average(cT8_).value() << nl
            << "    cT9:  min = " << min(cT9_).value() 
            << ", max = " << max(cT9_).value()
            << ", avg = " << average(cT9_).value() << nl
            << "    cT10: min = " << min(cT10_).value() 
            << ", max = " << max(cT10_).value()
            << ", avg = " << average(cT10_).value() << endl;
    }
}


void Foam::functionObjects::reynoldsStressPope::reconstructRPope()
{
    if (!Rptr_)
    {
        FatalErrorInFunction
            << "Reynolds stress tensor field " << RName_ << " not found in mesh"
            << exit(FatalError);
    }
    
    const volSymmTensorField& R = *Rptr_;
    
    if (writeTensors_ && log)
    {
        Info<< type() << " " << name() << ": "
            << "Reconstructing R_Pope from Pope coefficients and tensors" << endl;
    }
    
    // Reconstruct R_Pope = sum_{i=1}^{10} cTi * Ti
    forAll(mesh_.C(), cellI)
    {
        tensor R_pope_tensor = 
            cT1_[cellI] * T1_[cellI] +
            cT2_[cellI] * T2_[cellI] +
            cT3_[cellI] * T3_[cellI] +
            cT4_[cellI] * T4_[cellI] +
            cT5_[cellI] * T5_[cellI] +
            cT6_[cellI] * T6_[cellI] +
            cT7_[cellI] * T7_[cellI] +
            cT8_[cellI] * T8_[cellI] +
            cT9_[cellI] * T9_[cellI] +
            cT10_[cellI] * T10_[cellI];
        
        // Convert to symmTensor (taking only symmetric part)
        R_Pope_[cellI] = symm(R_pope_tensor);
    }
    
    // Calculate and display reconstruction error
    if (log)
    {
        scalarField errorField(mag(R_Pope_ - R));
        scalar maxError = gMax(errorField);
        scalar avgError = gAverage(errorField);
        
        Info<< type() << " " << name() << ": "
            << "Pope reconstruction error: max = " << maxError
            << ", avg = " << avgError << endl;
        
        // Print sample comparison for first cell
        if (mesh_.nCells() > 0)
        {
            Info<< "Sample comparison (cell 0):" << nl
                << "  Original R: " << R[0] << nl
                << "  Reconstructed R_Pope: " << R_Pope_[0] << nl
                << "  Difference: " << R_Pope_[0] - R[0] << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reynoldsStressPope::reynoldsStressPope
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    UName_("U"),
    RName_("R"),
    writeTensors_(true),
    useSymmetricEqns_(false),
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
    r_
    (
        IOobject
        (
            "r",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimVelocity/dimLength, Zero)
    ),
    T1_
    (
        IOobject
        (
            "T1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimVelocity/dimLength, Zero)
    ),
    T2_
    (
        IOobject
        (
            "T2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(sqr(dimVelocity/dimLength), Zero)
    ),
    T3_
    (
        IOobject
        (
            "T3",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(sqr(dimVelocity/dimLength), Zero)
    ),
    T4_
    (
        IOobject
        (
            "T4",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(sqr(dimVelocity/dimLength), Zero)
    ),
    T5_
    (
        IOobject
        (
            "T5",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(pow(dimVelocity/dimLength, 3), Zero)
    ),
    T6_
    (
        IOobject
        (
            "T6",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(pow(dimVelocity/dimLength, 3), Zero)
    ),
    T7_
    (
        IOobject
        (
            "T7",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(pow(dimVelocity/dimLength, 4), Zero)
    ),
    T8_
    (
        IOobject
        (
            "T8",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(pow(dimVelocity/dimLength, 4), Zero)
    ),
    T9_
    (
        IOobject
        (
            "T9",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(pow(dimVelocity/dimLength, 4), Zero)
    ),
    T10_
    (
        IOobject
        (
            "T10",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(pow(dimVelocity/dimLength, 5), Zero)
    ),
    Rptr_(nullptr),
    cT1_
    (
        IOobject
        (
            "cT1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength/dimTime, Zero)
    ),
    cT2_
    (
        IOobject
        (
            "cT2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength, Zero)
    ),
    cT3_
    (
        IOobject
        (
            "cT3",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength, Zero)
    ),
    cT4_
    (
        IOobject
        (
            "cT4",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength, Zero)
    ),
    cT5_
    (
        IOobject
        (
            "cT5",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength*dimTime, Zero)
    ),
    cT6_
    (
        IOobject
        (
            "cT6",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength*dimTime, Zero)
    ),
    cT7_
    (
        IOobject
        (
            "cT7",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength*dimTime*dimTime, Zero)
    ),
    cT8_
    (
        IOobject
        (
            "cT8",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength*dimTime*dimTime, Zero)
    ),
    cT9_
    (
        IOobject
        (
            "cT9",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength*dimTime*dimTime, Zero)
    ),
    cT10_
    (
        IOobject
        (
            "cT10",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength*dimLength*dimTime*dimTime*dimTime, Zero)
    ),
    R_Pope_
    (
        IOobject
        (
            "R_Pope",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor(sqr(dimVelocity), Zero)
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::reynoldsStressPope::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    
    // Read field names
    dict.readIfPresent("U", UName_);
    dict.readIfPresent("R", RName_);
    dict.readIfPresent("writeTensors", writeTensors_);
    dict.readIfPresent("useSymmetricEqns", useSymmetricEqns_);
    dict.readIfPresent("svdTolerance", svdTolerance_);
    
    // Check if fields exist
    if (!mesh_.foundObject<volVectorField>(UName_))
    {
        FatalErrorInFunction
            << "Velocity field " << UName_ << " not found in mesh"
            << exit(FatalError);
    }
    
    // Find Reynolds stress tensor field (it should exist in the simulation)
    Rptr_ = mesh_.findObject<volSymmTensorField>(RName_);
    if (!Rptr_)
    {
        FatalErrorInFunction
            << "Reynolds stress tensor field " << RName_ << " not found in mesh. "
            << "Make sure you are running a simulation that computes R "
            << "(e.g., LES, RSM, or post-processing from velocity fluctuations)."
            << exit(FatalError);
    }
    
    if (log)
    {
        Info<< type() << " " << name() << ":" << nl
            << "    U: " << UName_ << nl
            << "    R: " << RName_ << nl
            << "    writeTensors: " << writeTensors_ << nl
            << "    useSymmetricEqns: " << useSymmetricEqns_ << nl
            << "    SVD tolerance: " << svdTolerance_ << endl;
    }
    
    return true;
}


bool Foam::functionObjects::reynoldsStressPope::execute()
{
    if (log)
    {
        Info<< type() << " " << name() << " executing..." << endl;
    }
    
    // Calculate strain rate and rotation rate tensors
    calcStrainAndRotation();
    
    // Calculate Pope's tensors T1-T10
    calcPopeTensors();
    
    // Perform SVD decomposition to find Pope coefficients
    performSVD();
    
    // Reconstruct R_Pope from coefficients
    reconstructRPope();
    
    return true;
}


bool Foam::functionObjects::reynoldsStressPope::write()
{
    if (log)
    {
        Info<< type() << " " << name() << " writing fields:" << endl;
    }
    
    // Always write coefficients and reconstructed field
    cT1_.write(); cT2_.write(); cT3_.write(); cT4_.write(); cT5_.write();
    cT6_.write(); cT7_.write(); cT8_.write(); cT9_.write(); cT10_.write();
    R_Pope_.write();
    
    // Write intermediate tensors if requested
    if (writeTensors_)
    {
        s_.write();
        r_.write();
        T1_.write(); T2_.write(); T3_.write(); T4_.write(); T5_.write();
        T6_.write(); T7_.write(); T8_.write(); T9_.write(); T10_.write();
    }
    
    return true;
}


// ************************************************************************* //
