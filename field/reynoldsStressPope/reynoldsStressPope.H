/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Function object to calculate Pope's 10 tensors T1-T10 from velocity field,
    read Reynolds stress tensor R (as symmTensorField) from simulation, 
    perform SVD decomposition to find coefficients cT1-cT10 in the formula:
    R = cT1*T1 + cT2*T2 + ... + cT10*T10
    
    Pope's tensors are defined as:
    T1 = s
    T2 = s & r - r & s
    T3 = s & s - I * tr(s & s) / 3
    T4 = r & r - I * tr(r & r) / 3
    T5 = r & s & s - s & s & r
    T6 = r & r & s + s & r & r - I * 2 * tr(s & r & r) / 3
    T7 = r & s & r & r - r & r & s & r
    T8 = s & r & s & s - s & s & r & s
    T9 = r & r & s & s + s & s & r & r - I * 2 * tr(s & s & r & r) / 3
    T10 = r & s & s & r & r - r & r & s & s & r
    
    where:
    s = 0.5*(grad(U) + grad(U).T())  (strain rate tensor)
    r = 0.5*(grad(U) - grad(U).T())  (rotation rate tensor)

Usage
    Example of function object specification:
    \verbatim
    reynoldsStressPope
    {
        type            reynoldsStressPope;
        libs            ("libfieldFunctionObjects");
        enabled         true;
        writeControl    timeStep;
        writeInterval   1;
        
        // Name of Reynolds stress tensor field (must exist in simulation)
        R               R; // Default is "R"
        
        // Name of velocity field
        U               U; // Default is "U"
        
        // Optional: whether to write intermediate tensors
        writeTensors    true;
        
        // Optional: SVD tolerance
        svdTolerance    1e-6; // Default is 1e-6
        
        // Optional: use only 6 independent equations for symmetry
        useSymmetricEqns false; // Default is false (use 9 equations)
    }
    \endverbatim

SourceFiles
    reynoldsStressPope.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_reynoldsStressPope_H
#define functionObjects_reynoldsStressPope_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "OFstream.H"
#include "SVD.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class reynoldsStressPope Declaration
\*---------------------------------------------------------------------------*/

class reynoldsStressPope
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of velocity field
        word UName_;
        
        //- Name of Reynolds stress tensor field
        word RName_;
        
        //- Write intermediate tensor fields
        bool writeTensors_;
        
        //- Use only 6 independent equations (for symmetric tensor)
        bool useSymmetricEqns_;
        
        //- SVD tolerance
        scalar svdTolerance_;
        
        //- Strain rate tensor field
        volTensorField s_;
        
        //- Rotation rate tensor field
        volTensorField r_;
        
        //- Pope's tensors T1-T10
        volTensorField T1_, T2_, T3_, T4_, T5_, T6_, T7_, T8_, T9_, T10_;
        
        //- Reynolds stress tensor field (read from simulation)
        const volSymmTensorField* Rptr_;
        
        //- Pope coefficients cT1-cT10
        volScalarField cT1_, cT2_, cT3_, cT4_, cT5_;
        volScalarField cT6_, cT7_, cT8_, cT9_, cT10_;
        
        //- Reconstructed Reynolds stress tensor from Pope's tensors
        volSymmTensorField R_Pope_;


    // Private Member Functions

        //- Calculate strain rate and rotation rate tensors
        void calcStrainAndRotation();
        
        //- Calculate Pope's tensors T1-T10
        void calcPopeTensors();
        
        //- Perform SVD decomposition in each cell to find Pope coefficients
        void performSVD();
        
        //- Reconstruct R_Pope from Pope coefficients and tensors
        void reconstructRPope();


public:

    //- Runtime type information
    TypeName("reynoldsStressPope");


    // Constructors

        //- Construct from Time and dictionary
        reynoldsStressPope
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        reynoldsStressPope(const reynoldsStressPope&) = delete;

        //- No copy assignment
        void operator=(const reynoldsStressPope&) = delete;


    // Destructor
    virtual ~reynoldsStressPope() = default;


    // Member Functions

        //- Read the reynoldsStressPope data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the reynoldsStressPope
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
