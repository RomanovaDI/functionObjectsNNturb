/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Function object to calculate strain rate tensor, its powers, Reynolds stress
    tensor, perform SVD decomposition to find coefficients c_s, c_ss, c_sss,
    and reconstruct R_SVD = c_s*s + c_ss*ss + c_sss*sss.

    Solves: R = c_s * s + c_ss * ss + c_sss * sss using SVD in each cell.

Usage
    Example of function object specification:
    \verbatim
    reynoldsStressSVD
    {
        type            reynoldsStressSVD;
        libs            ("libfieldFunctionObjects");
        enabled         true;
        writeControl    timeStep;
        writeInterval   1;
        
        // Optional: turbulence kinetic energy field name
        k               k; // Default is "k"
        
        // Optional: whether to write intermediate tensors
        writeTensors    true;
        
        // Optional: SVD tolerance
        svdTolerance    1e-6; // Default is 1e-6
    }
    \endverbatim

SourceFiles
    reynoldsStressSVD.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_reynoldsStressSVD_H
#define functionObjects_reynoldsStressSVD_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "OFstream.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class reynoldsStressSVD Declaration
\*---------------------------------------------------------------------------*/

class reynoldsStressSVD
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of velocity field
        word UName_;
        
        //- Name of turbulent viscosity field
        word nutName_;
        
        //- Name of turbulence kinetic energy field
        word kName_;
        
        //- Write intermediate tensor fields
        bool writeTensors_;
        
        //- SVD tolerance
        scalar svdTolerance_;
        
        //- Strain rate tensor field
        volTensorField s_;
        
        //- Square of strain rate tensor
        volTensorField ss_;
        
        //- Cube of strain rate tensor
        volTensorField sss_;
        
        //- Reynolds stress tensor field
        volTensorField R_;
        
        //- Coefficient c_s field
        volScalarField c_s_;
        
        //- Coefficient c_ss field
        volScalarField c_ss_;
        
        //- Coefficient c_sss field
        volScalarField c_sss_;
        
        //- Reconstructed Reynolds stress tensor from SVD
        volTensorField R_SVD_;


    // Private Member Functions

        //- Calculate strain rate tensor
        void calcStrainRate();
        
        //- Calculate Reynolds stress tensor using Boussinesq hypothesis
        void calcReynoldsStress();
        
        //- Perform SVD decomposition in each cell
        void performSVD();
        
        //- Reconstruct R_SVD from coefficients
        void reconstructRSVD();


public:

    //- Runtime type information
    TypeName("reynoldsStressSVD");


    // Constructors

        //- Construct from Time and dictionary
        reynoldsStressSVD
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        reynoldsStressSVD(const reynoldsStressSVD&) = delete;

        //- No copy assignment
        void operator=(const reynoldsStressSVD&) = delete;


    // Destructor
    virtual ~reynoldsStressSVD() = default;


    // Member Functions

        //- Read the reynoldsStressSVD data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the reynoldsStressSVD
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
