/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
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

Class
    Foam::fv::AD_calaf

Group
    grpFvOptionsSources

Description
    To be done.

\*---------------------------------------------------------------------------*/

#ifndef AD_calaf_H
#define AD_calaf_H

#include "cellSetOption.H"
#include "writeFile.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class AD_calaf Declaration
\*---------------------------------------------------------------------------*/

// Define AD_calaf class, which inherits functionality from cellSetOption and writeFile.
class AD_calaf
:
    public cellSetOption,
    public functionObjects::writeFile
{
protected:

    // Variants of AD models
    enum forceMethodType
    {
        UNIFORM_THRUST, //!< "Uniform thrust distribution on the AD"
//        UNIFORM_THRUST_TANGENTIAL,           //!< "Uniform distribution of thrust and tangential force"
    };
    static const Enum<forceMethodType> forceMethodTypeNames;

    //- The type of the force computation method
    const enum forceMethodType forceMethod_;

    //- Start time for file output
    scalar writeFileStart_;

    //- End time for file output
    scalar writeFileEnd_;

    //- Actuator disk planar surface area [m2]
    scalar diskArea_;

    //- Surface-normal vector of the actuator disk pointing upstream
    vector diskDir_;

    //- Thrust coefficient
    scalar Ct_;

    //- Function to output file header information
    virtual void writeFileHeader(Ostream& os);


private:
    // Private Member Functions

    //- Select force computation method
    void calc
    (
        fvMatrix<vector>& eqn
    );

    //- Compute uniform thurst AD method
    void calcUniformThrustMethod
    (
        fvMatrix<vector>& eqn
    );


public:

    //- Runtime type information
    TypeName("AD_calaf");

    // Constructor: When initilizing a AD_calaf class, you need to specify these!
    AD_calaf
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- No copy construct
    AD_calaf(const AD_calaf&) = delete;

    //- No copy assignment
    void operator=(const AD_calaf&) = delete;

    //- Destructor
    virtual ~AD_calaf() = default;


    // ----------  Member Functions ----------------

    //- Normal disk direction
    const vector& diskDir() const
    {
        return diskDir_;
    }

    //- Disk area
    scalar diskArea() const
    {
        return diskArea_;
    }

    //- Source term to momentum equation
    virtual void addSup
    (
        fvMatrix<vector>& eqn,
        const label fieldi
    );

    //- IO. Read dictionary
    virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
