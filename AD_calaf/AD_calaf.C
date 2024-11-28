/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2018-2020 OpenCFD Ltd
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

#include "AD_calaf.H"
#include "geometricOneField.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(AD_calaf, 0);
    addToRunTimeSelectionTable(option, AD_calaf, dictionary);
}
}

// Define AD variant
const Foam::Enum
<
    Foam::fv::AD_calaf::forceMethodType
>
Foam::fv::AD_calaf::forceMethodTypeNames
({
//    { forceMethodType::UNIFORM_THRUST_AND_TANGENTIAL, "uniformThrustTangential" },
    { forceMethodType::UNIFORM_THRUST, "uniformThrust" },
});



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// Write header for the AD output file.
void Foam::fv::AD_calaf::writeFileHeader(Ostream& os)
{
    // Default output.
    writeFile::writeHeader(os, "AD based on disk-averaged velocities");
    writeFile::writeCommented(os, "Time");


    if (forceMethod_ == forceMethodType::UNIFORM_THRUST)
    {
        // Additional output for the uniform thrust AD.
        writeFile::writeCommented(os, "Uref");
        writeFile::writeCommented(os, "Cp");
        writeFile::writeCommented(os, "Ct");
        writeFile::writeCommented(os, "Udisk");
        writeFile::writeCommented(os, "CpStar");
        writeFile::writeCommented(os, "CtStar");
        writeFile::writeCommented(os, "T");
        writeFile::writeCommented(os, "P");
    }
    //    else if (forceMethod_ == forceMethodType::UNIFORM_THRUST_AND_TANGENTIAL)
    //    {
    //        // To be made: Could add torque and TSR.
    //        writeFile::writeCommented(os, "M");
    //        writeFile::writeCommented(os, "TSR");
    //    }

    os  << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::AD_calaf::AD_calaf
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    writeFile(mesh, name, modelType, coeffs_),
    forceMethod_
    (
        forceMethodTypeNames.getOrDefault
        (
            // Retrieve "variant" from the coeffs dictionary and if not found defaults to uniform thrust.
            "variant",
            coeffs_,
            forceMethodType::UNIFORM_THRUST
        )
    ),
    writeFileStart_(coeffs_.getOrDefault<scalar>("writeFileStart", 0)),
    writeFileEnd_(coeffs_.getOrDefault<scalar>("writeFileEnd", VGREAT)),
    diskArea_
    (
        coeffs_.getCheck<scalar>
        (
            "diskArea",
            scalarMinMax::ge(VSMALL)
        )
    ),
    diskDir_
    (
        coeffs_.getCheck<vector>
        (
            "diskDir",
            [&](const vector& vec){ return mag(vec) > VSMALL; }
        ).normalise()
    ),
    Ct_(coeffs_.get<scalar>("Ct"))

// Body of constructor (where we run some functions).
{
    fieldNames_.setSize(1, "U");

    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating actuation disk zone: " << this->name() << endl;

    Info<< "    - force computation method: "
        << forceMethodTypeNames[forceMethod_] << endl;

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Incompressible (add momentum source)
void Foam::fv::AD_calaf::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (V() > VSMALL)
    {
        calc(eqn);
    }
}

// Function that reads stuff from fvOptions file and checks that they are valid.
bool Foam::fv::AD_calaf::read(const dictionary& dict)
{
    if (cellSetOption::read(dict) && writeFile::read(dict))
    {
        dict.readIfPresent("writeFileStart", writeFileStart_);
        dict.readIfPresent("writeFileEnd", writeFileEnd_);
        dict.readIfPresent("diskArea", diskArea_);
        if (diskArea_ < VSMALL)
        {
            FatalErrorInFunction
                << "Actuator disk has zero area: "
                << "diskArea = " << diskArea_
                << exit(FatalIOError);
        }

        dict.readIfPresent("diskDir", diskDir_);
        diskDir_.normalise();
        if (mag(diskDir_) < VSMALL)
        {
            FatalErrorInFunction
                << "Actuator disk surface-normal vector is zero: "
                << "diskDir = " << diskDir_
                << exit(FatalIOError);
        }

        return true;
    }

    return false;
}




// Dispatcher function, which choose the AD variant.
void Foam::fv::AD_calaf::calc(
    fvMatrix<vector> &eqn)
{
    switch (forceMethod_)
    {

    case forceMethodType::UNIFORM_THRUST:
    {
        calcUniformThrustMethod(eqn);
        break;
    }

    default:
        break;
    }
}

// The uniform thrust AD.
void Foam::fv::AD_calaf::calcUniformThrustMethod(
    fvMatrix<vector> &eqn)
{
    // Velocity field.
    const vectorField &U = eqn.psi();
    // Source term field (to be calculated).
    vectorField &Usource = eqn.source();
    // Mesh.
    const scalarField &cellsV = mesh_.V();

    // Calculate disk-averaged quantities.
    vector Udisk(Zero);
    scalar totalV = 0.0;
    for (const auto &celli : cells_)
    {
        Udisk += U[celli] * cellsV[celli];
        totalV += cellsV[celli];
    }
    reduce(Udisk, sumOp<vector>());
    reduce(totalV, sumOp<scalar>());
    if (totalV < SMALL)
    {
        FatalErrorInFunction
            << "No cell in the actuator disk."
            << exit(FatalError);
    }
    Udisk /= totalV;
    const scalar magUdisk = mag(Udisk);
    if (mag(Udisk) < SMALL)
    {
        FatalErrorInFunction
            << "Velocity spatial-averaged on actuator disk is zero." << nl
            << "Please check if the initial U field is zero."
            << exit(FatalError);
    }

    // 1D mom theory (only depends on CT; CP is automatically calculated!)
    const scalar Ct = Ct_;
    if (Ct <= VSMALL)
    {
        FatalErrorInFunction
            << "Ct must be greater than zero." << nl
            << ", Ct = " << Ct
            << exit(FatalIOError);
    }
    const scalar a = 0.5 * (1.0 - sqrt(1 - Ct));
    const scalar CtStar = Ct * 1.0 / sqr(1.0 - a);
    const scalar Cp = 4 * a * sqr(1.0 - a);
    const scalar CpStar = Cp * 1.0 / pow3(1.0 - a);
    const scalar magUref = magUdisk / (1.0 - a);

    // Compute thrust and power
    const scalar T = 0.5 * diskArea_ * magSqr(Udisk & diskDir_) * CtStar;
    const scalar P = 0.5 * diskArea_ * pow3(mag(Udisk & diskDir_)) * CpStar;

    // Calculate momentum source term.
    for (const label celli : cells_)
    {
        Usource[celli] += (cellsV[celli] / totalV * T) * diskDir_;
    }

    // Write disk quantities to file.
    if (
        mesh_.time().timeOutputValue() >= writeFileStart_ && mesh_.time().timeOutputValue() <= writeFileEnd_)
    {
        Ostream &os = file();
        writeCurrentTime(os);

        // Output values of the current timestep.
        os << magUref << tab << Cp << tab << Ct << tab
           << magUdisk << tab << CpStar << tab << CtStar << tab << T << tab << P
           << endl;
    }
}

// ************************************************************************* //
