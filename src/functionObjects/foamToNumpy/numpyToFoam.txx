/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tanuj Ravi
    Copyright (C) 2026 Keysight Technologies
    Copyright (C) 2026 Andre Weiner
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "Pstream.H"
#include "fvMesh.H"
#include "numpyFileReader.H"

template<class GeoField>
void Foam::functionObjects::numpyToFoam::correctFieldBoundaryConditions
(
    GeoField& field
) const
{
    bool hasEquationCoupledPatch = false;

    for (const auto& patchField : field.boundaryField())
    {
        const word& patchType = patchField.type();

        if
        (
            patchType.find("epsilonWallFunction") != string::npos
         || patchType.find("omegaWallFunction") != string::npos
         || patchType == "fixedFluxPressure"
        )
        {
            hasEquationCoupledPatch = true;
            break;
        }
    }

    if (!hasEquationCoupledPatch)
    {
        field.correctBoundaryConditions();
        return;
    }

    Log << type() << ' ' << name()
        << ": retaining template values for equation-coupled patches of "
        << field.name() << endl;

    for (auto& patchField : field.boundaryFieldRef())
    {
        const word& patchType = patchField.type();

        if
        (
            patchType.find("epsilonWallFunction") != string::npos
         || patchType.find("omegaWallFunction") != string::npos
         || patchType == "fixedFluxPressure"
        )
        {
            continue;
        }

        patchField.setUpdated(false);
        patchField.initEvaluate(Pstream::defaultCommsType);
    }

    for (auto& patchField : field.boundaryFieldRef())
    {
        const word& patchType = patchField.type();

        if
        (
            patchType.find("epsilonWallFunction") == string::npos
         && patchType.find("omegaWallFunction") == string::npos
         && patchType != "fixedFluxPressure"
        )
        {
            patchField.evaluate(Pstream::defaultCommsType);
        }
    }
}


template<class GeoField>
bool Foam::functionObjects::numpyToFoam::writeField(const word& fieldName)
{
    GeoField* fieldPtr = mesh_.getObjectPtr<GeoField>(fieldName);

    if (!fieldPtr)
    {
        return false;
    }

    const word oldInstance(fieldPtr->instance());
    fieldPtr->instance() = time_.timeName();
    const bool written = fieldPtr->write();
    fieldPtr->instance() = oldInstance;
    return written;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

template<class GeoField>
GeoField& Foam::functionObjects::numpyToFoam::field
(
    const word& fieldName,
    HashPtrTable<GeoField>& ownedFields
)
{
    if (GeoField* fieldPtr = mesh_.getObjectPtr<GeoField>(fieldName))
    {
        return *fieldPtr;
    }

    const IOobject templateIO
    (
        fieldName,
        templateInstance_,
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    const GeoField templateField(templateIO, mesh_);

    GeoField* fieldPtr = new GeoField
    (
        IOobject
        (
            fieldName,
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        templateField
    );

    ownedFields.set(fieldName, fieldPtr);
    return *fieldPtr;
}


template<class GeoField>
void Foam::functionObjects::numpyToFoam::loadField
(
    const numpyDetail::numpyInputCatalog::snapshot& sample,
    const word& fieldName,
    HashPtrTable<GeoField>& ownedFields,
    const bool validateOnly
)
{
    if (zoneMode_)
    {
        for (const auto& zone : zones_)
        {
            numpyDetail::numpyFileReader reader
                (catalog_->fieldPath(sample, fieldName, Pstream::myProcNo(), zone.name));
            Field<typename GeoField::value_type> values;
            reader.readField(sample.index, zone.cells.size(), values);
            if (!validateOnly)
            {
                GeoField& target = field<GeoField>(fieldName, ownedFields);
                auto& internal = target.primitiveFieldRef();
                forAll(values, i) internal[zone.cells[i]] = values[i];
            }
        }
        if (!validateOnly)
        {
            GeoField& target = field<GeoField>(fieldName, ownedFields);
            if (correctBoundaryConditions_) correctFieldBoundaryConditions(target);
            if (!target.nOldTimes()) (void)target.oldTime();
        }
        return;
    }
    GeoField& target = field<GeoField>(fieldName, ownedFields);
    const fileName inputFile
    (
        catalog_->fieldPath(sample, fieldName, Pstream::myProcNo())
    );

    if (!isFile(inputFile))
    {
        FatalErrorInFunction
            << "Cannot find NumPy field file " << inputFile
            << exit(FatalError);
    }

    numpyDetail::numpyFileReader reader(inputFile);
    Field<typename GeoField::value_type> values;
    reader.readField(sample.index, target.size(), values);

    // Access-time tracking stores the preceding imported snapshot in
    // oldTime() after it has been allocated below.
    target.primitiveFieldRef() = values;

    if (correctBoundaryConditions_)
    {
        correctFieldBoundaryConditions(target);
    }

    if (!target.nOldTimes())
    {
        (void)target.oldTime();
    }
}


// ************************************************************************* //
