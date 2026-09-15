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

#include "numpyToFoam.H"
#include "IFstream.H"
#include "numpyFileReader.H"
#include "PstreamReduceOps.H"

#include "IOobject.H"
#include "ListOps.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(numpyToFoam, 0);
    addToRunTimeSelectionTable(functionObject, numpyToFoam, dictionary);
}
}


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::functionObjects::numpyToFoam::numpyToFoam
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldNames_(),
    templateInstance_("0"),
    writeFields_(false),
    correctBoundaryConditions_(false),
    catalog_(),
    lastTime_(-VGREAT),
    scalarFields_(),
    vectorFields_(),
    sphericalTensorFields_(),
    symmTensorFields_(),
    tensorFields_()
{
    read(dict);
}


Foam::functionObjects::numpyToFoam::~numpyToFoam()
{}


Foam::functionObjects::numpyToFoam::numpyToFoam
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, obr, dict),
    fieldNames_(),
    templateInstance_("0"),
    writeFields_(false),
    correctBoundaryConditions_(false),
    catalog_(),
    lastTime_(-VGREAT),
    scalarFields_(),
    vectorFields_(),
    sphericalTensorFields_(),
    symmTensorFields_(),
    tensorFields_()
{
    read(dict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::word Foam::functionObjects::numpyToFoam::fieldClass
(
    const word& fieldName
) const
{
    if (const regIOobject* objectPtr = mesh_.cfindObject<regIOobject>(fieldName))
    {
        return objectPtr->type();
    }

    IOobject templateIO
    (
        fieldName,
        templateInstance_,
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!templateIO.typeHeaderOk<regIOobject>(false))
    {
        FatalErrorInFunction
            << "Cannot read a field template for " << fieldName << " from "
            << templateInstance_ << nl
            << "Imported fields require a template for dimensions and "
            << "boundary conditions unless a compatible field is already "
            << "registered."
            << exit(FatalError);
    }

    return templateIO.headerClassName();
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::numpyToFoam::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    const bool newZoneMode = dict.found("cellZones");
    wordRes newZones;
    if (newZoneMode)
    {
        dict.readEntry("cellZones", newZones);
        if (newZones.empty())
            FatalIOErrorInFunction(dict) << "cellZones must not be empty" << exit(FatalIOError);
    }
    if (catalog_ && (newZoneMode != zoneMode_ || newZones != zoneSelection_))
        FatalIOErrorInFunction(dict) << "Input cellZones cannot change after construction"
            << exit(FatalIOError);
    zoneMode_ = newZoneMode;
    zoneSelection_ = newZones;
    wordList newFieldNames(dict.get<wordList>("fields"));
    inplaceUniqueSort(newFieldNames);

    const word newTemplateInstance
    (
        dict.getOrDefault<word>("templateInstance", "0")
    );
    const bool newWriteFields = dict.getOrDefault("writeFields", false);
    const bool newCorrectBoundaryConditions =
        dict.getOrDefault("correctBoundaryConditions", false);

    if
    (
        catalog_
     &&
        (
            newFieldNames != fieldNames_
         || newTemplateInstance != templateInstance_
        )
    )
    {
        FatalIOErrorInFunction(dict)
            << "Input fields and templates cannot change after construction"
            << exit(FatalIOError);
    }

    fieldNames_ = std::move(newFieldNames);
    templateInstance_ = newTemplateInstance;
    writeFields_ = newWriteFields;
    correctBoundaryConditions_ = newCorrectBoundaryConditions;

    if (!catalog_)
    {
        catalog_.reset
        (
            new numpyDetail::numpyInputCatalog(time_.globalPath(), dict)
        );
    }

    return true;
}


bool Foam::functionObjects::numpyToFoam::load(const scalar timeValue)
{
    if (mag(timeValue - lastTime_) <= SMALL)
    {
        return true;
    }

    if (timeValue < lastTime_)
    {
        FatalErrorInFunction
            << "Import time moved backwards from " << lastTime_
            << " to " << timeValue
            << exit(FatalError);
    }

    const numpyDetail::numpyInputCatalog::snapshot& sample =
        catalog_->find(timeValue);

    validateZones(sample);
    // Validate all selected payloads before modifying any registered field.
    for (int pass = zoneMode_ ? 0 : 1; pass < 2; ++pass)
    {
        const bool validateOnly = pass == 0;
        for (const word& fieldName : fieldNames_)
        {
            const word className(fieldClass(fieldName));

            if (className == volScalarField::typeName)
            {
                loadField(sample, fieldName, scalarFields_, validateOnly);
            }
            else if (className == volVectorField::typeName)
            {
                loadField(sample, fieldName, vectorFields_, validateOnly);
            }
            else if (className == volSphericalTensorField::typeName)
            {
                loadField(sample, fieldName, sphericalTensorFields_, validateOnly);
            }
            else if (className == volSymmTensorField::typeName)
            {
                loadField(sample, fieldName, symmTensorFields_, validateOnly);
            }
            else if (className == volTensorField::typeName)
            {
                loadField(sample, fieldName, tensorFields_, validateOnly);
            }
            else
            {
                FatalErrorInFunction
                    << "Unsupported class " << className
                    << " for imported field " << fieldName
                    << exit(FatalError);
            }
        }

        UPstream::barrier(UPstream::worldComm);
    }

    lastTime_ = timeValue;
    Log << type() << ' ' << name() << ": imported time "
        << Time::timeName(timeValue) << endl;
    return true;
}


bool Foam::functionObjects::numpyToFoam::execute()
{
    return load(time_.value());
}


void Foam::functionObjects::numpyToFoam::correctBoundaryConditions()
{
    for (const word& fieldName : fieldNames_)
    {
        const word className(fieldClass(fieldName));

        if (className == volScalarField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volScalarField>(fieldName)
            );
        }
        else if (className == volVectorField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volVectorField>(fieldName)
            );
        }
        else if (className == volSphericalTensorField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volSphericalTensorField>(fieldName)
            );
        }
        else if (className == volSymmTensorField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volSymmTensorField>(fieldName)
            );
        }
        else if (className == volTensorField::typeName)
        {
            correctFieldBoundaryConditions
            (
                mesh_.lookupObjectRef<volTensorField>(fieldName)
            );
        }
    }
}


bool Foam::functionObjects::numpyToFoam::write()
{
    if (!writeFields_)
    {
        return true;
    }

    for (const word& fieldName : fieldNames_)
    {
        const word className(fieldClass(fieldName));
        bool written = false;

        if (className == volScalarField::typeName)
        {
            written = writeField<volScalarField>(fieldName);
        }
        else if (className == volVectorField::typeName)
        {
            written = writeField<volVectorField>(fieldName);
        }
        else if (className == volSphericalTensorField::typeName)
        {
            written = writeField<volSphericalTensorField>(fieldName);
        }
        else if (className == volSymmTensorField::typeName)
        {
            written = writeField<volSymmTensorField>(fieldName);
        }
        else if (className == volTensorField::typeName)
        {
            written = writeField<volTensorField>(fieldName);
        }

        if (!written)
        {
            FatalErrorInFunction
                << "Cannot write imported field " << fieldName
                << exit(FatalError);
        }
    }

    return true;
}


Foam::scalarList Foam::functionObjects::numpyToFoam::times() const
{
    return catalog_->times();
}


void Foam::functionObjects::numpyToFoam::validateZones
(const numpyDetail::numpyInputCatalog::snapshot& sample)
{
    IFstream segmentStream(sample.batchPath.path()/"segmentInfo");
    dictionary segment(segmentStream);
    IFstream stateStream(sample.batchPath/"state");
    dictionary state(stateStream);
    const bool zoneData = segment.found("zoneLayoutVersion");
    if (zoneData != zoneMode_ || state.found("zoneLayoutVersion") != zoneData)
        FatalErrorInFunction << "Zone datasets require explicit cellZones; "
            << "whole-mesh datasets cannot be imported as zones" << exit(FatalError);
    if (!zoneData) return;
    if (segment.get<label>("zoneLayoutVersion") != 1
        || state.get<label>("zoneLayoutVersion") != 1)
        FatalErrorInFunction << "Unsupported zone layout version" << exit(FatalError);
    if (segment.get<word>("region") != mesh_.name()
        || segment.get<label>("nProcs") != Pstream::nProcs())
        FatalErrorInFunction << "Zone mesh region or processor count mismatch"
            << exit(FatalError);
    const wordList exported(state.get<wordList>("cellZones"));
    const wordList selected(numpyDetail::selectCellZoneNames(zoneSelection_, exported));
    zones_ = numpyDetail::cellZoneMappings(mesh_, zoneSelection_);
    wordList targetNames;
    for (const auto& zone : zones_) targetNames.push_back(zone.name);
    if (targetNames != selected)
        FatalErrorInFunction << "Source and target zone selections differ" << exit(FatalError);
    const fileName geometry(catalog_->geometryPath(sample));
    IFstream meshStream(geometry/("mesh_proc_" + Foam::name(Pstream::myProcNo())));
    dictionary metadata(meshStream);
    if (metadata.get<label>("nCells") != mesh_.nCells()
        || metadata.get<label>("meshRevision") != sample.meshRevision)
        FatalErrorInFunction << "Zone mesh size or mapping revision mismatch" << exit(FatalError);
    labelHashSet occupied;
    bool overlap = false;
    for (const auto& zone : zones_)
    {
        numpyDetail::numpyFileReader reader
        (
            geometry/"cellZones"/zone.name
            /("cellIds_proc_" + Foam::name(Pstream::myProcNo()) + ".npy")
        );
        const labelList ids(reader.readCellIds());
        if (ids != zone.cells)
            FatalErrorInFunction << "Cell addressing mismatch for zone " << zone.name
                << "; the original numbering and membership are required" << exit(FatalError);
        for (label id : ids) if (!occupied.insert(id)) overlap = true;
    }
    reduce(overlap, orOp<bool>());
    if (overlap)
        FatalErrorInFunction << "Imported cellZones overlap" << exit(FatalError);
    const dictionary& classes = state.subDict("fieldClasses");
    for (const word& name : fieldNames_)
        if (classes.get<word>(name) != fieldClass(name))
            FatalErrorInFunction << "Field class mismatch for " << name << exit(FatalError);
}

// ************************************************************************* //
