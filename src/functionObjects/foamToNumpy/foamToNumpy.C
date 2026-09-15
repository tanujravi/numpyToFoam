/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tanuj Ravi
    Copyright (C) 2026 Keysight Technologies
    Copyright (C) 2026 Andre Weiner
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

#include "foamToNumpy.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "volFields.H"

#include <algorithm>
#include <iomanip>
#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(foamToNumpy, 0);
    addToRunTimeSelectionTable(functionObject, foamToNumpy, dictionary);
}
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::foamToNumpy::foamToNumpy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSelection_(),
    dtype_(numpyDetail::dataType::FLOAT64),
    batchSize_(100),
    writeTimes_(true),
    writeCellCentres_(false),
    writeCellVolumes_(false),
    outputPath_(),
    segmentPath_(),
    batchPath_(),
    segmentReady_(false),
    batchOpen_(false),
    geometryDirty_(true),
    batchIndex_(0),
    batchCount_(0),
    meshRevision_(0),
    firstBatchTime_(VGREAT),
    lastBatchTime_(-VGREAT),
    lastOutputTime_(-VGREAT),
    batchFields_(),
    fieldWriters_(),
    timesWriter_()
{
    read(dict);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

std::vector<Foam::functionObjects::foamToNumpy::fieldInfo>
Foam::functionObjects::foamToNumpy::selectedFields() const
{
    std::vector<fieldInfo> info;
    HashSet<word> selected;

    appendFieldInfo<volScalarField>(info, selected);
    appendFieldInfo<volVectorField>(info, selected);
    appendFieldInfo<volSphericalTensorField>(info, selected);
    appendFieldInfo<volSymmTensorField>(info, selected);
    appendFieldInfo<volTensorField>(info, selected);

    std::sort
    (
        info.begin(),
        info.end(),
        [](const fieldInfo& a, const fieldInfo& b)
        {
            return a.name < b.name;
        }
    );

    if (info.empty())
    {
        FatalErrorInFunction
            << "No supported volume fields match " << fieldSelection_ << nl
            << "Available objects: "
            << mesh_.objectRegistry::sortedToc()
            << exit(FatalError);
    }

    return info;
}


void Foam::functionObjects::foamToNumpy::ensureSegment()
{
    if (segmentReady_)
    {
        return;
    }

    fileName segmentName
    (
        functionObject::postProcess
      ? time_.timeName()
      : Time::timeName(time_.startTime().value())
    );

    if (Pstream::master())
    {
        const fileName initialName(segmentName);
        label suffix = 0;

        while (isDir(outputPath_/segmentName))
        {
            segmentName = initialName + '_' + Foam::name(++suffix);
        }
    }

    Pstream::broadcast(segmentName);
    segmentPath_ = outputPath_/segmentName;

    if (Pstream::master())
    {
        mkDir(segmentPath_);
    }

    UPstream::barrier(UPstream::worldComm);
    segmentReady_ = true;
    writeSegmentInfo();

    Log << type() << ' ' << name() << ": writing to "
        << segmentPath_ << endl;
}


void Foam::functionObjects::foamToNumpy::beginBatch
(
    const std::vector<fieldInfo>& fields
)
{
    ensureSegment();

    batchPath_ = segmentPath_/indexedName("batch_", batchIndex_);

    if (Pstream::master())
    {
        mkDir(batchPath_);
    }

    UPstream::barrier(UPstream::worldComm);

    batchFields_ = fields;
    fieldWriters_.clear();

    const label datasets = zoneMode_ ? zones_.size() : 1;
    for (label zi = 0; zi < datasets; ++zi)
    {
        const fileName directory = zoneMode_
            ? batchPath_/"cellZones"/zones_[zi].name : batchPath_;
        if (Pstream::master()) mkDir(directory);
        UPstream::barrier(UPstream::worldComm);
        for (const fieldInfo& info : batchFields_)
        {
            std::vector<label> dimensions
                (1, zoneMode_ ? zones_[zi].cells.size() : info.localSize);
            if (info.nComponents > 1) dimensions.push_back(info.nComponents);
            const word key = zoneMode_ ? word(Foam::name(zi) + "_" + info.name) : info.name;
            fieldWriters_.set(key, new numpyDetail::numpyFileWriter
            (
                directory/(info.name + "_proc_" + Foam::name(Pstream::myProcNo()) + ".npy"),
                dimensions, dtype_
            ));
        }
    }

    if (writeTimes_ && Pstream::master())
    {
        timesWriter_.reset
        (
            new numpyDetail::numpyFileWriter
            (
                batchPath_/"times.npy",
                std::vector<label>(),
                numpyDetail::dataType::FLOAT64
            )
        );
    }

    batchCount_ = 0;
    firstBatchTime_ = VGREAT;
    lastBatchTime_ = -VGREAT;
    batchOpen_ = true;
    writeState(false);
}


void Foam::functionObjects::foamToNumpy::endBatch(bool sealed)
{
    if (!batchOpen_)
    {
        return;
    }

    UPstream::barrier(UPstream::worldComm);
    writeState(sealed);
    fieldWriters_.clear();
    timesWriter_.reset();
    batchFields_.clear();
    batchOpen_ = false;
    batchCount_ = 0;

    if (sealed)
    {
        ++batchIndex_;
    }
}


void Foam::functionObjects::foamToNumpy::writeState(bool sealed) const
{
    if (!Pstream::master() || !batchOpen_)
    {
        return;
    }

    const fileName stateFile(batchPath_/"state");
    const fileName temporaryFile(stateFile + ".tmp");

    {
        OFstream os(temporaryFile);

        IOobject::writeBanner(os);
        os  << "FoamFile" << nl
            << token::BEGIN_BLOCK << nl
            << "    version     2.0;" << nl
            << "    format      ascii;" << nl
            << "    class       dictionary;" << nl
            << "    object      state;" << nl
            << token::END_BLOCK << nl;
        IOobject::writeDivider(os) << nl;

        os  << "batch           " << batchIndex_ << token::END_STATEMENT << nl
            << "meshRevision    " << meshRevision_
            << token::END_STATEMENT << nl
            << "count           " << batchCount_ << token::END_STATEMENT << nl
            << "sealed          " << Switch(sealed)
            << token::END_STATEMENT << nl;

        if (batchCount_)
        {
            os  << "firstTime       " << firstBatchTime_
                << token::END_STATEMENT << nl
                << "lastTime        " << lastBatchTime_
                << token::END_STATEMENT << nl;
        }

        wordList fieldNames(batchFields_.size());
        forAll(fieldNames, i)
        {
            fieldNames[i] = batchFields_[i].name;
        }

        os  << "fields          " << fieldNames
            << token::END_STATEMENT << nl;

        os << "fieldClasses" << nl << token::BEGIN_BLOCK << nl;
        for (const fieldInfo& info : batchFields_)
        {
            os  << "    " << info.name << ' ' << info.className
                << token::END_STATEMENT << nl;
        }
        os << token::END_BLOCK << nl;
        if (zoneMode_)
        {
            wordList names;
            for (const auto& zone : zones_) names.push_back(zone.name);
            os << "zoneLayoutVersion 1;" << nl << "cellZones " << names << ";" << nl;
        }
    }

    if (!mv(temporaryFile, stateFile))
    {
        FatalErrorInFunction
            << "Cannot commit batch state " << stateFile
            << exit(FatalError);
    }
}


void Foam::functionObjects::foamToNumpy::writeSegmentInfo() const
{
    if (!Pstream::master())
    {
        return;
    }

    OFstream os(segmentPath_/"segmentInfo");

    IOobject::writeBanner(os);
    os  << "FoamFile" << nl
        << token::BEGIN_BLOCK << nl
        << "    version     2.0;" << nl
        << "    format      ascii;" << nl
        << "    class       dictionary;" << nl
        << "    object      segmentInfo;" << nl
        << token::END_BLOCK << nl;
    IOobject::writeDivider(os) << nl;

    if (zoneMode_)
    {
        os << "zoneLayoutVersion 1;" << nl
           << "region " << mesh_.name() << ";" << nl
           << "nProcs " << Pstream::nProcs() << ";" << nl;
    }
    os  << "startTime       " << time_.startTime().value()
        << token::END_STATEMENT << nl
        << "firstOutput     " << time_.value() << token::END_STATEMENT << nl
        << "dataType        " << numpyDetail::dataTypeName(dtype_)
        << token::END_STATEMENT << nl
        << "order           fortran" << token::END_STATEMENT << nl
        << "batchSize       " << batchSize_ << token::END_STATEMENT << nl
        << "nativeNumpyCore " << Switch(numpyDetail::hasNativeCore())
        << token::END_STATEMENT << nl;
}


void Foam::functionObjects::foamToNumpy::writeGeometry()
{
    if (!geometryDirty_ || (!zoneMode_ && !writeCellCentres_ && !writeCellVolumes_))
    {
        geometryDirty_ = false;
        return;
    }

    ensureSegment();
    const fileName geometryPath
    (
        segmentPath_/indexedName("geometry_", meshRevision_)
    );

    if (Pstream::master())
    {
        mkDir(geometryPath);
    }

    UPstream::barrier(UPstream::worldComm);

    const word procSuffix
    (
        "_proc_" + Foam::name(Pstream::myProcNo()) + ".npy"
    );

    if (zoneMode_)
    {
        // Geometry construction can communicate across processor patches.
        // Every rank must initialise it, even if all its selected zones are empty.
        const vectorField* centres = writeCellCentres_
            ? &mesh_.C().primitiveField() : nullptr;
        const scalarField* volumes = writeCellVolumes_ ? &mesh_.V() : nullptr;
        OFstream metadata(geometryPath/("mesh_proc_" + Foam::name(Pstream::myProcNo())));
        metadata << "nCells " << mesh_.nCells() << ";" << nl
            << "meshRevision " << meshRevision_ << ";" << nl;
        for (const auto& zone : zones_)
        {
            const fileName directory(geometryPath/"cellZones"/zone.name);
            if (Pstream::master()) mkDir(directory);
            UPstream::barrier(UPstream::worldComm);
            numpyDetail::writeCellIds(directory/("cellIds" + procSuffix), zone.cells);
            if (writeCellCentres_)
            {
                vectorField values(zone.cells.size());
                forAll(values, i) values[i] = (*centres)[zone.cells[i]];
                numpyDetail::numpyFileWriter writer
                    (directory/("cellCentres" + procSuffix), {values.size(), 3}, dtype_);
                writer.appendField(values);
            }
            if (writeCellVolumes_)
            {
                scalarField values(zone.cells.size());
                forAll(values, i) values[i] = (*volumes)[zone.cells[i]];
                numpyDetail::numpyFileWriter writer
                    (directory/("cellVolumes" + procSuffix), {values.size()}, dtype_);
                writer.appendField(values);
            }
        }
        geometryDirty_ = false;
        return;
    }
    if (writeCellCentres_)
    {
        const vectorField& centres = mesh_.C().primitiveField();
        numpyDetail::numpyFileWriter writer
        (
            geometryPath/("cellCentres" + procSuffix),
            {centres.size(), vector::nComponents},
            dtype_
        );
        writer.appendField(centres);
    }

    if (writeCellVolumes_)
    {
        const scalarField& volumes = mesh_.V();
        numpyDetail::numpyFileWriter writer
        (
            geometryPath/("cellVolumes" + procSuffix),
            {volumes.size()},
            dtype_
        );
        writer.appendField(volumes);
    }

    geometryDirty_ = false;
}


void Foam::functionObjects::foamToNumpy::meshChanged()
{
    endBatch(true);
    ++meshRevision_;
    geometryDirty_ = true;
}


Foam::word Foam::functionObjects::foamToNumpy::indexedName
(
    const char* prefix,
    label index
)
{
    std::ostringstream os;
    os << prefix << std::setw(6) << std::setfill('0') << index;
    return word(os.str());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

bool Foam::functionObjects::foamToNumpy::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (dict.found("storageOrder"))
    {
        FatalIOErrorInFunction(dict)
            << "The storageOrder entry is not supported. "
            << "foamToNumpy always writes Fortran-order arrays."
            << exit(FatalIOError);
    }

    const bool newZoneMode = dict.found("cellZones");
    wordRes newZones;
    if (newZoneMode)
    {
        dict.readEntry("cellZones", newZones);
        if (newZones.empty())
            FatalIOErrorInFunction(dict) << "cellZones must not be empty" << exit(FatalIOError);
    }
    if (segmentReady_ && (newZoneMode != zoneMode_ || newZones != zoneSelection_))
        FatalIOErrorInFunction(dict) << "cellZones cannot change after the first snapshot"
            << exit(FatalIOError);
    zoneMode_ = newZoneMode;
    zoneSelection_ = newZones;
    wordRes newFieldSelection;
    dict.readEntry("fields", newFieldSelection);
    newFieldSelection.uniq();

    const numpyDetail::dataType newDtype = numpyDetail::parseDataType
    (
        dict.getOrDefault<word>("dataType", "float64")
    );

    const label newBatchSize = dict.getOrDefault<label>("batchSize", 100);
    if (newBatchSize <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "batchSize must be greater than zero"
            << exit(FatalIOError);
    }

    const bool newWriteTimes = dict.getOrDefault("writeTimes", true);
    const bool newWriteCellCentres =
        dict.getOrDefault("writeCellCentres", false);
    const bool newWriteCellVolumes =
        dict.getOrDefault("writeCellVolumes", false);

    fileName configuredPath;
    fileName newOutputPath;
    if (dict.readIfPresent("outputDir", configuredPath))
    {
        newOutputPath =
        (
            configuredPath.isAbsolute()
          ? configuredPath
          : time_.globalPath()/configuredPath
        );
    }
    else
    {
        newOutputPath = time_.globalPath()/functionObject::outputPrefix;

        if (mesh_.name() != polyMesh::defaultRegion)
        {
            newOutputPath /= mesh_.name();
        }

        newOutputPath /= name();
    }

    newOutputPath.clean();

    if
    (
        segmentReady_
     &&
        (
            newFieldSelection != fieldSelection_
         || newDtype != dtype_
         || newBatchSize != batchSize_
         || newWriteTimes != writeTimes_
         || newWriteCellCentres != writeCellCentres_
         || newWriteCellVolumes != writeCellVolumes_
         || newOutputPath != outputPath_
        )
    )
    {
        FatalIOErrorInFunction(dict)
            << "Output settings cannot be changed after the first snapshot"
            << exit(FatalIOError);
    }

    fieldSelection_ = newFieldSelection;
    dtype_ = newDtype;
    batchSize_ = newBatchSize;
    writeTimes_ = newWriteTimes;
    writeCellCentres_ = newWriteCellCentres;
    writeCellVolumes_ = newWriteCellVolumes;
    outputPath_ = newOutputPath;
    return true;
}


bool Foam::functionObjects::foamToNumpy::execute()
{
    return true;
}


bool Foam::functionObjects::foamToNumpy::write()
{
    const scalar outputTime = time_.value();
    scalar minTime = outputTime;
    scalar maxTime = outputTime;
    reduce(minTime, minOp<scalar>());
    reduce(maxTime, maxOp<scalar>());

    if (mag(maxTime - minTime) > SMALL)
    {
        FatalErrorInFunction
            << "Inconsistent output times across processors: "
            << minTime << " and " << maxTime
            << exit(FatalError);
    }

    if (mag(outputTime - lastOutputTime_) <= SMALL)
    {
        WarningInFunction
            << "Skipping duplicate output at time " << time_.timeName()
            << endl;
        return true;
    }

    if (outputTime < lastOutputTime_)
    {
        FatalErrorInFunction
            << "Output time moved backwards from " << lastOutputTime_
            << " to " << outputTime
            << exit(FatalError);
    }

    if (zoneMode_)
    {
        auto mappings = numpyDetail::cellZoneMappings(mesh_, zoneSelection_);
        bool changed = mappings != zones_;
        reduce(changed, orOp<bool>());
        if (changed && segmentReady_) meshChanged();
        zones_ = std::move(mappings);
    }
    const std::vector<fieldInfo> fields(selectedFields());

    if (!batchOpen_)
    {
        beginBatch(fields);
    }
    else
    {
        if (fields.size() != batchFields_.size())
        {
            FatalErrorInFunction
                << "Selected field set changed within batch " << batchIndex_
                << exit(FatalError);
        }

        for (std::size_t i = 0; i < fields.size(); ++i)
        {
            if
            (
                fields[i].name != batchFields_[i].name
             || fields[i].className != batchFields_[i].className
             || fields[i].localSize != batchFields_[i].localSize
            )
            {
                FatalErrorInFunction
                    << "Field schema changed for " << fields[i].name
                    << " within batch " << batchIndex_
                    << exit(FatalError);
            }
        }
    }

    writeGeometry();

    for (const fieldInfo& info : batchFields_)
    {
        bool written = false;

        if (info.className == volScalarField::typeName)
        {
            written = appendField<volScalarField>(info.name);
        }
        else if (info.className == volVectorField::typeName)
        {
            written = appendField<volVectorField>(info.name);
        }
        else if (info.className == volSphericalTensorField::typeName)
        {
            written = appendField<volSphericalTensorField>(info.name);
        }
        else if (info.className == volSymmTensorField::typeName)
        {
            written = appendField<volSymmTensorField>(info.name);
        }
        else if (info.className == volTensorField::typeName)
        {
            written = appendField<volTensorField>(info.name);
        }

        if (!written)
        {
            FatalErrorInFunction
                << "Cannot access selected volume field " << info.name
                << exit(FatalError);
        }
    }

    UPstream::barrier(UPstream::worldComm);

    if (timesWriter_)
    {
        timesWriter_->appendValue(outputTime);
    }

    UPstream::barrier(UPstream::worldComm);

    if (!batchCount_)
    {
        firstBatchTime_ = outputTime;
    }

    ++batchCount_;
    lastBatchTime_ = outputTime;
    lastOutputTime_ = outputTime;
    writeState(false);

    Log << type() << ' ' << name() << ": committed time "
        << time_.timeName() << " to batch " << batchIndex_
        << " (" << batchCount_ << '/' << batchSize_ << ')' << endl;

    if (batchCount_ == batchSize_)
    {
        endBatch(true);
    }

    return true;
}


bool Foam::functionObjects::foamToNumpy::end()
{
    endBatch(true);
    return true;
}


void Foam::functionObjects::foamToNumpy::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        meshChanged();
    }
}


void Foam::functionObjects::foamToNumpy::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        meshChanged();
    }
}


void Foam::functionObjects::foamToNumpy::readUpdate
(
    const polyMesh::readUpdateState state
)
{
    if (state != polyMesh::UNCHANGED)
    {
        meshChanged();
    }
}


// ************************************************************************* //
