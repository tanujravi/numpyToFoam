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

#include "numpyInputCatalog.H"

#include "IFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "numpyFileReader.H"
#include "numpyCellZones.H"

#include <algorithm>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::functionObjects::numpyDetail::numpyInputCatalog::appendSegment
(
    const fileName& segmentPath
)
{
    const fileNameList entries(readDir(segmentPath, fileName::DIRECTORY));

    for (const fileName& entry : entries)
    {
        if (!entry.starts_with("batch_"))
        {
            continue;
        }

        const fileName batchPath(segmentPath/entry);
        const fileName statePath(batchPath/"state");
        const fileName timesPath(batchPath/"times.npy");

        if (!isFile(statePath) || !isFile(timesPath))
        {
            continue;
        }

        IFstream stateStream(statePath);
        dictionary state(stateStream);
        const label count = state.get<label>("count");
        const label meshRevision =
            state.getOrDefault<label>("meshRevision", 0);

        const numpyFileReader timesReader(timesPath);
        const scalarField batchTimes(timesReader.readScalarArray());

        if (count < 0 || count > batchTimes.size())
        {
            FatalErrorInFunction
                << "Invalid committed count " << count << " in " << statePath
                << ". The times array contains " << batchTimes.size()
                << " values."
                << exit(FatalError);
        }

        for (label timei = 0; timei < count; ++timei)
        {
            snapshots_.push_back
            ({
                batchTimes[timei],
                batchPath,
                timei,
                meshRevision
            });
        }
    }
}


void Foam::functionObjects::numpyDetail::numpyInputCatalog::normalise()
{
    // Stable sorting preserves segment order for duplicate time values. The
    // final duplicate below therefore belongs to the later configured segment.
    std::stable_sort
    (
        snapshots_.begin(),
        snapshots_.end(),
        [](const snapshot& a, const snapshot& b)
        {
            return a.time < b.time;
        }
    );

    std::vector<snapshot> unique;
    unique.reserve(snapshots_.size());

    for (const snapshot& sample : snapshots_)
    {
        if
        (
            !unique.empty()
         && mag(sample.time - unique.back().time)
          <= SMALL*max(scalar(1), mag(sample.time))
        )
        {
            unique.back() = sample;
        }
        else
        {
            unique.push_back(sample);
        }
    }

    snapshots_.swap(unique);
}


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::functionObjects::numpyDetail::numpyInputCatalog::numpyInputCatalog
(
    const fileName& casePath,
    const dictionary& dict
)
:
    inputPath_(dict.getOrDefault<fileName>("inputDir", "postProcessing")),
    snapshots_()
{
    if (!inputPath_.isAbsolute())
    {
        inputPath_ = casePath/inputPath_;
    }
    inputPath_.clean();

    wordList segments;

    if (!dict.readIfPresent("segments", segments))
    {
        word segment;

        if (dict.readIfPresent("segment", segment))
        {
            segments.resize(1);
            segments.front() = segment;
        }
        else
        {
            const fileNameList entries
            (
                readDir(inputPath_, fileName::DIRECTORY)
            );

            for (const fileName& entry : entries)
            {
                if (isFile(inputPath_/entry/"segmentInfo"))
                {
                    segments.push_back(entry);
                }
            }

            if (segments.size() > 1)
            {
                FatalIOErrorInFunction(dict)
                    << "Multiple NumPy restart segments exist in "
                    << inputPath_ << ": " << segments << nl
                    << "Specify their precedence explicitly with 'segments'. "
                    << "Later entries replace duplicate times."
                    << exit(FatalIOError);
            }
        }
    }

    if (segments.empty())
    {
        FatalIOErrorInFunction(dict)
            << "No NumPy segments found in " << inputPath_
            << exit(FatalIOError);
    }

    for (const word& segment : segments)
    {
        const fileName segmentPath(inputPath_/segment);

        if (!isDir(segmentPath))
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find NumPy segment " << segmentPath
                << exit(FatalIOError);
        }

        appendSegment(segmentPath);
    }

    normalise();

    if (snapshots_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "No committed NumPy snapshots found in " << inputPath_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalarList
Foam::functionObjects::numpyDetail::numpyInputCatalog::times() const
{
    scalarList values(snapshots_.size());

    forAll(values, timei)
    {
        values[timei] = snapshots_[timei].time;
    }

    return values;
}


const Foam::functionObjects::numpyDetail::numpyInputCatalog::snapshot&
Foam::functionObjects::numpyDetail::numpyInputCatalog::find
(
    const scalar timeValue
) const
{
    for (const snapshot& sample : snapshots_)
    {
        if
        (
            mag(sample.time - timeValue)
         <= SMALL*max(scalar(1), mag(timeValue))
        )
        {
            return sample;
        }
    }

    FatalErrorInFunction
        << "No committed NumPy snapshot exists at time " << timeValue << nl
        << "Available times: " << times()
        << exit(FatalError);

    return snapshots_.front();
}


Foam::fileName
Foam::functionObjects::numpyDetail::numpyInputCatalog::fieldPath
(
    const snapshot& sample,
    const word& fieldName,
    const label proci,
    const word& zone
) const
{
    const fileName directory = zone.empty() ? sample.batchPath
        : sample.batchPath/"cellZones"/zone;
    return
        directory
       /(fieldName + "_proc_" + Foam::name(proci) + ".npy");
}


Foam::fileName Foam::functionObjects::numpyDetail::numpyInputCatalog::geometryPath
(const snapshot& sample) const
{
    return sample.batchPath.path()/revisionName(sample.meshRevision);
}

// ************************************************************************* //
