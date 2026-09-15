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

Class
    Foam::functionObjects::foamToNumpy

Group
    grpFieldFunctionObjects

Description
    Writes internal volume fields to appendable, Fortran-order NumPy files.

    Output is divided into fixed-size batches. Each solver or postProcess
    invocation uses a start-time segment, which keeps restarted trajectories
    separate from data produced before a restart.

Usage
    Minimal example using \c system/controlDict.functions:
    \verbatim
    numpyExport
    {
        type                foamToNumpy;
        libs                (numpyFunctionObjects);

        fields              (p U);
        dataType            float64;
        batchSize           100;

        writeTimes          true;
        writeCellCentres    false;
        writeCellVolumes    false;

        writeControl        adjustableRunTime;
        writeInterval       0.1;
    }
    \endverbatim

    The inherited function-object controls determine when snapshots are
    written. In post-processing mode, the selected postProcess times override
    these controls as usual.

    The entries are:
    \table
      Property          | Description                         | Type | Req'd | Dflt
      fields            | Volume field names or expressions   | wordRes | yes | -
      cellZones         | Optional cell-zone names/expressions | wordRes | no | -
      dataType          | NumPy precision: float32/float64    | word | no | float64
      batchSize         | Maximum snapshots in a batch        | label | no | 100
      writeTimes        | Write batch time arrays             | bool | no | true
      writeCellCentres  | Write cell centres by mesh revision | bool | no | false
      writeCellVolumes  | Write cell volumes by mesh revision | bool | no | false
      outputDir         | Alternative output directory        | fileName | no | -
    \endtable

    Only Fortran-order output is supported. Boundary fields are not written.

SourceFiles
    foamToNumpy.C
    foamToNumpy.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_foamToNumpy_H
#define Foam_functionObjects_foamToNumpy_H

#include "HashPtrTable.H"
#include "fvMeshFunctionObject.H"
#include "numpyFileWriter.H"
#include "polyMesh.H"
#include "wordRes.H"
#include "numpyCellZones.H"

#include <memory>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class foamToNumpy Declaration
\*---------------------------------------------------------------------------*/

class foamToNumpy
:
    public fvMeshFunctionObject
{
    // Private Data Types

        struct fieldInfo
        {
            word name;
            word className;
            label localSize;
            direction nComponents;
        };


    // Private Data

        //- Selected volume fields
        wordRes fieldSelection_;
        wordRes zoneSelection_;
        std::vector<numpyDetail::cellZoneAddressing> zones_;
        bool zoneMode_ = false;

        //- Output data type
        numpyDetail::dataType dtype_;

        //- Maximum outputs per batch
        label batchSize_;

        //- Write the batch times array
        bool writeTimes_;

        //- Write cell centres for each mesh revision
        bool writeCellCentres_;

        //- Write cell volumes for each mesh revision
        bool writeCellVolumes_;

        //- Root output directory
        fileName outputPath_;

        //- Current invocation segment directory
        fileName segmentPath_;

        //- Current batch directory
        fileName batchPath_;

        //- Has the invocation segment been created
        bool segmentReady_;

        //- Is a batch currently open
        bool batchOpen_;

        //- Geometry must be written before the next snapshot
        bool geometryDirty_;

        //- Batch index within the invocation segment
        label batchIndex_;

        //- Committed outputs in the current batch
        label batchCount_;

        //- Mesh revision within the invocation segment
        label meshRevision_;

        //- First and last committed batch times
        scalar firstBatchTime_;
        scalar lastBatchTime_;

        //- Last output time in this invocation
        scalar lastOutputTime_;

        //- Schema for the open batch
        std::vector<fieldInfo> batchFields_;

        //- Per-field, processor-local NumPy writers
        HashPtrTable<numpyDetail::numpyFileWriter> fieldWriters_;

        //- Master-only batch time writer
        std::unique_ptr<numpyDetail::numpyFileWriter> timesWriter_;


    // Private Member Functions

        //- Add selected fields of a given volume-field type
        template<class GeoField>
        void appendFieldInfo
        (
            std::vector<fieldInfo>& info,
            HashSet<word>& selected
        ) const;

        //- Discover and validate selected fields
        std::vector<fieldInfo> selectedFields() const;

        //- Append a field if it has the requested geometric type
        template<class GeoField>
        bool appendField(const word& fieldName);

        template<class Type>
        void appendValues(const word&, const Field<Type>&);

        //- Create a unique start-time segment for this invocation
        void ensureSegment();

        //- Create writers for a new batch
        void beginBatch(const std::vector<fieldInfo>& fields);

        //- Close the current batch and optionally mark it sealed
        void endBatch(bool sealed);

        //- Write the authoritative commit state for the current batch
        void writeState(bool sealed) const;

        //- Write invocation metadata
        void writeSegmentInfo() const;

        //- Write geometry for the current mesh revision
        void writeGeometry();

        //- Finish a batch and invalidate geometry for a mesh change
        void meshChanged();

        //- Return a zero-padded batch or revision name
        static word indexedName(const char* prefix, label index);

        //- No copy construct
        foamToNumpy(const foamToNumpy&) = delete;

        //- No copy assignment
        void operator=(const foamToNumpy&) = delete;


public:

    //- Runtime type information
    TypeName("foamToNumpy");


    // Constructors

        //- Construct from Time and dictionary
        foamToNumpy
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~foamToNumpy() = default;


    // Member Functions

        //- Read the function-object dictionary
        virtual bool read(const dictionary& dict);

        //- Execute is a no-op; snapshots are captured by write()
        virtual bool execute();

        //- Append the current snapshot
        virtual bool write();

        //- Seal a partially filled final batch
        virtual bool end();

        //- Handle mesh topology changes
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Handle mesh point motion
        virtual void movePoints(const polyMesh& mesh);

        //- Handle mesh changes read from disk
        virtual void readUpdate(const polyMesh::readUpdateState state);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamToNumpy.txx"

#endif

// ************************************************************************* //
