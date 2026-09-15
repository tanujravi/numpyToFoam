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

Class
    Foam::functionObjects::numpyDetail::numpyInputCatalog

Description
    Index of committed snapshots in foamToNumpy segment and batch output.

SourceFiles
    numpyInputCatalog.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_numpyInputCatalog_H
#define Foam_functionObjects_numpyInputCatalog_H

#include "dictionary.H"
#include "fileName.H"
#include "scalarList.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace numpyDetail
{

/*---------------------------------------------------------------------------*\
                    Class numpyInputCatalog Declaration
\*---------------------------------------------------------------------------*/

class numpyInputCatalog
{
public:

    // Public Data Types

        struct snapshot
        {
            scalar time;
            fileName batchPath;
            label index;
            label meshRevision;
        };


private:

    // Private Data

        //- Root input directory
        fileName inputPath_;

        //- Snapshots sorted by time
        std::vector<snapshot> snapshots_;


    // Private Member Functions

        //- Append committed snapshots from one segment
        void appendSegment(const fileName& segmentPath);

        //- Sort by time and let later segments replace duplicate times
        void normalise();


public:

    // Constructors

        //- Construct from an input dictionary
        numpyInputCatalog
        (
            const fileName& casePath,
            const dictionary& dict
        );


    // Member Functions

        //- Root input directory
        const fileName& inputPath() const noexcept
        {
            return inputPath_;
        }

        //- Available time values
        scalarList times() const;

        //- Find a snapshot matching the time value
        const snapshot& find(scalar timeValue) const;

        //- Path for a processor-local field in a snapshot
        fileName fieldPath
        (
            const snapshot& sample,
            const word& fieldName,
            label proci,
            const word& zone = word::null
        ) const;

        fileName geometryPath(const snapshot& sample) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numpyDetail
} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
