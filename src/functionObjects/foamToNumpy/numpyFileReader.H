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
    Foam::functionObjects::numpyDetail::numpyFileReader

Description
    Reader for little-endian float32/float64 fields and int64 cell mappings.

    The reader intentionally supports only the representation emitted by the
    NumPy function objects. In particular, C-order arrays are rejected.

SourceFiles
    numpyFileReader.C
    numpyFileReader.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_numpyFileReader_H
#define Foam_functionObjects_numpyFileReader_H

#include "Field.H"
#include "fileName.H"
#include "foamNumpyCoreAdapter.H"

#include <cstddef>
#include <fstream>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace numpyDetail
{

/*---------------------------------------------------------------------------*\
                       Class numpyFileReader Declaration
\*---------------------------------------------------------------------------*/

class numpyFileReader
{
    // Private Data

        //- Input file
        fileName path_;

        //- Numeric data type
        dataType dtype_;
        //- True only for dedicated int64 addressing arrays
        bool integer_;

        //- Array shape
        std::vector<label> shape_;

        //- Byte offset of the binary payload
        std::streamoff dataStart_;


    // Private Member Functions

        //- Read a scalar payload value
        scalar readValue(std::ifstream& is) const;

        //- No copy construct
        numpyFileReader(const numpyFileReader&) = delete;

        //- No copy assignment
        void operator=(const numpyFileReader&) = delete;


public:

    // Constructors

        //- Read and validate a NumPy header
        explicit numpyFileReader(const fileName& path);


    // Member Functions

        //- Input path
        const fileName& path() const noexcept
        {
            return path_;
        }

        //- Array shape
        const std::vector<label>& shape() const noexcept
        {
            return shape_;
        }

        //- Number of snapshots in the final dimension
        label nSnapshots() const noexcept
        {
            return shape_.empty() ? 0 : shape_.back();
        }

        //- Read a one-dimensional scalar array
        scalarField readScalarArray() const;

        labelList readCellIds() const;

        //- Read one field snapshot
        template<class Type>
        void readField
        (
            label snapshoti,
            label expectedSize,
            Field<Type>& values
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace numpyDetail
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "numpyFileReader.txx"

#endif

// ************************************************************************* //
