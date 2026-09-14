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

#include "error.H"
#include "pTraits.H"

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::numpyDetail::numpyFileReader::readField
(
    const label snapshoti,
    const label expectedSize,
    Field<Type>& values
) const
{
    if (integer_)
    {
        FatalErrorInFunction << "Expected floating-point field in " << path_
            << exit(FatalError);
    }
    const direction nComponents = pTraits<Type>::nComponents;
    const bool scalarType = (nComponents == 1);
    const std::size_t expectedRank = scalarType ? 2u : 3u;

    if (shape_.size() != expectedRank)
    {
        FatalErrorInFunction
            << "Expected rank " << expectedRank << " for "
            << pTraits<Type>::typeName << " field in " << path_ << nl
            << "Found shape with rank " << shape_.size()
            << exit(FatalError);
    }

    if (shape_[0] != expectedSize)
    {
        FatalErrorInFunction
            << "Entity count mismatch in " << path_ << nl
            << "Expected " << expectedSize << " but found " << shape_[0]
            << exit(FatalError);
    }

    if (!scalarType && shape_[1] != nComponents)
    {
        FatalErrorInFunction
            << "Component count mismatch in " << path_ << nl
            << "Expected " << nComponents << " but found " << shape_[1]
            << exit(FatalError);
    }

    if (snapshoti < 0 || snapshoti >= nSnapshots())
    {
        FatalErrorInFunction
            << "Snapshot index " << snapshoti << " is outside [0, "
            << nSnapshots() << ") in " << path_
            << exit(FatalError);
    }

    const std::size_t valuesPerSnapshot =
        static_cast<std::size_t>(expectedSize)*nComponents;
    const std::streamoff byteOffset =
        dataStart_
      + static_cast<std::streamoff>(snapshoti*valuesPerSnapshot)
       *static_cast<std::streamoff>(dataTypeSize(dtype_));

    std::ifstream is(path_.c_str(), std::ios::binary);
    is.seekg(byteOffset, std::ios::beg);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot seek to snapshot " << snapshoti << " in " << path_
            << exit(FatalError);
    }

    values.resize(expectedSize);
    values = Zero;

    for (direction cmpti = 0; cmpti < nComponents; ++cmpti)
    {
        forAll(values, entityi)
        {
            setComponent(values[entityi], cmpti) = readValue(is);
        }
    }

    if (!is.good())
    {
        FatalErrorInFunction
            << "Unexpected end of NumPy payload in " << path_
            << exit(FatalError);
    }
}


// ************************************************************************* //
