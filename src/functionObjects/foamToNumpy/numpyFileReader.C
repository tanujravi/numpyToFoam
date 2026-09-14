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

#include "numpyFileReader.H"

#include "error.H"

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

std::string dictionaryValue(const std::string& header, const char* key)
{
    const std::size_t keyPos = header.find(key);

    if (keyPos == std::string::npos)
    {
        return std::string();
    }

    const std::size_t colon = header.find(':', keyPos);
    if (colon == std::string::npos)
    {
        return std::string();
    }

    const std::size_t begin = header.find_first_not_of(" \t'\"", colon + 1);
    const std::size_t end = header.find_first_of(" ,}'\"", begin);

    return header.substr(begin, end - begin);
}


std::vector<Foam::label> parseShape(const std::string& header)
{
    const std::size_t keyPos = header.find("shape");
    const std::size_t begin = header.find('(', keyPos);
    const std::size_t end = header.find(')', begin);

    if
    (
        keyPos == std::string::npos
     || begin == std::string::npos
     || end == std::string::npos
    )
    {
        return std::vector<Foam::label>();
    }

    std::vector<Foam::label> shape;
    std::stringstream entries(header.substr(begin + 1, end - begin - 1));
    std::string entry;

    while (std::getline(entries, entry, ','))
    {
        const std::size_t first = entry.find_first_not_of(" \t");

        if (first != std::string::npos)
        {
            const auto value = std::stoll(entry.substr(first));

            if (value < 0 || value > Foam::labelMax)
            {
                return std::vector<Foam::label>();
            }

            shape.push_back(static_cast<Foam::label>(value));
        }
    }

    return shape;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::functionObjects::numpyDetail::numpyFileReader::numpyFileReader
(
    const fileName& path
)
:
    path_(path),
    dtype_(dataType::FLOAT64),
    integer_(false),
    shape_(),
    dataStart_(0)
{
    std::ifstream is(path_.c_str(), std::ios::binary);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open NumPy file " << path_
            << exit(FatalError);
    }

    char magic[6] = {};
    is.read(magic, sizeof(magic));

    static const char expectedMagic[6] =
        {'\x93', 'N', 'U', 'M', 'P', 'Y'};

    if (!is.good() || !std::equal(magic, magic + 6, expectedMagic))
    {
        FatalErrorInFunction
            << path_ << " is not a NumPy file"
            << exit(FatalError);
    }

    const unsigned char major = static_cast<unsigned char>(is.get());
    const unsigned char minor = static_cast<unsigned char>(is.get());
    (void)minor;

    std::uint32_t headerLength = 0;

    if (major == 1)
    {
        const std::uint32_t b0 = static_cast<unsigned char>(is.get());
        const std::uint32_t b1 = static_cast<unsigned char>(is.get());
        headerLength = b0 | (b1 << 8);
    }
    else if (major == 2 || major == 3)
    {
        for (unsigned int bytei = 0; bytei < 4; ++bytei)
        {
            headerLength |=
                static_cast<std::uint32_t>
                (
                    static_cast<unsigned char>(is.get())
                ) << (8*bytei);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported NumPy format version " << int(major)
            << " in " << path_
            << exit(FatalError);
    }

    std::string header(headerLength, ' ');
    is.read(&header[0], headerLength);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read NumPy header from " << path_
            << exit(FatalError);
    }

    const std::string dtype(dictionaryValue(header, "descr"));
    if (dtype == "<f8" || dtype == "=f8")
    {
        dtype_ = dataType::FLOAT64;
    }
    else if (dtype == "<f4" || dtype == "=f4")
    {
        dtype_ = dataType::FLOAT32;
    }
    else if (dtype == "<i8" || dtype == "=i8")
    {
        integer_ = true;
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported NumPy dtype '" << dtype << "' in " << path_
            << ". Expected little-endian float32, float64, or int64 addressing."
            << exit(FatalError);
    }

    const std::string fortranOrder
    (
        dictionaryValue(header, "fortran_order")
    );

    if (fortranOrder != "True")
    {
        FatalErrorInFunction
            << "NumPy file " << path_ << " is not Fortran-order. "
            << "C-order input is not supported."
            << exit(FatalError);
    }

    shape_ = parseShape(header);
    if (shape_.empty())
    {
        FatalErrorInFunction
            << "Cannot parse NumPy shape in " << path_
            << exit(FatalError);
    }

    dataStart_ = is.tellg();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::scalar
Foam::functionObjects::numpyDetail::numpyFileReader::readValue
(
    std::ifstream& is
) const
{
    if (dtype_ == dataType::FLOAT64)
    {
        double value = 0;
        is.read(reinterpret_cast<char*>(&value), sizeof(value));
        return value;
    }

    float value = 0;
    is.read(reinterpret_cast<char*>(&value), sizeof(value));
    return value;
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalarField
Foam::functionObjects::numpyDetail::numpyFileReader::readScalarArray() const
{
    if (integer_ || shape_.size() != 1)
    {
        FatalErrorInFunction
            << "Expected a one-dimensional scalar array in " << path_
            << exit(FatalError);
    }

    scalarField values(shape_[0]);
    std::ifstream is(path_.c_str(), std::ios::binary);
    is.seekg(dataStart_, std::ios::beg);

    for (scalar& value : values)
    {
        value = readValue(is);
    }

    if (!is.good())
    {
        FatalErrorInFunction
            << "Unexpected end of NumPy payload in " << path_
            << exit(FatalError);
    }

    return values;
}


Foam::labelList
Foam::functionObjects::numpyDetail::numpyFileReader::readCellIds() const
{
    if (!integer_ || shape_.size() != 1)
    {
        FatalErrorInFunction << "Expected one-dimensional int64 cell IDs in "
            << path_ << exit(FatalError);
    }
    labelList ids(shape_[0]);
    std::ifstream is(path_.c_str(), std::ios::binary);
    is.seekg(dataStart_);
    for (label& id : ids)
    {
        // Decode little-endian explicitly, independent of host byte order.
        std::uint64_t value = 0;
        for (unsigned b = 0; b < 8; ++b)
            value |= std::uint64_t(static_cast<unsigned char>(is.get())) << (8*b);
        if (!is.good() || value > std::uint64_t(labelMax))
        {
            FatalErrorInFunction << "Invalid or truncated cell ID in " << path_
                << exit(FatalError);
        }
        id = label(value);
    }
    if (is.peek() != std::ifstream::traits_type::eof())
    {
        FatalErrorInFunction << "Trailing cell ID payload in " << path_
            << exit(FatalError);
    }
    return ids;
}

// ************************************************************************* //
