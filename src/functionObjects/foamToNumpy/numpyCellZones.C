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

#include "numpyCellZones.H"
#include "cellZoneMesh.H"
#include "Pstream.H"
#include "foamNumpyCoreCompat.H"
#include <fstream>
#include <cstdint>
#include <iomanip>
#include <sstream>

namespace Foam
{
namespace functionObjects
{
namespace numpyDetail
{

wordList selectCellZoneNames(const wordRes& selection, const wordList& available)
{
    if (selection.empty())
    {
        FatalErrorInFunction
            << "cellZones must not be empty" << exit(FatalError);
    }

    wordList result;
    for (const auto& pattern : selection)
    {
        bool matched = false;
        for (const word& name : available)
        {
            if (pattern.match(name))
            {
                result.push_back(name);
                matched = true;
            }
        }
        if (!matched)
        {
            FatalErrorInFunction
                << "Unmatched cellZones selection " << pattern
                << "; available zones: " << available << exit(FatalError);
        }
    }
    inplaceUniqueSort(result);
    return result;
}


std::vector<cellZoneAddressing> cellZoneMappings
(const fvMesh& mesh, const wordRes& selection)
{
    const wordList names(selectCellZoneNames(selection, mesh.cellZones().names()));
    wordList masterNames(names);
    Pstream::broadcast(masterNames);
    if (masterNames != names)
    {
        FatalErrorInFunction
            << "Inconsistent cell zones across processors" << exit(FatalError);
    }
    std::vector<cellZoneAddressing> result;
    for (const word& name : names)
    {
        labelList cells(mesh.cellZones()[mesh.cellZones().findZoneID(name)]);
        Foam::sort(cells);
        forAll(cells, i)
        {
            if (cells[i] < 0 || cells[i] >= mesh.nCells()
                || (i && cells[i] == cells[i-1]))
            {
                FatalErrorInFunction
                    << "Invalid cell addressing in zone " << name
                    << exit(FatalError);
            }
        }
        result.push_back({name, std::move(cells)});
    }
    return result;
}


void writeCellIds(const fileName& path, const labelList& cells)
{
    std::ofstream os(path.c_str(), std::ios::binary);
    foamNumpyCoreCompat::writeHeader
        (os, "<i8", "(" + std::to_string(cells.size()) + ",)", true);
    for (const label cell : cells)
    {
        const std::uint64_t value(cell);
        for (unsigned b = 0; b < 8; ++b)
        {
            os.put(char((value >> (8*b)) & 255));
        }
    }
    os.flush();
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write cell IDs " << path << exit(FatalError);
    }
}


word revisionName(label revision)
{
    std::ostringstream os;
    os << "geometry_" << std::setw(6) << std::setfill('0') << revision;
    return word(os.str());
}

} // End namespace numpyDetail
} // End namespace functionObjects
} // End namespace Foam
