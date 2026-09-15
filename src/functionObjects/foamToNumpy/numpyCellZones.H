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

// Shared processor-local addressing for volume cell-zone exchange.
#ifndef numpyCellZones_H
#define numpyCellZones_H
#include "fvMesh.H"
#include "wordRes.H"
#include "labelList.H"
#include <vector>
namespace Foam
{
namespace functionObjects
{
namespace numpyDetail
{

struct cellZoneAddressing
{
    word name;
    labelList cells;
    bool operator==(const cellZoneAddressing& other) const
    {
        return name == other.name && cells == other.cells;
    }
};
//- Resolve every selection and return sorted unique zone names.
wordList selectCellZoneNames(const wordRes&, const wordList&);
//- Validate and sort processor-local addressing, including empty zones.
std::vector<cellZoneAddressing> cellZoneMappings(const fvMesh&, const wordRes&);
//- Write a one-dimensional Fortran-order little-endian int64 array.
void writeCellIds(const fileName&, const labelList&);
word revisionName(label);

} // End namespace numpyDetail
} // End namespace functionObjects
} // End namespace Foam
#endif
