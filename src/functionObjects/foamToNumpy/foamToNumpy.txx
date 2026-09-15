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

#include "fvMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
void Foam::functionObjects::foamToNumpy::appendFieldInfo
(
    std::vector<fieldInfo>& info,
    HashSet<word>& selected
) const
{
    HashSet<word> candidates;

    for (const word& fieldName : mesh_.sortedNames<GeoField>(fieldSelection_))
    {
        candidates.insert(fieldName);
    }

    if (functionObject::postProcess)
    {
        const IOobjectList objects(mesh_, time_.timeName());

        for
        (
            const word& fieldName
          : objects.names<GeoField>(fieldSelection_, true)
        )
        {
            candidates.insert(fieldName);
        }
    }

    for (const word& fieldName : candidates.sortedToc())
    {
        if (selected.insert(fieldName))
        {
            const GeoField* fieldPtr = mesh_.cfindObject<GeoField>(fieldName);

            if (fieldPtr)
            {
                info.push_back
                ({
                    fieldName,
                    GeoField::typeName,
                    fieldPtr->size(),
                    pTraits<typename GeoField::value_type>::nComponents
                });
            }
            else
            {
                const GeoField field
                (
                    IOobject
                    (
                        fieldName,
                        time_.timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

                info.push_back
                ({
                    fieldName,
                    GeoField::typeName,
                    field.size(),
                    pTraits<typename GeoField::value_type>::nComponents
                });
            }
        }
    }
}


template<class GeoField>
bool Foam::functionObjects::foamToNumpy::appendField
(
    const word& fieldName
)
{
    const GeoField* fieldPtr = mesh_.cfindObject<GeoField>(fieldName);

    if (!fieldPtr)
    {
        if (!functionObject::postProcess)
        {
            return false;
        }

        IOobject io
        (
            fieldName,
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<GeoField>(false))
        {
            return false;
        }

        const GeoField field(io, mesh_);
        appendValues(fieldName, field.primitiveField());
        return true;
    }

    appendValues(fieldName, fieldPtr->primitiveField());
    return true;
}


template<class Type>
void Foam::functionObjects::foamToNumpy::appendValues
(const word& name, const Field<Type>& values)
{
    if (!zoneMode_)
    {
        fieldWriters_[name]->appendField(values);
        return;
    }
    label zi = 0;
    for (const auto& zone : zones_)
    {
        Field<Type> subset(zone.cells.size());
        forAll(subset, i) subset[i] = values[zone.cells[i]];
        fieldWriters_[Foam::name(zi++) + "_" + name]->appendField(subset);
    }
}

// ************************************************************************* //
