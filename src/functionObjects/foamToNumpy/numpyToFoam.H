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
    Foam::functionObjects::numpyToFoam

Description
    Loads volume-field snapshots written by foamToNumpy and registers the
    fields for use by subsequent function objects.

    Existing registered fields are updated in place. Otherwise, the field is
    constructed from its OpenFOAM template file, which supplies dimensions
    and boundary conditions. No native field output is produced unless
    writeFields is enabled.

SourceFiles
    numpyToFoam.C
    numpyToFoam.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_numpyToFoam_H
#define Foam_functionObjects_numpyToFoam_H

#include "HashPtrTable.H"
#include "fvMeshFunctionObject.H"
#include "numpyInputCatalog.H"
#include "volFieldsFwd.H"
#include "wordList.H"
#include "numpyCellZones.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class numpyToFoam Declaration
\*---------------------------------------------------------------------------*/

class numpyToFoam
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Exact field names to import
        wordList fieldNames_;
        bool zoneMode_ = false;
        wordRes zoneSelection_;
        std::vector<numpyDetail::cellZoneAddressing> zones_;


        //- Instance containing field templates
        word templateInstance_;

        //- Write imported fields in OpenFOAM format
        bool writeFields_;

        //- Correct patch fields after loading internal values
        bool correctBoundaryConditions_;

        //- Input snapshot catalogue
        std::unique_ptr<numpyDetail::numpyInputCatalog> catalog_;

        //- Last imported time
        scalar lastTime_;

        //- Volume fields owned by this function object
        HashPtrTable<volScalarField> scalarFields_;
        HashPtrTable<volVectorField> vectorFields_;
        HashPtrTable<volSphericalTensorField> sphericalTensorFields_;
        HashPtrTable<volSymmTensorField> symmTensorFields_;
        HashPtrTable<volTensorField> tensorFields_;


    // Private Member Functions

        //- Validate source layout and current mesh before assigning values
        void validateZones(const numpyDetail::numpyInputCatalog::snapshot&);

        //- Determine a registered or template field class
        word fieldClass(const word& fieldName) const;

        //- Find or construct a volume field from its template
        template<class GeoField>
        GeoField& field
        (
            const word& fieldName,
            HashPtrTable<GeoField>& ownedFields
        );

        //- Import one typed field
        template<class GeoField>
        void loadField
        (
            const numpyDetail::numpyInputCatalog::snapshot& sample,
            const word& fieldName,
            HashPtrTable<GeoField>& ownedFields,
            bool validateOnly = false
        );

        //- Correct patches unless they require an active model calculation
        template<class GeoField>
        void correctFieldBoundaryConditions(GeoField& field) const;

        //- Write a registered imported field at the current time
        template<class GeoField>
        bool writeField(const word& fieldName);

        //- No copy construct
        numpyToFoam(const numpyToFoam&) = delete;

        //- No copy assignment
        void operator=(const numpyToFoam&) = delete;


public:

    //- Runtime type information
    TypeName("numpyToFoam");


    // Constructors

        //- Construct from Time and dictionary
        numpyToFoam
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct for an object registry
        numpyToFoam
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~numpyToFoam();


    // Member Functions

        //- Read the function-object dictionary
        virtual bool read(const dictionary& dict);

        //- Load the snapshot at the current time
        virtual bool execute();

        //- Optionally write imported fields in native format
        virtual bool write();

        //- Explicitly load a snapshot (used by numpyPostProcess)
        bool load(scalar timeValue);

        //- Correct patch fields after the driver environment is available
        void correctBoundaryConditions();

        //- Available committed input times
        scalarList times() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "numpyToFoam.txx"

#endif

// ************************************************************************* //
