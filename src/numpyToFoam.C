/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 AUTHOR,AFFILIATION
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

Application
    numpyToFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OSspecific.H"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "readData.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "functionObjectList.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);
    
    const IOdictionary dict
    (
        IOobject::selectIO
        (
            IOobject
            (
                "numpyToFoamDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            ),
            args.getOrDefault<fileName>("dict", "")
        )
    );
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),   // constant/
            mesh,                 // register into mesh database
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const bool is_write(dict.get<bool>("write"));
    wordList fields(dict.get<wordList>("fields"));
    const bool is_postProcess(dict.get<bool>("postProcess"));
    word order(dict.get<word>("order"));
    const dictionary& timeDict = dict.subDict("time");
    const scalar t_start = timeDict.get<scalar>("startTime");
    const scalar t_end   = timeDict.get<scalar>("endTime");
    const scalar dt      = timeDict.get<scalar>("deltaT");
    label procNo = Pstream::myProcNo();
    fileName caseDir = runTime.path();     // .../case/processorN
    fileName rootCaseDir = caseDir.path(); // .../case
    fileName dataDir = rootCaseDir/"data";

    fileNameList files = readDir(dataDir, fileName::FILE);
    
    label nSteps = floor((t_end - t_start)/dt + 0.5) + 1;
    scalarList timeValues(nSteps);
    for (label i = 0; i < nSteps; ++i)
    {
        timeValues[i] = t_start + i*dt;
    }
    timeValues[nSteps-1] = t_end;

    Info << "Generated times: " << timeValues << endl;

    
    HashTable<NpyMeta> metaByField; // fieldName -> meta
    forAll(fields, fieldInd)
    {
        word fieldName = fields[fieldInd];
        word fname = fieldName + "_proc_" + Foam::name(procNo) + ".npy";

        fileName fullPath = dataDir/fname;

        if (!isFile(fullPath))
        {
            FatalErrorInFunction
                << "Missing file: " << fullPath << exit(FatalError);
        }

        NpyMeta meta = readNpyMeta(fullPath);

        metaByField.insert(fieldName, meta);

        // ---- PRINT INFO ----
        Info<< "Rank " << Pstream::myProcNo()
            << " | Loaded file: " << meta.file
            << " | Shape: (";

        forAll(meta.shape, d)
        {
            Info<< meta.shape[d];
            if (d != meta.shape.size()-1)
                Info<< " x ";
        }

        Info<< ")"
            << " | Type: " << (meta.is_f8 ? "float64" : "float32")
            << nl;
    }

    runTime.setTime(timeValues[0], 0);

    HashTable<autoPtr<volScalarField>> scalarTemplates, scalarFields;
    HashTable<autoPtr<volVectorField>> vectorTemplates, vectorFields;
    HashTable<autoPtr<volSymmTensorField>> symmTensorTemplates, symmTensorFields;

    forAll(fields, fieldInd)
    {
        const word fieldName = fields[fieldInd];
        const NpyMeta& meta  = metaByField[fieldName];

        if (meta.shape.size() == 2)
        {
            scalarTemplates.insert
            (
                fieldName,
                autoPtr<volScalarField>
                (
                    new volScalarField
                    (
                        IOobject(fieldName, "0", mesh, IOobject::MUST_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
                        mesh
                    )
                )
            );

            scalarFields.insert
            (
                fieldName,
                autoPtr<volScalarField>
                (
                    new volScalarField
                    (
                        IOobject(fieldName, runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, IOobject::REGISTER),
                        scalarTemplates[fieldName]()
                    )
                )
            );
        }
        else if (meta.shape.size() == 3)
        {
            vectorTemplates.insert
            (
                fieldName,
                autoPtr<volVectorField>
                (
                    new volVectorField
                    (
                        IOobject(fieldName, "0", mesh, IOobject::MUST_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
                        mesh
                    )
                )
            );

            vectorFields.insert
            (
                fieldName,
                autoPtr<volVectorField>
                (
                    new volVectorField
                    (
                        IOobject(fieldName, runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, IOobject::REGISTER),
                        vectorTemplates[fieldName]()
                    )
                )
            );
        }
        else if (meta.shape.size() == 4)
        {
            symmTensorTemplates.insert
            (
                fieldName,
                autoPtr<volSymmTensorField>
                (
                    new volSymmTensorField
                    (
                        IOobject(fieldName, "0", mesh, IOobject::MUST_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
                        mesh
                    )
                )
            );

            symmTensorFields.insert
            (
                fieldName,
                autoPtr<volSymmTensorField>
                (
                    new volSymmTensorField
                    (
                        IOobject(fieldName, runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, IOobject::REGISTER),
                        symmTensorTemplates[fieldName]()
                    )
                )
            );
        }
    }

    autoPtr<functionObjectList> functionsPtr;

    if (is_postProcess)
    {
        functionsPtr.reset(new functionObjectList(runTime, true));
        functionsPtr->start();
    }

    forAll(timeValues, timeIndex)

    {
        runTime.setTime(timeValues[timeIndex], timeIndex);

        forAll(fields, fieldInd)
        {   
            const word fieldName = fields[fieldInd];
            const NpyMeta& meta = metaByField[fieldName];

            // Now you can read based on meta.shape.size():
            if (meta.shape.size() == 2)
            {
                scalarField snapshot;
                readScalarSnapshot(meta, timeIndex, snapshot);
                volScalarField& scalar_ = scalarFields[fieldName]();

                forAll (scalar_.internalFieldRef(), cellI) 
                { 
                    scalar_.internalFieldRef()[cellI] = snapshot[cellI]; 
                }
                scalar_.correctBoundaryConditions();
                if (is_write) scalar_.write();                    
            }

            else if (meta.shape.size() == 3)
            {
                vectorField snapshot;
                readVectorSnapshot(meta, timeIndex, snapshot);
                volVectorField& vector_ = vectorFields[fieldName]();

                forAll (vector_.internalFieldRef(), cellI) 
                { 
                    vector_.internalFieldRef()[cellI] = snapshot[cellI]; 
                }
                vector_.correctBoundaryConditions();
                if (is_write) vector_.write();

            }

            else if (meta.shape.size() == 4)
            {
                symmTensorField snapshot;
                readSymmTensorSnapshot(meta, timeIndex, snapshot);
                volSymmTensorField& tensor_ = symmTensorFields[fieldName]();
                forAll (tensor_.internalFieldRef(), cellI) 
                { 
                    tensor_.internalFieldRef()[cellI] = snapshot[cellI]; 
                }
                tensor_.correctBoundaryConditions();
                if (is_write) tensor_.write();                    

                }
                else
                {
                    FatalErrorInFunction
                        << "Unsupported dims " << meta.shape.size()
                        << " for " << meta.file << exit(FatalError);
                }

        }

        if (is_postProcess)
        {
            functionsPtr->execute();
        }
        
    }
    if (is_postProcess)
    {
        functionsPtr->end();
    }

    return 0;
}


// ************************************************************************** //
