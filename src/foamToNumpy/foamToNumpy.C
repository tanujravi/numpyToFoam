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
    foamToNumpy

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"
#include "argList.H"
#include "writeData.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption("dict", "file", "Alternative foamToNumpyDict");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);
    const word dictName("foamToNumpyDict");

    fileName dictPath;

    if (args.readIfPresent("dict", dictPath))
    {
        // Dictionary specified on the command-line ...
 
        if (isDir(dictPath))
        {
            dictPath /= dictName;
        }
    }
    else
    {
        // Assume dictionary is to be found in the system directory
 
        dictPath = runTime.system()/dictName;
    }

    
    IOobject DictIO
    (
        dictPath,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (!DictIO.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorInFunction
            << DictIO.objectPath() << nl
            << exit(FatalError);
    }

    Info<< "Reading foamToNumpy settings from "
        << DictIO.objectRelPath() << endl;
        
        
    const IOdictionary dict(DictIO);

    wordList fields;

    word dataTypeWord_field = "float64";
    npyType dtype_field = parseNpyType(dataTypeWord_field);
    
    if (dict.found("fields"))
    {
        const dictionary& fieldDict = dict.subDict("fields");
        fields = fieldDict.get<wordList>("names");
        dataTypeWord_field = fieldDict.getOrDefault<word>("dataType", "float64");
        dtype_field = parseNpyType(dataTypeWord_field);
    }


    const dictionary& timeDict = dict.subDict("time");
    const scalar t_start = timeDict.get<scalar>("startTime");
    const scalar t_end   = timeDict.get<scalar>("endTime");
    const label every      = timeDict.get<label>("every");
    
    bool writeCellCentre = false;
    bool writeCellVolumes = false;
    bool writeWriteTimes = false;
    word dataTypeWord_export = "float64";

    npyType dtype_export = parseNpyType(dataTypeWord_export);

    if (dict.found("exportData"))
    {
        const dictionary& exportDict = dict.subDict("exportData");
        writeCellCentre = exportDict.getOrDefault<bool>("cellCentre", false);
        writeCellVolumes = exportDict.getOrDefault<bool>("cellVolumes", false);
        writeWriteTimes = exportDict.getOrDefault<bool>("writeTimes", false);
        dataTypeWord_export = exportDict.getOrDefault<word>("dataType", "float64");
        dtype_export = parseNpyType(dataTypeWord_export);
    }

    const word storageOrderWord = dict.getOrDefault<word>("storageOrder", "C");
    const bool fortranOrder = parseFortranOrder(storageOrderWord);
    fileName caseDir = runTime.path();     // .../case/processorN
    fileName rootCaseDir = caseDir.path(); // .../case
    fileName dataSubDir = dict.getOrDefault<fileName>("dataDir", "data");
    fileName dataDir;
    
    if (dataSubDir.isAbsolute())
    {
        dataDir = dataSubDir;
    }
    else
    {
        dataDir = rootCaseDir/dataSubDir;
    }

    if (Pstream::master() && !isDir(dataDir))
    {
        mkDir(dataDir);
    }
    Pstream::barrier(UPstream::worldComm);

    instantList allTimes = runTime.times();
    instantList selectedTimes;

    label count = 0;
    if (every <= 0)
    {
        FatalErrorInFunction
            << "'every' must be greater than 0" << nl
            << exit(FatalError);
    }

    forAll(allTimes, i)
    {
        const instant& t = allTimes[i];

        // Skip constant
        if (t.name() == runTime.constant())
        {
            continue;
        }

        const scalar timeValue = t.value();

        // Keep only times inside [startTime, endTime]
        if (timeValue >= t_start && timeValue <= t_end)
        {
            if (count % every == 0)
            {
                selectedTimes.append(t);
            }
            ++count;
        }
    }

    if (selectedTimes.empty())
    {
        FatalErrorInFunction
            << "No times selected in range [" << t_start
            << ", " << t_end << "]"
            << exit(FatalError);
    }

    runTime.setTime(selectedTimes[0], 0);
    
    //List<fieldMeta> metadata(fields.size());
    DynamicList<fieldMeta> metadata;
    
    Info<< "Inspecting fields at initial selected time "
        << runTime.timeName() << nl;
    forAll(fields, i)
    {
        const word& fieldName = fields[i];
        const fileName fieldDir = dataDir/fieldName;

        bool skipField = false;
        if (Pstream::master())
        {
            if (isDir(fieldDir))
            {
                WarningInFunction
                    << "Output directory already exists for field '" << fieldName << "'"
                    << nl << "Skipping field." << nl
                    << "Directory: " << fieldDir << nl;
                skipField = true;
            }
            else
            {
                mkDir(fieldDir);
            }
        }
        Pstream::broadcast(skipField);
        Pstream::barrier(UPstream::worldComm);
        
        if (skipField)
        {
            continue;
        }

        IOobject io
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        if (!io.typeHeaderOk<regIOobject>(false))
        {
            FatalErrorInFunction
                << "Cannot read header for field " << fieldName
                << " at time " << runTime.timeName()
                << exit(FatalError);
        }

        const word cls = io.headerClassName();
        
        fieldMeta meta;
        meta.name = fieldName;
        meta.className = cls;
        meta.kind = classifyField(cls);

        if (meta.kind == fieldKind::SCALAR)
        {
            volScalarField fld(io, mesh);
            meta.nCells = fld.size();
            meta.nComp = 1;
        }
        else if (meta.kind == fieldKind::VECTOR)
        {
            volVectorField fld(io, mesh);
            meta.nCells = fld.size();
            meta.nComp = 3;
        }
        else if (meta.kind == fieldKind::SYMM_TENSOR)
        {
            volSymmTensorField fld(io, mesh);
            meta.nCells = fld.size();
            meta.nComp = 6;
        }
        else if (meta.kind == fieldKind::TENSOR)
        {
            volTensorField fld(io, mesh);
            meta.nCells = fld.size();
            meta.nComp = 9;
        }
        else
        {
            FatalErrorInFunction
                << "Unhandled field type for " << fieldName
                << " with class " << cls
                << exit(FatalError);
        }


        meta.outFile =
            fieldDir/
            (
                fieldName
              + "_proc_"
              + Foam::name(Pstream::myProcNo())
              + ".npy"
            );
        
        metadata.append(meta);

        Info<< "Field " << meta.name
            << " : class=" << meta.className
            << ", nCells=" << meta.nCells
            << ", nComp=" << meta.nComp
            << nl;
    }

    // ------------------------------------------------------------
    // write phase
    // ------------------------------------------------------------

    // actual per-field write loop
    forAll(metadata, i)
    {
        const fieldMeta& meta = metadata[i];
        const std::vector<std::size_t> shape =
            makeShape(selectedTimes.size(), meta.nCells, meta.nComp);

        const fileName fieldPattern =
            (dataDir/meta.name)/
            (
                meta.name
            + "_proc_*.npy"
            );

        Info<< "Writing field " << meta.name
            << " to " << fieldPattern << nl;            

        if (meta.kind == fieldKind::SCALAR)
        {
            npyWriter<volScalarField> writer
            (
                meta,
                shape,
                dtype_field,
                fortranOrder,
                selectedTimes.size()
            );

            forAll(selectedTimes, timei)
            {
                runTime.setTime(selectedTimes[timei], timei);

                IOobject io
                (
                    meta.name,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                );

                volScalarField fld(io, mesh);
                writer.write(fld, timei);
            }

            writer.flush();
        }
        else if (meta.kind == fieldKind::VECTOR)
        {
            npyWriter<volVectorField> writer
            (
                meta,
                shape,
                dtype_field,
                fortranOrder,
                selectedTimes.size()
            );

            forAll(selectedTimes, timei)
            {
                runTime.setTime(selectedTimes[timei], timei);

                IOobject io
                (
                    meta.name,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                );

                volVectorField fld(io, mesh);
                writer.write(fld, timei);
            }

            writer.flush();
        }
        else if (meta.kind == fieldKind::SYMM_TENSOR)
        {
            npyWriter<volSymmTensorField> writer
            (
                meta,
                shape,
                dtype_field,
                fortranOrder,
                selectedTimes.size()
            );

            forAll(selectedTimes, timei)
            {
                runTime.setTime(selectedTimes[timei], timei);

                IOobject io
                (
                    meta.name,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                );

                volSymmTensorField fld(io, mesh);
                writer.write(fld, timei);
            }

            writer.flush();
        }
        else if (meta.kind == fieldKind::TENSOR)
        {
            npyWriter<volTensorField> writer
            (
                meta,
                shape,
                dtype_field,
                fortranOrder,
                selectedTimes.size()
            );

            forAll(selectedTimes, timei)
            {
                runTime.setTime(selectedTimes[timei], timei);

                IOobject io
                (
                    meta.name,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                );

                volTensorField fld(io, mesh);
                writer.write(fld, timei);
            }

            writer.flush();
        }
    }

    if (writeWriteTimes && Pstream::master())
    {
        const fileName timesDir = dataDir/"times";

        
        if (isDir(timesDir))
        {
            WarningInFunction
                << "Output directory already exists for times export." << nl
                << "Skipping times export." << nl
                << "Directory: " << timesDir << nl;
        }        
        
        else
        {
            
            mkDir(timesDir);                                                                                                                
            
            scalarField timesFld(selectedTimes.size());

            forAll(selectedTimes, i)
            {
                timesFld[i] = selectedTimes[i].value();
            }

            fieldMeta timesMeta;
            timesMeta.name = "times";
            timesMeta.className = "scalarField";
            timesMeta.kind = fieldKind::SCALAR;
            timesMeta.nCells = timesFld.size();
            timesMeta.nComp = 1;
            timesMeta.outFile = timesDir/"times.npy";
            Info<< "Writing " << timesMeta.name
                << " to " << timesMeta.outFile << nl;

            const std::vector<std::size_t> timesShape =
                makeShape(1, timesMeta.nCells, timesMeta.nComp);

            npyWriter<scalarField> timesWriter
            (
                timesMeta,
                timesShape,
                dtype_export,
                fortranOrder,
                1
            );

            timesWriter.write(timesFld, 0);
            timesWriter.flush();
        }
    }

    if (writeCellCentre)
    {
        const fileName centreDir = dataDir/"cellCentre";
        bool skipCellCentre = false;

        if (Pstream::master())
        {
            if (isDir(centreDir))
            {
                WarningInFunction
                    << "Output directory already exists for cellCentre export." << nl
                    << "Skipping cellCentre export." << nl
                    << "Directory: " << centreDir << nl;
                skipCellCentre = true;
            }
            else
            {
                mkDir(centreDir);
            }
        }

        Pstream::broadcast(skipCellCentre);
        Pstream::barrier(UPstream::worldComm);

        if (!skipCellCentre)
        {
            const volVectorField& C = mesh.C();

            fieldMeta centreMeta;
            centreMeta.name = "cellCentre";
            centreMeta.className = volVectorField::typeName;
            centreMeta.kind = fieldKind::VECTOR;
            centreMeta.nCells = C.size();
            centreMeta.nComp = 3;
            centreMeta.outFile =
                centreDir/
                (
                    "cellCentre_proc_"
                + Foam::name(Pstream::myProcNo())
                + ".npy"
                );
                
            const fileName fieldPattern =
                (dataDir/centreMeta.name)/
                (
                    centreMeta.name
                + "_proc_*.npy"
                );

            Info<< "Writing field " << centreMeta.name
                << " to " << fieldPattern << nl;


            const std::vector<std::size_t> shape =
                makeShape(1, centreMeta.nCells, centreMeta.nComp);

            npyWriter<volVectorField> writer
            (
                centreMeta,
                shape,
                dtype_export,
                fortranOrder,
                1
            );

            writer.write(C, 0);
            writer.flush();
        }
    }

    if (writeCellVolumes)
    {
        const fileName volDir = dataDir/"cellVolumes";
        bool skipCellVolumes = false;

        if (Pstream::master())
        {
            if (isDir(volDir))
            {
                WarningInFunction
                    << "Output directory already exists for cellVolumes export." << nl
                    << "Skipping cellVolumes export." << nl
                    << "Directory: " << volDir << nl;
                skipCellVolumes = true;
            }
            else
            {
                mkDir(volDir);
            }
        }

        Pstream::broadcast(skipCellVolumes);
        Pstream::barrier(UPstream::worldComm);

        if (!skipCellVolumes)
        {
            const scalarField& V = mesh.V();

            fieldMeta volMeta;
            volMeta.name = "cellVolumes";
            volMeta.className = "scalarField";
            volMeta.kind = fieldKind::SCALAR;
            volMeta.nCells = V.size();
            volMeta.nComp = 1;
            volMeta.outFile =
                volDir/
                (
                    "cellVolumes_proc_"
                + Foam::name(Pstream::myProcNo())
                + ".npy"
                );

            const fileName fieldPattern =
                (dataDir/volMeta.name)/
                (
                    volMeta.name
                + "_proc_*.npy"
                );

            Info<< "Writing field " << volMeta.name
                << " to " << fieldPattern << nl;


            const std::vector<std::size_t> shape =
                makeShape(1, volMeta.nCells, volMeta.nComp);

            npyWriter<scalarField> writer
            (
                volMeta,
                shape,
                dtype_export,
                fortranOrder,
                1
            );

            writer.write(V, 0);
            writer.flush();
        }
    }
}
    // ************************************************************************* //
