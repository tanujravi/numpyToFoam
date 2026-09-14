// Exercise same-size membership revisions and registered-field old-time tracking.
#include "fvCFD.H"
#include "cellZoneMesh.H"
#include "foamToNumpy.H"
#include "numpyToFoam.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    using namespace Foam::functionObjects;
    const label zonei = mesh.cellZones().findZoneID("left");
    const labelList original(mesh.cellZones()[zonei]);
    labelList changed(original);
    for (label& cell : changed) cell += original.size();

    runTime.setTime(0.1, 1);
    volScalarField q
    (
        IOobject("q", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    dictionary out;
    out.add("fields", wordList({"q"}));
    out.add("cellZones", wordList({"left"}));
    out.add("outputDir", fileName("postProcessing/lifecycle"));
    out.add("batchSize", 100);
    foamToNumpy exporter("zoneLifecycle", runTime, out);
    exporter.write();
    const scalarField first(q.primitiveField());

    runTime.setTime(0.2, 2);
    mesh.cellZones()[zonei] = changed;
    q.primitiveFieldRef() += 100;
    exporter.write();
    const scalarField second(q.primitiveField());
    runTime.setTime(0.3, 3);
    exporter.readUpdate(polyMesh::TOPO_CHANGE);
    exporter.write();
    exporter.end();

    dictionary in;
    in.add("fields", wordList({"q"}));
    in.add("cellZones", wordList({"left"}));
    in.add("inputDir", fileName("postProcessing/lifecycle"));
    in.add("segment", word("0"));
    numpyToFoam importer("zoneLifecycleImport", runTime, in);

    runTime.setTime(0.1, 3);
    mesh.cellZones()[zonei] = original;
    q.primitiveFieldRef() = -99;
    importer.load(0.1);
    const scalarField importedFirst(q.primitiveField());
    forAll(q, cell)
    {
        const scalar expected = original.found(cell) ? first[cell] : -99;
        if (q[cell] != expected)
            FatalErrorInFunction << "First partial import mismatch" << exit(FatalError);
    }

    runTime.setTime(0.2, 4);
    mesh.cellZones()[zonei] = changed;
    importer.load(0.2);
    forAll(q, cell)
    {
        const scalar expected = changed.found(cell) ? second[cell] : importedFirst[cell];
        if (q[cell] != expected || q.oldTime()[cell] != importedFirst[cell])
            FatalErrorInFunction << "Outside-cell or old-time mismatch" << exit(FatalError);
    }
    Info << "Zone lifecycle checks passed" << endl;
    return 0;
}
