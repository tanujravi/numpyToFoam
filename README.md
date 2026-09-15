# OpenFOAM NumPy function objects

This project provides native OpenFOAM function objects for exchanging field
snapshots with NumPy while a solver is running or during post-processing:

- `foamToNumpy`: volume-field export
- `numpyToFoam`: volume-field import
- `finiteAreaToNumpy`: finite-area face/edge-field export
- `finiteAreaToFoam`: finite-area face/edge-field import
- `numpyPostProcess`: a driver that imports NumPy snapshots and immediately
  runs ordinary OpenFOAM function objects without staging native field files

The standalone `foamToNumpy` and `numpyToFoam` applications are retained for
legacy workflows. New development targets the function objects and
`numpyPostProcess`.

All new NumPy arrays use Fortran order. C-order output is not available and
C-order input is rejected.

## Compatibility and build

The supported release matrix is OpenFOAM v2506, v2512, and v2606. The
implementation also builds against the current development line. On
development (API 2607 and newer), it uses OpenFOAM's
`Foam::fileFormats::numpyCore`. On v2606 and earlier releases, a private
compatibility implementation supplies the same small header-writing interface;
a development installation is therefore not needed to compile a release
version.

From an initialized OpenFOAM environment:

```bash
./Allwmake
```

This builds `libnumpyFunctionObjects`, `numpyPostProcess`, and the two legacy
applications. Use `./Allwclean` to clean all targets.

## Volume export

Add the following to `system/controlDict`, or to a separate dictionary passed
to `postProcess -dict`:

```text
functions
{
    numpyExport
    {
        type                foamToNumpy;
        libs                (numpyFunctionObjects);
        fields              (p U "T.*");
        dataType            float64;
        batchSize           100;
        writeTimes          true;
        writeCellCentres    true;
        writeCellVolumes    true;
        outputDir           "postProcessing/numpyExport";
        writeControl        adjustableRunTime;
        writeInterval       0.1;
    }
}
```

`fields` accepts OpenFOAM word/regular-expression selections. Supported field
classes are volume scalar, vector, spherical-tensor, symmetric-tensor, and
tensor fields. Only internal values are exported.

The function object's output control is independent of the solver's general
write control. Thus an export interval of 0.1 s works with a solver write
interval of 1 s without creating native time directories every 0.1 s.

Offline usage is the normal OpenFOAM workflow:

```bash
postProcess -dict system/controlDict.numpy -time '0.1:1'
mpirun -np 4 postProcess -parallel \
    -dict system/controlDict.numpy -time '0.1:1'
```

## Batched, crash-consistent layout

For volume fields the output layout is:

```text
postProcessing/numpyExport/
└── 0/                         # invocation/restart segment
    ├── segmentInfo
    ├── batch_000000/
    │   ├── state
    │   ├── times.npy
    │   ├── p_proc_0.npy
    │   └── U_proc_0.npy
    ├── batch_000001/
    └── geometry_000000/
        ├── cellCentres_proc_0.npy
        └── cellVolumes_proc_0.npy
```

Shapes are `(nCells, nOutputs)` for scalar fields and
`(nCells, nComponents, nOutputs)` for non-scalar fields. A batch contains at
most `batchSize` snapshots. The atomic `state` dictionary contains the
authoritative committed count; consumers must ignore uncommitted trailing data
after an interrupted write.

Each invocation creates a new non-overwriting segment. For example, if NumPy
output is every 0.1 s, native solver output is every 1 s, and a run crashes at
1.2 s before restarting from 1 s, the segments contain:

```text
0/  -> 0.1 ... 1.0 1.1 1.2
1/  ->             1.1 1.2 1.3 ...
```

Both histories remain recoverable. Import can name several `segments` in
precedence order; later segments replace duplicate times.

## Volume import function object

`numpyToFoam` imports the batched output and registers fields for subsequent
function objects:

```text
functions
{
    numpyImport
    {
        type                numpyToFoam;
        libs                (numpyFunctionObjects);
        inputDir            "postProcessing/numpyExport";
        segment             "0";
        // segments         (0 1);  // later entries replace duplicates
        fields              (p U);
        templateInstance    "0";
        writeFields         false;
        correctBoundaryConditions false;
        executeControl      timeStep;
        writeControl        timeStep;
    }
}
```

An existing registered field is updated in place. Otherwise, its file in
`templateInstance` supplies the field class, dimensions, and boundary
conditions. Boundary values are not present in the NumPy schema. Patch-field
correction is disabled by default because model-coupled conditions such as
wall functions and `fixedFluxPressure` must be updated by their owning models.
Set `correctBoundaryConditions true` only when every imported field can be
corrected independently. `writeFields` is opt-in.

## Cell-zone subsets

Volume function objects can exchange separate arrays for named cell zones:

```text
// In a foamToNumpy function-object dictionary:
fields          (p U);
cellZones       (fluidCore "porous.*");

// In numpyToFoam, or numpyPostProcessDict.input:
inputDir        "postProcessing/numpyExport";
segment         "0";
fields          (p U);
cellZones       (fluidCore);
templateInstance "0";
```

Omitting `cellZones` retains whole-mesh exchange. Every supplied name or regular
expression must match at least one zone; an empty selection is an error.
Existing empty zones are valid, including zones with no cells on a particular
processor. Zone imports require explicit `cellZones` and may select only some
of the exported zones. Whole-mesh datasets cannot be read as zone datasets.

Each exported zone contains the same selected fields. Overlapping export zones
are allowed and duplicate their shared values. **Zones selected together for
import must be disjoint**, even when shared values agree. Imported values replace
only selected internal cells. Other cells retain registered field values, or
values from `templateInstance` when a field is first created. Subsequent solver,
boundary, and `numpyPostProcess` model corrections still operate on full fields
and can change values outside the imported zones.

Zone output retains shared batch times and a single atomic commit state:

```text
<segment>/
├── segmentInfo                 # zoneLayoutVersion 1, region, nProcs
├── batch_000000/
│   ├── state                   # count, meshRevision, cellZones, fieldClasses
│   ├── times.npy
│   └── cellZones/fluidCore/
│       ├── p_proc_0.npy
│       └── U_proc_0.npy
└── geometry_000000/
    ├── mesh_proc_0             # nCells, meshRevision
    └── cellZones/fluidCore/
        ├── cellIds_proc_0.npy  # mandatory, even with geometry output disabled
        ├── cellCentres_proc_0.npy
        └── cellVolumes_proc_0.npy
```

Cell IDs are one-dimensional, signed 64-bit integer arrays in ascending
processor-local cell order. Fields and optional geometry use exactly that row
order. Field shapes are `(nZoneCells, nOutputs)` for scalar and spherical-tensor
fields and `(nZoneCells, nComponents, nOutputs)` otherwise. Geometry retains the
existing trailing singleton snapshot dimension. All arrays explicitly specify
Fortran order; scripts rewriting even one-dimensional ID arrays must preserve
that header flag. Floating-point `dataType` applies only to fields and geometry.

Topology callbacks and changes to selected zone membership close the current
batch and produce a new mapping revision, including when membership changes
without changing its size. Import validates region, processor count, local mesh
size, sorted cell IDs, zone membership, field classes, and payloads before
assigning zone values. Restart precedence and committed-count rules also apply
to zone data. Selection and output settings cannot change after export starts.

Import requires the original mesh cell numbering and parallel decomposition;
it performs no interpolation, redistribution, zone creation, or mesh advancement.
Membership checks cannot establish equivalence of arbitrary different meshes.
When importing changing mappings, the caller must provide the matching current
mesh. `numpyPostProcess` does not reconstruct or advance mesh topology.

Storage scales with the sum of selected zone sizes. Export and import process
one zone/field buffer at a time; import preflights payloads before rereading them
for assignment to avoid partially updating fields on malformed input. This
feature applies to the volume function objects and driver, not the legacy
standalone applications or finite-area fields.

## `numpyPostProcess`

The driver imports fields and runs the `functions` dictionary in the same
object registry. It avoids the NumPy -> native OpenFOAM -> function-object
read/write cycle.

For solver environments, the driver constructs thermo, transport, and
turbulence objects from `templateInstance` before importing a snapshot. It
then updates the model-owned fields, corrects independently evaluable patch
fields, and revalidates the turbulence model. Effective viscosity and its wall
functions are therefore refreshed for every imported time. Conditions owned
by a discretised equation (`epsilonWallFunction`, `omegaWallFunction`, and
`fixedFluxPressure`) retain their template values because their update requires
equation data which does not exist during function-object post-processing. The
`none` environment has no model context and is intended for independently
correctable fields.

```text
input
{
    inputDir            "postProcessing/numpyExport";
    segment             "0";
    fields              (p U);
    templateInstance    "0";
}

environment
{
    type                pimpleFoam;
}

functions
{
    importedForces
    {
        type            forces;
        libs            (forces);
        patches         (body);
        rho             rhoInf;
        rhoInf          1;
        CofR            (0 0 0);
        writeFields     true;
    }
}
```

Run it serially or on the original decomposition:

```bash
numpyPostProcess -time '0.1:1'
mpirun -np 4 numpyPostProcess -parallel -time '0.1:1'
```

The environment adapters create the solver-owned objects required by
function objects:

| Environment | Solver aliases | Objects provided |
| --- | --- | --- |
| `none` | — | Imported fields only |
| `incompressible` | `simpleFoam`, `pimpleFoam` | `phi`, single-phase transport, optional turbulence |
| `compressible` | `rhoSimpleFoam`, `rhoPimpleFoam` | fluid thermo, `rho`, `phi`, optional compressible turbulence |
| `buoyantCompressible` | `buoyantPimpleFoam` | rho thermo, `rho`, `phi`, gravity, optional compressible turbulence |

Laminar and turbulent cases share an adapter. The OpenFOAM runtime model
selector reads `turbulenceProperties`; no duplicate laminar/turbulent adapter
hierarchy is needed.

Setting `writeFields true` in `input` makes the driver write the imported
volume fields after correction. The default remains an entirely in-memory
workflow.

## Finite-area fields

`finiteAreaToNumpy` and `finiteAreaToFoam` support scalar, vector,
spherical-tensor, symmetric-tensor, and tensor area and edge fields. The
layout adds the finite-area mesh name:

```text
postProcessing/areaExport/<area>/<segment>/batch_000000/
```

Example export:

```text
areaExport
{
    type                finiteAreaToNumpy;
    libs                (numpyFunctionObjects);
    area                region0;
    areaFields          (Ts hs);
    edgeFields          (phis phi2s);
    outputDir           "postProcessing/areaExport";
    batchSize           20;
    writeAreaCentres    true;
    writeFaceAreas      true;
    writeEdgeCentres    true;
}
```

`writeFaceAreas` writes the positive finite-area face measures from
`faMesh::S()` as `faceAreas_proc_<rank>.npy`. Together with
`writeAreaCentres`, these values support area-weighted reductions without
reading the native finite-area mesh.

Example import in `numpyPostProcessDict`:

```text
finiteAreaInput
{
    inputDir            "postProcessing/areaExport";
    segment             "0";
    area                region0;
    areaFields          (Ts hs);
    edgeFields          (phis phi2s);
    templateInstance    "0";
}
```

`input` and `finiteAreaInput` are independently optional, but at least one is
required. When both are present their committed time sets must match. Multiple
finite-area meshes can be selected with `areas (...)`.

The buoyant adapter links the standard finite-area region models. This allows
the `hotRoomWithThermalShell` setup to instantiate its thermal-shell model,
import `Ts` and `hs`, and run standard consumers such as `areaWrite` directly.
The v2506 tutorial uses the older model-qualified names `Ts_ceilingShell` and
`h_ceilingShell`; the integration suite detects this layout and tests those
fields instead.

## Verification

The extended integration suite uses OpenFOAM tutorial cases and performs
export -> in-memory import -> re-export -> SHA-256/array equality checks:

- the release-native parallel `icoFoam` cavity tutorial for volume batching,
  online cadence, volume round-trip, and incompressible forces;
- `rhoSimpleFoam/squareBend` for compressible forces, including pressure,
  viscous, and total-force comparisons at two times after wall-viscosity
  regeneration;
- `buoyantPimpleFoam/hotRoomWithThermalShell` for turbulent and laminar
  buoyant-compressible environments, the thermal-shell path, and finite-area
  face-area geometry;
- `finiteArea/liquidFilmFoam/cylinder` for area and edge fields and their
  finite-area geometry.

The cavity mesh, initial fields, and numerical dictionaries are copied from
the active release's `$FOAM_TUTORIALS`. The repository retains only a small
test overlay containing the NumPy dictionaries, four-way decomposition, and
the transport-model selector required by the post-processing environment.

After building, run:

```bash
cd unittest
./Allrun-functionObjects
```

For isolated local testing against all supported release versions, use the
Apptainer runner from an ordinary host shell. In addition to Apptainer, the
runner requires `mksquashfs` to build images and SquashFS FUSE support to mount
cached SIF images:

```bash
sudo apt-get install squashfs-tools squashfuse
```

`fuse2fs` is not required because the test images use SquashFS, not EXT3.

On Ubuntu 23.10 and newer, an Apptainer installation built from source also
needs an AppArmor profile permitting its starter to create user namespaces.
GitHub release `.deb` packages install this profile automatically. For a source
build installed under the default `/usr/local` prefix, create
`/etc/apparmor.d/apptainer` with:

```text
# Permit unprivileged user namespace creation for Apptainer starter
abi <abi/4.0>,
include <tunables/global>

profile apptainer /usr/local/libexec/apptainer/bin/starter{,-suid} flags=(unconfined) {
    userns,

    include if exists <local/apptainer>
}
```

Reload AppArmor after installing the profile:

```bash
sudo systemctl reload apparmor
```

See the official
[Apptainer installation instructions](https://github.com/apptainer/apptainer/blob/main/INSTALL.md#apparmor-profile-ubuntu-2310)
for other installation prefixes and the alternative system-wide setting.

After satisfying these host requirements, run:

```bash
cd unittest
./unittest.sh
```

It builds and caches NumPy-enabled OpenFOAM v2506, v2512, and v2606 images in
`unittest/of_versions`, then builds the current source and runs the complete
integration suite independently in each image. Select a subset with, for
example, `OPENFOAM_VERSIONS="2512 2606" ./unittest.sh`. Set
`REBUILD_IMAGES=1` to rebuild the cached images. If rootless image creation is
not configured, use `APPTAINER_BUILD_WITH_SUDO=1 ./unittest.sh`; executing
cached images does not require root. Per-version build and test logs are written
below `unittest/run`.

The suite has been run successfully with v2506, v2512, v2606, and the current
development line. GitHub Actions builds the function-object library and driver
and runs this complete suite with v2506, v2512, and v2606.

## Tutorials

The examples under `tutorials/pod_time` and `tutorials/pod_parametric` use
online `foamToNumpy` function objects and `numpyPostProcess`. Their Python
scripts contain only the reduced-order calculation and batch assembly; solver
model setup and post-processing remain in OpenFOAM dictionaries.

## Legacy applications

The original standalone applications and their dictionaries remain in
`src/foamToNumpy`, `src/numpyToFoam`, and the existing unit tests. They use the
old, non-batched schema and should be considered compatibility code. Once
downstream workflows have migrated, they can be removed independently of the
function-object library. Their remaining read/write paths also reject C-order
arrays.
