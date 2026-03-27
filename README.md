# numpyToFoam

## General Description
`numpyToFoam` is a custom OpenFOAM utility for reconstructing OpenFOAM field data from NumPy (`.npy`) arrays. It supports scalar, vector, symmetric tensor, and full tensor fields.

The field values over the computational cells at a given time instant constitute a **snapshot**. A collection of such snapshots over multiple time steps forms the input dataset used by the utility.

These snapshots are processed sequentially in time, avoiding the need to store large datas

---

## Dependencies
- OpenFOAM v2406
- A working C++ compiler available in the OpenFOAM environment (e.g. `g++`)

---

## Compilation and Installation
In a working OpenFOAM environment:

```bash
cd src
wmake
```

The utility will be compiled and made available within your OpenFOAM environment.

---

## Usage

### Requirements

To use the utility, the following are required:

### 1. OpenFOAM Case
- A decomposed OpenFOAM case (parallel setup required)
- Mesh must be available in `processor*/` directories
- Optional: `0/` directory with initial and boundary conditions

### 2. NumPy Snapshot Files
- Data must be stored as `.npy` files
- Naming convention:

```bash
fieldName_proc_i.npy
```

where:
- `fieldName` is the OpenFOAM field name, for example `p` or `U`
- `i` is the processor index

If there are:
- `p` processors
- `n` fields

then the total number of `.npy` files is `n × p`.

### Supported Shapes

| Field Type        | Shape                     |
|-------------------|---------------------------|
| Scalar            | `(nCells, nTimeSteps)`    |
| Vector            | `(nCells, 3, nTimeSteps)` |
| Symmetric Tensor  | `(nCells, 6, nTimeSteps)` |
| Tensor            | `(nCells, 9, nTimeSteps)` |

### Storage Format
- Both **row-major (C-order)** and **column-major (Fortran-order)** storage formats are supported
- The utility automatically detects the array storage format from the `.npy` header
- **Recommended:** use column-major format for faster snapshot access

### 3. `numpyToFoamDict`
A `numpyToFoamDict` file must be present in the `system/` directory.

This dictionary controls:
- the location of the `.npy` files (`dataDir`)
- the list of fields to reconstruct
- the time range and time-step size used for writing OpenFOAM fields

Example:

```plaintext
dataDir    data;

fields     (p U);

time
{
    startTime   0;
    endTime     1;
    deltaT      0.1;
}
```

The number of generated time steps must match the number of snapshots stored in the `.npy` files.

---

## Running the Utility

```bash
mpirun -np 2 numpyToFoam -parallel
```

Replace `2` with the number of processors used for case decomposition.

---

## Limitations
1. **Parallel-only execution**
   - Serial execution is not supported

2. **Precision conversion**
   - If the input `.npy` files are stored in single precision (`float32`), the values are converted to OpenFOAM `scalar` precision during reading
   - In most OpenFOAM builds, this means conversion to **double precision**

3. **No in-built post-processing**
   - Derived quantities cannot be computed directly within the utility
   - Post-processing must be performed separately after reconstruction