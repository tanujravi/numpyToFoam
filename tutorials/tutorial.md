# Tutorial

## Mode Visualization and Reconstruction Using `numpyToFoam`

One useful application of `numpyToFoam` is the visualization of POD modes computed from a snapshot matrix, as well as the visualization of low-order reconstructed fields. In this tutorial, the utility is used together with a cavity-flow simulation to reconstruct pressure fields and write POD modes back into OpenFOAM format.

The tutorial is based on the flow-across-cavity case provided in the folder `of_cavity`.

---

## 1. Run the OpenFOAM Simulation

In a working OpenFOAM environment, run:

```bash
cd of_cavity
./Allrun
```

This runs the simulation and stores the decomposed field data in the processor directories.

---

## 2. Read Snapshot Data from the OpenFOAM Case

To load the OpenFOAM field data into Python, the [`flowtorch`](https://github.com/FlowModelingControl/flowtorch) library is used.

```python
from flowtorch.data import FOAMDataloader
import numpy as np
import os

folder_path = "of_cavity"
loader = FOAMDataloader(folder_path, distributed=True)
times = loader.write_times

p_snap = np.asarray(loader.load_snapshot("p", times[1:]), dtype=np.float64)
```

Here, `p_snap` is the pressure snapshot matrix, where each column corresponds to one time instant.

---

## 3. Compute POD Modes and Low-Order Reconstruction

The POD modes are computed using the singular value decomposition (SVD) of the snapshot matrix. A low-order reconstruction is then obtained using the dominant modes.

```python
U, D, VT = np.linalg.svd(p_snap, full_matrices=False)

rank = 15
Ur = U[:, :rank]
Dr = D[:rank]
VTr = VT[:rank, :]

p_reconstructed = Ur @ np.diag(Dr) @ VTr
```

In this decomposition:

- `U` contains the spatial POD modes
- `D` contains the singular values
- `VT` contains the temporal coefficients

The matrix `p_reconstructed` is the rank-`15` approximation of the original pressure snapshot matrix.

---

## 4. Save Reconstruction and Modes as NumPy Files

The reconstructed fields and POD modes must be split according to the number of processors used in the OpenFOAM simulation and then saved as `.npy` files.

```python
n_proc = 4

os.makedirs("reconstruction_data", exist_ok=True)

p_parts = np.array_split(p_reconstructed, n_proc, axis=0)
for i, part in enumerate(p_parts):
    np.save(f"reconstruction_data/p_proc_{i}.npy", part)

n_modes = 5
os.makedirs("mode_data", exist_ok=True)

modes = U[:, :n_modes]
mode_parts = np.array_split(modes, n_proc, axis=0)

for i, part in enumerate(mode_parts):
    np.save(f"mode_data/p_modes_proc_{i}.npy", part)
```

Two datasets are now created:

- `reconstruction_data/` containing reconstructed pressure snapshots
- `mode_data/` containing the first few POD modes

For reconstructed fields, the `.npy` arrays contain both spatial and temporal information. For POD modes, each column corresponds to one mode, and the modes are later written as separate OpenFOAM time steps.

---

## 5. Write the Reconstructed Fields Back to OpenFOAM

Create a copy of the original case and clean the processor data so that only the mesh and `0/` fields remain.

```bash
cp -r of_cavity of_cavity_reconstructed
cd of_cavity_reconstructed
./Clean_proc_data
mv ../reconstruction_data .
```

The script `Clean_proc_data` removes processor time directories while retaining the decomposed mesh and the `0/` directory.

Use the following `system/numpyToFoamDict`:

```plaintext
dataDir       "reconstruction_data";
fields        (p);
write         true;

time
{
    startTime     0.01;
    endTime       0.5;
    deltaT        0.01;
}
```

Now run:

```bash
mpirun -np 4 numpyToFoam -parallel
```

This writes the reconstructed pressure snapshots into the corresponding OpenFOAM time directories.

---

## 6. Write the POD Modes to OpenFOAM

The POD modes can be written in the same way by creating a separate case:

```bash
cp -r of_cavity of_cavity_modes
cd of_cavity_modes
./Clean_proc_data
mv ../mode_data .
```

Use the following `system/numpyToFoamDict`:

```plaintext
dataDir       "mode_data";
fields        (p_modes);
write         true;

time
{
    startTime     1;
    endTime       5;
    deltaT        1;
}
```

Since POD modes do not have a physical time dimension, the time settings are used here only as labels. Each written time directory corresponds to one POD mode.

For example:

- time `1` → mode 1
- time `2` → mode 2
- time `3` → mode 3
- and so on

Run:

```bash
mpirun -np 4 numpyToFoam -parallel
```

The POD modes are now available as OpenFOAM fields and can be visualized directly in ParaView.

---

## Summary

This tutorial shows how `numpyToFoam` can be used to:

- write low-order reconstructed fields from reduced-order models back into OpenFOAM format
- write POD modes as OpenFOAM fields for visualization
