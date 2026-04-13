#!/bin/bash
#------------------------------------------------------------------------------

# definition of path variables
main_folder="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
echo $main_folder
run_folder="$main_folder/run"
simulation_base_name="$main_folder/of_cavity"
versions_folder_name="$main_folder/of_versions"

SRC_FOLDER="$main_folder/../src"

# create run directory if necessary
echo "Starting test"

mkdir $run_folder
versions=(2112 2206 2312 2412 2512)
for version in "${versions[@]}"; do

    foamToNumpy_build_ok=0
    numpyToFoam_build_ok=0
    allrun_ok=0
    foamToNumpy_run_ok=0
    clean_proc_data_ok=0
    numpyToFoam_run_ok=0

    cd $run_folder
    version_folder="$run_folder/of$version"
    mkdir $version_folder
    cd $version_folder
    image="$versions_folder_name/openfoam$version.sif"
    source="source /usr/lib/openfoam/openfoam$version/etc/bashrc"
    
    log_foamToNumpy_build="$version_folder/log.foamToNumpy.build"
    log_numpyToFoam_build="$version_folder/log.numpyToFoam.build"
    log_allrun="$version_folder/log.Allrun"
    log_foamToNumpy_run="$version_folder/log.foamToNumpy"
    log_numpyToFoam_run="$version_folder/log.numpyToFoam"

    : > "$log_foamToNumpy_build"
    : > "$log_numpyToFoam_build"
    : > "$log_allrun"
    : > "$log_foamToNumpy_run"
    : > "$log_numpyToFoam_run"

    cp -r $SRC_FOLDER .
    FOAMTONUMPY="$version_folder/src/foamToNumpy"
    cd $FOAMTONUMPY
    apptainer exec $image bash -lc "$source && wmake" \
            > "$log_foamToNumpy_build" 2>&1
    foamToNumpy_build_ok=$(( $? == 0 ))

    NUMPYTOFOAM="$version_folder/src/numpyToFoam"
    cd $NUMPYTOFOAM
    apptainer exec $image bash -lc "$source && wmake" \
            > "$log_numpyToFoam_build" 2>&1
    numpyToFoam_build_ok=$(( $? == 0 ))

    cd $version_folder
    cp -r $simulation_base_name .
    
    cd $version_folder/of_cavity

    apptainer exec $image bash -lc "$source && ./Allrun" > "$log_allrun" 2>&1
    allrun_ok=$(( $? == 0 ))

    mpirun -np 4 apptainer exec $image bash -lc "$source && foamToNumpy -parallel" > "$log_foamToNumpy_run" 2>&1
    foamToNumpy_run_ok=$(( $? == 0 ))

    apptainer exec $image bash -lc "$source && ./Clean_proc_data"

    mpirun -np 4 apptainer exec $image bash -lc "$source && numpyToFoam -parallel" > "$log_numpyToFoam_run" 2>&1
    numpyToFoam_run_ok=$(( $? == 0 ))

    echo "----------------------------------------"
    echo "Version: $version"
    echo "foamToNumpy build      : $foamToNumpy_build_ok"
    echo "numpyToFoam build      : $numpyToFoam_build_ok"
    echo "Allrun                 : $allrun_ok"
    echo "foamToNumpy run        : $foamToNumpy_run_ok"
    echo "numpyToFoam run        : $numpyToFoam_run_ok"
    echo "----------------------------------------"

done


