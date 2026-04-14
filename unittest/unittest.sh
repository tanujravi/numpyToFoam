#!/bin/bash
#------------------------------------------------------------------------------

# definition of path variables
main_folder="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
echo "$main_folder"

run_folder="$main_folder/run"
simulation_base_name="$main_folder/of_cavity"
versions_folder_name="$main_folder/of_versions"
SRC_FOLDER="$main_folder/../src"

echo "Starting test"
mkdir -p "$run_folder"

# ----------------------------------------------------------------------
# function: make checksum manifest for processor files
# ----------------------------------------------------------------------
makeChecksumManifest() {
    case_dir="$1"
    out_file="$2"

    (
        cd "$case_dir" || exit 1

        find processor* -type f \
            ! -path "*/constant/polyMesh/*" \
            ! -path "*/uniform/*" \
            ! -name "phi" \
            -print0 \
        | sort -z \
        | xargs -0 -r md5sum
    ) > "$out_file"
}

# ----------------------------------------------------------------------
# function: compare two checksum manifests
# ----------------------------------------------------------------------
compareChecksumManifest() {
    ref_file="$1"
    new_file="$2"
    diff_file="$3"

    ref_hashes="${ref_file}.hashes"
    new_hashes="${new_file}.hashes"

    awk '{print $1}' "$ref_file" > "$ref_hashes"
    awk '{print $1}' "$new_file" > "$new_hashes"

    diff -u "$ref_hashes" "$new_hashes" > "$diff_file" 2>&1
    status=$?

    rm -f "$ref_hashes" "$new_hashes"
    return $status
}

versions=(2112 2206 2312 2412 2512)

mkdir -p "$versions_folder_name"

for version in "${versions[@]}"; do
    sif_file="$versions_folder_name/openfoam$version.sif"

    if [ ! -f "$sif_file" ]; then
        echo "Building $sif_file"
        cd "$versions_folder_name" || exit 1
        sudo apptainer build "openfoam$version.sif" "docker://opencfd/openfoam-default:$version"
    else
        echo "$sif_file already exists, skipping build"
    fi
done

for version in "${versions[@]}"; do

    foamToNumpy_build_ok=0
    numpyToFoam_build_ok=0
    allrun_ok=0
    foamToNumpy_run_ok=0
    clean_proc_data_ok=0
    numpyToFoam_run_ok=0
    checksum_ok=0

    cd "$run_folder" || exit 1

    version_folder="$run_folder/of$version"
    mkdir -p "$version_folder"
    cd "$version_folder" || exit 1

    image="$versions_folder_name/openfoam$version.sif"
    source="source /usr/lib/openfoam/openfoam$version/etc/bashrc"

    log_foamToNumpy_build="$version_folder/log.foamToNumpy.build"
    log_numpyToFoam_build="$version_folder/log.numpyToFoam.build"
    log_allrun="$version_folder/log.Allrun"
    log_foamToNumpy_run="$version_folder/log.foamToNumpy"
    log_clean_proc_data="$version_folder/log.Clean_proc_data"
    log_numpyToFoam_run="$version_folder/log.numpyToFoam"

    checksum_ref="$version_folder/checksum.reference.md5"
    checksum_final="$version_folder/checksum.final.md5"
    checksum_diff="$version_folder/checksum.diff"

    : > "$log_foamToNumpy_build"
    : > "$log_numpyToFoam_build"
    : > "$log_allrun"
    : > "$log_foamToNumpy_run"
    : > "$log_clean_proc_data"
    : > "$log_numpyToFoam_run"

    cp -r "$SRC_FOLDER" .

    FOAMTONUMPY="$version_folder/src/foamToNumpy"
    cd "$FOAMTONUMPY" || exit 1
    apptainer exec "$image" bash -lc "$source && wmake" \
        > "$log_foamToNumpy_build" 2>&1
    foamToNumpy_build_ok=$(( $? == 0 ))

    NUMPYTOFOAM="$version_folder/src/numpyToFoam"
    cd "$NUMPYTOFOAM" || exit 1
    apptainer exec "$image" bash -lc "$source && wmake" \
        > "$log_numpyToFoam_build" 2>&1
    numpyToFoam_build_ok=$(( $? == 0 ))

    cd "$version_folder" || exit 1
    cp -r "$simulation_base_name" .

    case_dir="$version_folder/of_cavity"
    cd "$case_dir" || exit 1

    apptainer exec "$image" bash -lc "$source && ./Allrun" \
        > "$log_allrun" 2>&1
    allrun_ok=$(( $? == 0 ))

    # save reference checksums after original simulation
    if [ "$allrun_ok" -eq 1 ]; then
        makeChecksumManifest "$case_dir" "$checksum_ref"
    fi

    mpirun -np 4 apptainer exec "$image" bash -lc "$source && foamToNumpy -parallel" \
        > "$log_foamToNumpy_run" 2>&1
    foamToNumpy_run_ok=$(( $? == 0 ))

    apptainer exec "$image" bash -lc "$source && ./Clean_proc_data" \
        > "$log_clean_proc_data" 2>&1
    clean_proc_data_ok=$(( $? == 0 ))

    mpirun -np 4 apptainer exec "$image" bash -lc "$source && numpyToFoam -parallel" \
        > "$log_numpyToFoam_run" 2>&1
    numpyToFoam_run_ok=$(( $? == 0 ))

    # compare checksums after numpyToFoam
    if [ "$numpyToFoam_run_ok" -eq 1 ] && [ -f "$checksum_ref" ]; then
        makeChecksumManifest "$case_dir" "$checksum_final"

        if compareChecksumManifest "$checksum_ref" "$checksum_final" "$checksum_diff"; then
            checksum_ok=1
        else
            checksum_ok=0
        fi
    fi

    echo "----------------------------------------"
    echo "Version: $version"
    echo "foamToNumpy build      : $foamToNumpy_build_ok"
    echo "numpyToFoam build      : $numpyToFoam_build_ok"
    echo "Allrun                 : $allrun_ok"
    echo "foamToNumpy run        : $foamToNumpy_run_ok"
    echo "Clean_proc_data        : $clean_proc_data_ok"
    echo "numpyToFoam run        : $numpyToFoam_run_ok"
    echo "checksum match         : $checksum_ok"
    echo "----------------------------------------"

done