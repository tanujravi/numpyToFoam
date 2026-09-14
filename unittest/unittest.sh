#!/usr/bin/env bash
#------------------------------------------------------------------------------

set -uo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd -- "$script_dir/.." && pwd)"
image_dir="$script_dir/of_versions"
run_dir="$script_dir/run"
definition="$script_dir/openfoam-test.def"
build_jobs="${JOBS:-$(nproc)}"
read -r -a versions <<< "${OPENFOAM_VERSIONS:-2506 2512 2606}"

if ! command -v apptainer >/dev/null 2>&1; then
    echo "Apptainer was not found on PATH" >&2
    exit 1
fi

if ! command -v mksquashfs >/dev/null 2>&1; then
    echo "mksquashfs was not found; install the squashfs-tools package" >&2
    exit 1
fi

mkdir -p "$image_dir" "$run_dir"

overall_status=0

apptainer_build=(apptainer build)
if [[ "${APPTAINER_BUILD_WITH_SUDO:-0}" == 1 ]]; then
    apptainer_build=(sudo apptainer build)
fi

for version in "${versions[@]}"; do
    image="$image_dir/openfoam${version}-test.sif"
    image_log="$image_dir/openfoam${version}-test.image.log"
    version_dir="$run_dir/of${version}"
    work_dir="$version_dir/work"
    build_log="$version_dir/build.log"
    test_log="$version_dir/test.log"

    mkdir -p "$version_dir"

    if [[ ! -f "$image" || "${REBUILD_IMAGES:-0}" == 1 ]]; then
        echo "Building OpenFOAM $version test image"
        if ! "${apptainer_build[@]}" --force \
            --build-arg "openfoam_version=$version" \
            "$image" "$definition" >"$image_log" 2>&1; then
            echo "OpenFOAM $version: image build failed; see $image_log" >&2
            tail -n 80 "$image_log" >&2
            overall_status=1
            continue
        fi
    else
        echo "Using cached image $image"
    fi

    rm -rf -- "$work_dir"
    mkdir -p "$work_dir/unittest"
    cp -a "$repo_dir/src" "$repo_dir/Allwmake" "$work_dir/"
    find "$work_dir/src" -type d -path '*/Make/linux*' -prune \
        -exec rm -rf -- {} +
    cp -a \
        "$script_dir/Allrun-functionObjects" \
        "$script_dir/check_cell_zones.py" \
        "$script_dir/zoneLifecycle.C" \
        "$script_dir/check_forces.py" \
        "$script_dir/check_function_object.py" \
        "$script_dir/check_numpy_roundtrip.py" \
        "$script_dir/integration" \
        "$script_dir/prepare_cavity" \
        "$work_dir/unittest/"

    echo "Building with OpenFOAM $version"
    if ! apptainer exec --cleanenv "$image" bash -lc "
        source /usr/lib/openfoam/openfoam${version}/etc/bashrc
        set -eo pipefail
        export FOAM_USER_APPBIN=\"$version_dir/platforms/\$WM_OPTIONS/bin\"
        export FOAM_USER_LIBBIN=\"$version_dir/platforms/\$WM_OPTIONS/lib\"
        export PATH=\"\$FOAM_USER_APPBIN:\$PATH\"
        export LD_LIBRARY_PATH=\"\$FOAM_USER_LIBBIN:\$LD_LIBRARY_PATH\"
        export WM_NCOMPPROCS='$build_jobs'
        mkdir -p \"\$FOAM_USER_APPBIN\" \"\$FOAM_USER_LIBBIN\"
        cd '$work_dir'
        ./Allwmake
    " >"$build_log" 2>&1; then
        echo "OpenFOAM $version: build failed; see $build_log" >&2
        tail -n 80 "$build_log" >&2
        overall_status=1
        continue
    fi

    echo "Running the integration suite with OpenFOAM $version"
    if ! apptainer exec --cleanenv "$image" bash -lc "
        source /usr/lib/openfoam/openfoam${version}/etc/bashrc
        set -eo pipefail
        export FOAM_USER_APPBIN=\"$version_dir/platforms/\$WM_OPTIONS/bin\"
        export FOAM_USER_LIBBIN=\"$version_dir/platforms/\$WM_OPTIONS/lib\"
        export PATH=\"\$FOAM_USER_APPBIN:\$PATH\"
        export LD_LIBRARY_PATH=\"\$FOAM_USER_LIBBIN:\$LD_LIBRARY_PATH\"
        export OMPI_MCA_rmaps_base_oversubscribe=1
        cd '$work_dir'
        ./unittest/Allrun-functionObjects
    " >"$test_log" 2>&1; then
        echo "OpenFOAM $version: tests failed; see $test_log" >&2
        tail -n 120 "$test_log" >&2
        overall_status=1
        continue
    fi

    echo "OpenFOAM $version: build and tests passed"
done

if (( overall_status != 0 )); then
    echo "One or more OpenFOAM versions failed" >&2
else
    echo "All requested OpenFOAM versions passed"
fi

exit "$overall_status"

#------------------------------------------------------------------------------
