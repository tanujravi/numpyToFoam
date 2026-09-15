#!/usr/bin/env python3
"""End-to-end zone tests against real OpenFOAM meshes and independent NumPy gathers."""
from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile
import numpy as np

HERE = Path(__file__).resolve().parent
FIELDS = {'q': ('scalar', 1), 'qv': ('vector', 3), 'qs': ('sphericalTensor', 1),
          'qy': ('symmTensor', 6), 'qt': ('tensor', 9)}
HEADER = 'FoamFile { version 2.0; format ascii; class dictionary; object test; }\n'


def run(args, case, failure=None):
    result = subprocess.run(list(map(str, args)), cwd=case, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
    if failure is None:
        assert result.returncode == 0, (args, result.stdout[-10000:])
    else:
        assert result.returncode != 0 and failure in result.stdout, (args, result.stdout[-10000:])
    return result.stdout


def save(path, values):
    # np.save uses C order for arrays contiguous in both orders, including 1-D
    # arrays. Write the explicitly Fortran-only format required by the library.
    with path.open('wb') as stream:
        np.lib.format.write_array_header_1_0(stream, {
            'descr': np.lib.format.dtype_to_descr(values.dtype),
            'fortran_order': True, 'shape': values.shape})
        stream.write(values.tobytes(order='F'))


def export_object(name, root, zones='', dtype='float64'):
    return f'''{name} {{ type foamToNumpy; libs (numpyFunctionObjects);
        fields ({' '.join(FIELDS)}); outputDir "postProcessing/{root}";
        batchSize 2; dataType {dtype}; writeCellCentres true;
        writeCellVolumes true; writeControl timeStep; {zones} }}'''


def driver(case, root, output, zones, dtype='float64', segments='segment "0.1";', extra=''):
    (case/'system/zoneDriver').write_text(HEADER + f'''
input {{ inputDir "postProcessing/{root}"; {segments}
    fields ({' '.join(FIELDS)}); templateInstance "0"; {zones} }}
environment {{ type none; }}
functions {{ {extra}
{export_object('reexport', output, '', dtype)}
}}''')


def setup(case):
    run([HERE/'prepare_cavity', case], HERE)
    run(['blockMesh'], case)
    owner = (case/'constant/polyMesh/owner').read_text()
    n = int(re.search(r'nCells:\s*(\d+)', owner)[1])
    zones = {'left': list(range(0, n//10)), 'right': list(range(n-n//10, n)),
             'overlap': list(range(0, n//20)), 'empty': []}
    text = 'FoamFile { version 2.0; format ascii; class regIOobject; object cellZones; }\n'
    text += str(len(zones)) + '\n(\n'
    for name, ids in zones.items():
        text += f'{name} {{ type cellZone; cellLabels List<label> {len(ids)} ( '
        text += ' '.join(map(str, ids)) + ' ); }\n'
    (case/'constant/polyMesh/cellZones').write_text(text + ')\n')
    for ti, time in enumerate(['0', '0.1', '0.2', '0.3']):
        (case/time).mkdir(exist_ok=True)
        for name, (kind, nc) in FIELDS.items():
            values = []
            for cell in range(n):
                components = [str(-99 if ti == 0 else 1000*ti + cell*10 + c) for c in range(nc)]
                values.append(components[0] if kind == 'scalar' else '('+' '.join(components)+')')
            text = f'''FoamFile {{ version 2.0; format ascii; class vol{kind[0].upper()+kind[1:]}Field; object {name}; }}
dimensions [0 0 0 0 0 0 0];
internalField nonuniform List<{kind}> {n} ( {' '.join(values)} );
boundaryField {{ movingWall {{ type zeroGradient; }} fixedWalls {{ type zeroGradient; }}
frontAndBack {{ type empty; }} }}
'''
            (case/time/name).write_text(text)


def check(case, parallel):
    prefix = ['mpirun', '-np', '4'] if parallel else []
    flag = ['-parallel'] if parallel else []
    if parallel:
        run(['decomposePar', '-time', '0:0.3'], case)
    def process(failure=None):
        run(prefix+['numpyPostProcess']+flag+['-dict', 'zoneDriver'], case, failure)
    for dtype in ['float64', 'float32']:
        root = 'zones'+dtype
        full = 'full'+dtype
        (case/'system/zoneExport').write_text(HEADER+'functions {\n'+
            export_object('full', full, '', dtype)+ '\n'+
            export_object('zones', root, 'cellZones (".*");', dtype)+'\n}')
        run(prefix+['postProcess']+flag+['-dict', 'system/zoneExport', '-time', '0.1:0.3'], case)
        source = case/'postProcessing'/root/'0.1'
        reference = case/'postProcessing'/full/'0.1'
        empty_rank = False
        for mapping in source.glob('geometry_*/cellZones/*/cellIds_proc_*.npy'):
            zone = mapping.parent.name
            rank = mapping.stem.split('_')[-1]
            ids = np.load(mapping)
            assert ids.dtype == np.dtype('<i8') and ids.ndim == 1
            assert np.all(ids[1:] > ids[:-1])
            if zone != 'empty' and not ids.size:
                empty_rank = True
            for batch in source.glob('batch_*'):
                for name in FIELDS:
                    a = np.load(batch/'cellZones'/zone/f'{name}_proc_{rank}.npy')
                    b = np.load(reference/batch.name/f'{name}_proc_{rank}.npy')
                    assert a.flags.f_contiguous
                    np.testing.assert_array_equal(a, b[ids])
            for geometry in ['cellCentres', 'cellVolumes']:
                a = np.load(mapping.parent/f'{geometry}_proc_{rank}.npy')
                b = np.load(reference/mapping.parents[2].name/f'{geometry}_proc_{rank}.npy')
                np.testing.assert_array_equal(a, b[ids])
        if parallel:
            assert empty_rank
        for selection in ['left right empty', '"le.*"']:
            output = root+('All' if 'right' in selection else 'Subset')
            driver(case, root, output, f'cellZones ({selection});', dtype,
                   extra=export_object('zoneReexport', output+'Zones', f'cellZones ({selection});', dtype))
            process()
            target = case/'postProcessing'/output/'0.1'
            for batch in source.glob('batch_*'):
                for rank in range(4 if parallel else 1):
                    selected = ['left', 'right', 'empty'] if 'right' in selection else ['left']
                    ids = np.concatenate([np.load(source/'geometry_000000'/'cellZones'/z/f'cellIds_proc_{rank}.npy') for z in selected])
                    for name in FIELDS:
                        actual = np.load(target/batch.name/f'{name}_proc_{rank}.npy')
                        expected = np.full_like(actual, -99)
                        original = np.load(reference/batch.name/f'{name}_proc_{rank}.npy')
                        expected[ids] = original[ids]
                        np.testing.assert_array_equal(actual, expected)
                        for zone in selected:
                            relative = Path(batch.name)/'cellZones'/zone/f'{name}_proc_{rank}.npy'
                            reexport = case/'postProcessing'/(output+'Zones')/'0.1'/relative
                            assert reexport.read_bytes() == (source/relative).read_bytes()
        driver(case, root, 'bad', 'cellZones (left overlap);')
        process('Imported cellZones overlap')
        driver(case, root, 'bad', '')
        process('require explicit cellZones')
        driver(case, root, 'bad', 'cellZones ();')
        process('cellZones must not be empty')
        driver(case, root, 'bad', 'cellZones (missing);')
        process('Unmatched cellZones selection')
        driver(case, full, 'bad', 'cellZones (left);')
        process('cannot be imported as zones')
        # Exercise existing registered fields: a full import precedes a zone
        # importer in the driver's function list; outside values must survive.
        extra = f'''patch {{ type numpyToFoam; libs (numpyFunctionObjects);
            inputDir "postProcessing/{root}"; segment "0.1"; fields ({' '.join(FIELDS)});
            cellZones (left); executeControl timeStep; }}'''
        driver(case, full, root+'Existing', '', dtype, extra=extra)
        process()
        for original in reference.glob('batch_*/*_proc_*.npy'):
            actual = case/'postProcessing'/(root+'Existing')/'0.1'/original.relative_to(reference)
            np.testing.assert_array_equal(np.load(actual), np.load(original))
        if dtype == 'float64':
            # Corrupt one file at a time and ensure failures are explicit.
            driver(case, root, 'bad', 'cellZones (left);')
            mapping = source/'geometry_000000/cellZones/left/cellIds_proc_0.npy'
            payload = source/'batch_000000/cellZones/left/q_proc_0.npy'
            backups = {f: f.read_bytes() for f in [mapping, payload]}
            metadata = source/'geometry_000000/mesh_proc_0'
            metadata_text = metadata.read_text()
            metadata.write_text(re.sub(r'nCells\s+\d+', 'nCells 99999', metadata_text))
            process('mesh size or mapping revision mismatch')
            metadata.write_text(metadata_text)
            state = source/'batch_000000/state'
            state_text = state.read_text()
            state.write_text(state_text.replace('volScalarField', 'volVectorField'))
            process('Field class mismatch')
            state.write_text(state_text)
            ids = np.load(mapping)
            for values, message in [(ids[::-1], 'addressing mismatch'),
                                    (np.array([-1], dtype='<i8'), 'Invalid or truncated cell ID'),
                                    (np.array([2**40], dtype='<i8'), 'Invalid or truncated cell ID'),
                                    (ids.astype('<f8'), 'Expected one-dimensional int64')]:
                save(mapping, values)
                process(message)
                mapping.write_bytes(backups[mapping])
            original = np.load(payload)
            np.save(payload, np.ascontiguousarray(original))
            process('not Fortran-order')
            save(payload, original.astype('<i8'))
            process('Expected floating-point field')
            save(payload, original[:-1])
            process('Entity count mismatch')
            # Truncate first snapshot so preflight cannot consume it.
            payload.write_bytes(backups[payload][:130])
            process('Unexpected end of NumPy payload')
            payload.unlink()
            process('Cannot open NumPy file')
            payload.write_bytes(backups[payload])
            segment = source/'segmentInfo'
            data = segment.read_text()
            segment.write_text(re.sub(r'nProcs\s+\d+', 'nProcs 99', data))
            process('processor count mismatch')
            segment.write_text(data)
            # A later restart overrides times; trailing uncommitted values
            # must not become snapshots.
            restart = source.parent/'restart'
            shutil.copytree(source, restart)
            for array in restart.glob('batch_*/cellZones/left/*_proc_*.npy'):
                values = np.load(array)
                save(array, np.concatenate([values+7, values[..., :1]+999], axis=-1))
            driver(case, root, root+'Restart', 'cellZones (left);', segments='segments ("0.1" restart);')
            process()
            for batch in source.glob('batch_*'):
                for rank in range(4 if parallel else 1):
                    ids = np.load(source/f'geometry_000000/cellZones/left/cellIds_proc_{rank}.npy')
                    actual = np.load(case/f'postProcessing/{root}Restart/0.1/{batch.name}/q_proc_{rank}.npy')
                    orig = np.load(reference/batch.name/f'q_proc_{rank}.npy')
                    np.testing.assert_array_equal(actual[ids], orig[ids]+7)
                    assert actual.shape[-1] == orig.shape[-1]
    print('Cell-zone tests passed:', 'parallel' if parallel else 'serial', flush=True)


def lifecycle(case, build):
    build.mkdir()
    (build/'Make').mkdir()
    shutil.copy(HERE/'zoneLifecycle.C', build)
    (build/'Make/files').write_text('zoneLifecycle.C\nEXE = ' + str(build/'zoneLifecycle') + '\n')
    (build/'Make/options').write_text(
        'EXE_INC = -I' + str(HERE.parent/'src/functionObjects/foamToNumpy') +
        ' -I$(LIB_SRC)/finiteVolume/lnInclude -I$(LIB_SRC)/meshTools/lnInclude\n'
        'EXE_LIBS = -L$(FOAM_USER_LIBBIN) -lnumpyFunctionObjects -lfiniteVolume -lmeshTools\n')
    run(['wmake'], build)
    run([build/'zoneLifecycle'], case)
    segment = case/'postProcessing/lifecycle/0'
    batches = sorted(segment.glob('batch_*'))
    assert len(batches) == 3
    maps = []
    for batch in batches:
        state = (batch/'state').read_text()
        revision = int(re.search(r'meshRevision\s+(\d+)', state)[1])
        assert int(re.search(r'count\s+(\d+)', state)[1]) == 1
        mapping = segment/f'geometry_{revision:06d}/cellZones/left/cellIds_proc_0.npy'
        maps.append(np.load(mapping))
    assert len(maps[0]) == len(maps[1]) and not np.array_equal(maps[0], maps[1])
    np.testing.assert_array_equal(maps[1], maps[2])
    print('Cell-zone lifecycle and old-time tests passed', flush=True)


def main():
    with tempfile.TemporaryDirectory(prefix='numpy-cell-zones-') as tmp:
        case = Path(tmp)/'cavity'
        setup(case)
        lifecycle(case, Path(tmp)/"lifecycle-build")
        check(case, False)
        # Keep mesh and native fields; isolate serial output from parallel.
        shutil.rmtree(case/'postProcessing')
        check(case, True)

if __name__ == '__main__':
    main()
