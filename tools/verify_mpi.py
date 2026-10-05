#!/usr/bin/env python3
"""Local MPI integration checks; requires a Release build, mpicc, mpirun and Python 3."""
import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / 'src/TomoStar.Cli/bin/Release/net10.0/tomostar'


def run(args, env, expected=0):
    result = subprocess.run([str(x) for x in args], env=env, cwd=ROOT,
                            capture_output=True, text=True, timeout=120)
    if result.returncode != expected:
        raise AssertionError(f'{args}: exit {result.returncode}\n{result.stdout}\n{result.stderr}')
    return result


def volume(path):
    blob = path.read_bytes()
    header_length = struct.unpack_from('<i', blob, 8)[0]
    offset = ((12 + header_length + 63) // 64) * 64
    return blob[offset:]


with tempfile.TemporaryDirectory(prefix='tomostar-mpi-') as scratch:
    scratch = Path(scratch)
    native = scratch / 'native'
    env = os.environ.copy()
    run(['bash', ROOT / 'native/mpi/build.sh', native], env)
    env['LD_LIBRARY_PATH'] = str(native) + ':' + env.get('LD_LIBRARY_PATH', '')
    data = scratch / 'data'
    run([EXE, 'synth', '--out', data, '--stations', '4', '--events', '8', '--no-opencl',
         '--set', 'Synthetic.HorizontalSpacingKm=12', '--set', 'Synthetic.VerticalSpacingKm=6'], env)
    common = [data, '--grid', data / 'grid.json', '--model', data / 'model1d.txt',
              '--no-opencl', '--threads', '1', '--name', 'check', '--quiet']
    serial = scratch / 'serial'
    run([EXE, 'invert', *common, '--iterations', '2', '--out', serial], env)
    reference = {p.name: volume(p) for p in (serial / 'volumes').glob('*.qvol')}
    assert reference
    for ranks in (1, 2, 4):
        output = scratch / f'mpi-{ranks}'
        launch = ['mpirun', '-np', str(ranks), EXE, '--mpi']
        run([*launch, 'invert', *common, '--iterations', '2', '--out', output], env)
        actual = {p.name: volume(p) for p in (output / 'volumes').glob('*.qvol')}
        assert reference == actual, f'Volume mismatch with {ranks} ranks'
        for name in ('iterations.csv', 'residuals.csv', 'hypocentres.csv'):
            assert (serial / name).read_bytes() == (output / name).read_bytes(), (ranks, name)
        # Standalone table dispatch, independently of the ray worker path.
        relocated = scratch / f'relocated-{ranks}'
        run([*launch, 'relocate', *common, '--out', relocated], env)
        if ranks == 1:
            relocation_reference = (relocated / 'events_relocated.csv').read_bytes()
        else:
            assert relocation_reference == (relocated / 'events_relocated.csv').read_bytes()
        print(f'PASS: {ranks} MPI ranks; inversion and relocation identical')
    # CLI errors stop workers; accidental multi-rank execution without --mpi is rejected.
    run(['mpirun', '-np', '2', EXE, '--mpi', 'unknown-command'], env, expected=2)
    run(['mpirun', '-np', '2', EXE, 'version'], env, expected=2)
    # Exercise multiple dispatches through the root-only script interpreter.
    script = scratch / 'study.tomo'
    script.write_text(f'invert {data} --grid {data}/grid.json --model {data}/model1d.txt '
                      f'--no-opencl --threads 1 --iterations 1 --out {scratch}/script-output\n'
                      f'invert {data} --grid {data}/grid.json --model {data}/model1d.txt '
                      f'--no-opencl --threads 1 --iterations 1 --out {scratch}/script-output-2\n')
    run(['mpirun', '-np', '2', EXE, '--mpi', 'run', script], env)
    # Straight rays and extended domains use distinct forward branches. Compare both.
    cropped = json.loads((data / 'grid.json').read_text())
    cropped['MinLon'] += 0.1
    cropped['MaxLon'] -= 0.1
    crop_path = scratch / 'cropped.json'
    crop_path.write_text(json.dumps(cropped))
    for label, extra in (('straight', ['--straight']),
                         ('extended', ['--grid', crop_path, '--set', 'OutsideData=WhenRaysCross'])):
        outputs = []
        for ranks in (0, 2):
            output = scratch / f'{label}-{ranks}'
            launch = [EXE] if ranks == 0 else ['mpirun', '-np', str(ranks), EXE, '--mpi']
            run([*launch, 'invert', *common, *extra, '--iterations', '1', '--out', output], env)
            outputs.append({p.name: volume(p) for p in (output / 'volumes').glob('*.qvol')})
        assert outputs[0] == outputs[1], label
        print(f'PASS: {label} forward matches serial execution')
    # A worker-side filesystem error must reach rank zero without a collective deadlock.
    failed = scratch / 'failed'
    (failed / 'work').mkdir(parents=True)
    (failed / 'work/mpi-rank-1').write_text('not a directory')
    result = run(['mpirun', '-np', '2', EXE, '--mpi', 'invert', *common,
                  '--iterations', '1', '--out', failed], env, expected=1)
    assert 'rank 1:' in result.stderr, result.stderr
    # More ranks than stations exercises workers with no arrivals.
    few = scratch / 'few'
    few.mkdir()
    (few / 'stations.csv').write_text('\n'.join((data / 'stations.csv').read_text().splitlines()[:3]) + '\n')
    run(['mpirun', '-np', '4', EXE, '--mpi', 'invert', data, '--stations', few / 'stations.csv',
         '--grid', data / 'grid.json', '--model', data / 'model1d.txt', '--no-opencl',
         '--threads', '1', '--iterations', '1', '--out', scratch / 'few-output'], env)
    print('PASS: worker shutdown, launcher guard, scripts, worker errors and empty partitions')
