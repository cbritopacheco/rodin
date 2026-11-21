#!/usr/bin/env python3
"""Compile and execute component consumers against an existing Rodin installation."""
import argparse
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--components', nargs='+', default=['Geometry'])
    parser.add_argument('--cmake-arg', action='append', default=[])
    args = parser.parse_args()
    source = Path(__file__).resolve().parent / 'components'
    with tempfile.TemporaryDirectory(prefix='rodin-installed-components-') as directory:
        for component in args.components:
            build = Path(directory) / component.replace('::', '-')
            command = ['cmake', '-S', str(source), '-B', str(build),
                       '-DRodin_ROOT=' + str(args.prefix.resolve()),
                       '-DRODIN_COMPONENT=' + component, *args.cmake_arg]
            if component != 'Kokkos::PETSc':
                command.append('-DCMAKE_DISABLE_FIND_PACKAGE_PETSc=TRUE')
            subprocess.run(command, check=True, timeout=120)
            subprocess.run(['cmake', '--build', str(build), '--parallel', '2'],
                           check=True, timeout=300)
            subprocess.run([str(build / 'consumer')], cwd=build, check=True, timeout=60)


if __name__ == '__main__':
    main()
