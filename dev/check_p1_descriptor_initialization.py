#!/usr/bin/env python3
"""Check the initialization-cost contract in an unstripped benchmark binary.

This opt-in diagnostic uses GNU/LLVM-compatible demangled nm output. It is not
a CTest timing gate. Build RodinReferenceElementBenchmarks first: its real and
complex P1 cases retain all 33 reference descriptors per scalar type. Guard
variables for valid basis descriptors imply runtime static initialization;
constant-initialized descriptors must have no such guards. Numerical behavior
is qualified separately by the benchmark guards and finite-element tests.
"""

import argparse
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--nm", default="nm")
    args = parser.parse_args()
    symbols = subprocess.check_output(
        [args.nm, "-C", str(args.binary)], text=True).splitlines()
    descriptors = [line for line in symbols
                   if "P1Element<" in line and "::getBasis(" in line
                   and "::basis" in line]
    guards = [line for line in descriptors if "guard variable for " in line]
    objects = [line for line in descriptors if "guard variable for " not in line]
    # Eight geometries have 1+2+3+4+4+8+5+6 vertices, for each scalar field.
    expected = 2 * (1 + 2 + 3 + 4 + 4 + 8 + 5 + 6)
    print(f"P1 basis objects: {len(objects)}; runtime guards: {len(guards)}")
    if len(objects) != expected:
        raise SystemExit("Expected all real/complex P1 descriptors; wrong or stripped binary")
    if guards:
        raise SystemExit("Valid P1 basis descriptors require runtime initialization")


if __name__ == "__main__":
    main()
