"""Frozen, serial point-count-matched witness comparison; no production toggles."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import statistics
import subprocess
import time


class Comparison:
    def __init__(self, root, witness):
        self.repo = Path(__file__).resolve().parents[2]
        self.root = Path(root).resolve()
        self.witness = Path(witness).resolve()
        self.build = self.repo / 'build-p1-3d-clang19'
        self.env = dict(os.environ, OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4',
                        VECLIB_MAXIMUM_THREADS='1', OMP_MAX_ACTIVE_LEVELS='1')
        self.active = None

    def stop(self, signum, frame):
        if self.active is not None and self.active.poll() is None:
            os.killpg(self.active.pid, signal.SIGTERM)
            try:
                self.active.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(self.active.pid, signal.SIGKILL)
                self.active.wait()
        raise KeyboardInterrupt(f'Stopped by signal {signum}')

    @staticmethod
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def save(self, name, value):
        temporary = self.root / (name + '.new')
        temporary.write_text(json.dumps(value, indent=2))
        temporary.replace(self.root / name)

    def execute(self, command, folder, env, timeout=600, memory_limit=8 * 1024**3):
        started = time.monotonic()
        with (folder / 'run.log').open('w') as log:
            process = subprocess.Popen(['/usr/bin/time', '-l', '/usr/bin/nice', '-n', '19', *command],
                cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            self.active = process
            reason = None
            peak = 0
            while process.poll() is None:
                time.sleep(2)
                snapshot = subprocess.run(['ps', '-axo', 'pid,pgid,rss'], capture_output=True,
                                          text=True, check=True)
                rss = sum(int(fields[2]) * 1024 for line in snapshot.stdout.splitlines()[1:]
                          if len(fields := line.split()) == 3 and int(fields[1]) == process.pid)
                peak = max(peak, rss)
                if time.monotonic() - started > timeout:
                    reason = 'timeout'
                elif rss > memory_limit:
                    reason = 'memory-limit'
                if reason:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                    break
            code = process.wait()
            self.active = None
        return dict(returncode=code, failure=reason, wall_seconds=time.monotonic()-started,
                    sampled_peak_rss=peak)

    def prepare(self):
        self.root.mkdir(parents=True, exist_ok=False)
        for name in ('examples', 'include/Rodin/Adaptation/SWIFT', 'sets'):
            (self.root / name).mkdir(parents=True)
        shutil.copy2(self.witness, self.root / 'Witness')
        self.witness = self.root / 'Witness'
        for name in ('Options.h', 'LobedSphereLevelSet.h', 'Reconstruction.cpp'):
            shutil.copy2(self.repo / 'examples/Adaptation/SWIFT' / name, self.root / 'examples' / name)
        source = self.repo / 'src/Rodin/Adaptation'
        # Freeze all adaptation headers, including quoted relative includes.
        shutil.copytree(source, self.root / 'include/Rodin/Adaptation', dirs_exist_ok=True)
        shutil.copy2(self.repo / 'src/Rodin/Adaptation.h', self.root / 'include/Rodin/Adaptation.h')
        path = self.root / 'include/Rodin/Adaptation/SWIFT/QualityLattice.h'
        text = path.read_text().replace('#include <cassert>',
            '#include <cassert>\n#include <cstdlib>\n#include <fstream>\n#include <stdexcept>\n#include <Rodin/QF/GaussLobatto.h>')
        text = text.replace('Real getWeight(size_t) const override { return m_weight; }',
            'Real getWeight(size_t i) const override\n'
            '      { return m_referenceWeights.empty() ? m_weight : m_referenceWeights[i]; }')
        text = text.replace('      Real m_weight;', '      std::vector<Real> m_referenceWeights;\n      Real m_weight;')
        marker = '        m_weight = volume / m_points.size();'
        assert text.count(marker) == 1
        injection = '''        // Experimental control reproduces the former Lobatto sampling policy.
        const auto& lobatto = QF::GaussLobatto::get(geometry, std::max<size_t>(2, (m+4)/2));
        m_points.clear();
        for (size_t q = 0; q < lobatto.getSize(); ++q)
        {
          m_points.push_back(lobatto.getPoint(q));
          m_referenceWeights.push_back(lobatto.getWeight(q));
        }
        if (const char* filename = std::getenv("SWIFT_COVERING_POINTS"))
        {
          m_referenceWeights.clear();
          std::ifstream input(filename);
          size_t fileDimension = 0, count = 0;
          input >> fileDimension >> count;
          if (!input || fileDimension != dimension || count != m_points.size())
            throw std::runtime_error("Covering witness dimension/count mismatch");
          for (auto& point : m_points)
          {
            for (size_t axis = 0; axis < dimension; ++axis)
              input >> point[axis];
          }
          input >> std::ws;
          if (!input || !input.eof())
            throw std::runtime_error("Invalid covering witness file");
        }
'''
        path.write_text(text.replace(marker, injection + marker))
        example = self.root / 'examples/Reconstruction.cpp'
        text = example.read_text().replace('#include <Rodin/Adaptation.h>',
            '#include <Rodin/Adaptation/SWIFT/Problem.h>\n#include <Rodin/Adaptation.h>')
        text = text.replace('  try\n  {', '''  try
  {
    if (argc == 2 && std::string(argv[1]) == "--dump-quality-points")
    {
      std::cout << std::setprecision(17);
      for (const size_t dimension : {2, 3})
      {
        for (const size_t degree : {1, 2})
        {
          const auto geometry = dimension == 2 ? Polytope::Type::Triangle : Polytope::Type::Tetrahedron;
          const auto& rule = QF::GaussLobatto::get(geometry, degree == 1 ? 3 : 10);
          std::cout << dimension << ' ' << degree << ' ' << rule.getSize();
          for (size_t q = 0; q < rule.getSize(); ++q)
          {
            for (size_t axis = 0; axis < dimension; ++axis)
              std::cout << ' ' << rule.getPoint(q)[axis];
          }
          std::cout << '\\n';
        }
      }
      return 0;
    }''')
        marker = '    // Apply the solution to a separate mesh, including all higher-order geometry nodes.'
        assert text.count(marker) == 1
        text = text.replace('    if (!report.qualityBudgetSatisfied)\n      return 1;', '')
        text = text.replace(marker, '''    std::cout << "solve_setup_seconds=" << report.tSetup
      << " assembly_seconds=" << report.tAssembly << " solve_seconds=" << report.tSolve
      << " outer_search_seconds=" << report.tLineSearch
      << " inner_assembly_seconds=" << report.tInnerAssembly
      << " inner_solve_seconds=" << report.tInnerSolve
      << " inner_search_seconds=" << report.tInnerLineSearch << '\\n';
    // Independent boundary-inclusive audit; never used by the optimizer.
    const auto auditStart = std::chrono::steady_clock::now();
    const auto derivative = Jacobian(u.getSolution());
    Real auditJ = std::numeric_limits<Real>::infinity(), auditQ = 1;
    size_t violations = 0, samples = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto& lattice = Adaptation::SWIFT::QualityLattice::get(cell->getGeometry(), 20);
      // This lattice must remain independent of the experimental file override.
      for (size_t q = 0; q < lattice.getSize(); ++q)
      {
        const Point point(*cell, lattice.getPoint(q));
        Adaptation::CellDeformation state(dimension);
        state.setDisplacementGradient(derivative.getValue(point));
        const Real j = state.getJacobian();
        const Real Q = state.isAdmissible() ? state.getRelativeDistortion()
          : std::numeric_limits<Real>::infinity();
        auditJ = std::min(auditJ, j);
        auditQ = std::max(auditQ, Q);
        violations += j < options.parameters.model.jacobian - 1e-10
          || Q > options.parameters.model.distortion + 1e-8;
        ++samples;
      }
    }
    std::cout << "audit_min_j=" << auditJ << " audit_max_Q=" << auditQ
      << " audit_violations=" << violations << " audit_samples=" << samples
      << " audit_seconds=" << std::chrono::duration<double>(
          std::chrono::steady_clock::now()-auditStart).count() << '\\n';
    return report.qualityBudgetSatisfied ? 0 : 1;

''' + marker).replace('#include <filesystem>', '#include <filesystem>\n#include <chrono>')
        example.write_text(text)
        # A separately named, unmodified lattice supplies the independent audit.
        original = (source / 'SWIFT/QualityLattice.h').read_text()
        original = original.replace('QUALITYLATTICE_H', 'AUDITLATTICE_H').replace('QualityLattice', 'AuditLattice')
        (self.root / 'include/Rodin/Adaptation/SWIFT/AuditLattice.h').write_text(original)
        example.write_text(example.read_text().replace('#include <chrono>',
            '#include <chrono>\n#include <Rodin/Adaptation/SWIFT/AuditLattice.h>')
            .replace('SWIFT::QualityLattice::get(cell->getGeometry(), 20)',
                     'SWIFT::AuditLattice::get(cell->getGeometry(), 20)'))
        manifest = dict(planned=64, threads=4, warmups=0, repeats=0,
            point_counts={'triangle': {'1': 7, '2': 91}, 'tetrahedron': {'1': 15, '2': 820}},
            hinge=[10, 100], degrees=[1, 2], policies=['lobatto', 'covering'],
            cases=[dict(dimension=d, n=n, lobes=l) for d, ns, ls in
                [(2, [8, 16], [0, 6]), (3, [4, 6], [0, 4])] for n in ns for l in ls],
            witness_binary=str(self.witness), witness_hash=self.digest(self.witness),
            source_hashes={str(p.relative_to(self.root)): self.digest(p)
                for base in ['include', 'examples'] for p in (self.root/base).rglob('*') if p.is_file()},
            note='Offline equal-count covering generation; no production alternate path. '
                 'Final independent subdivision-20 audit. Concurrent timings not isolated.')
        self.save('manifest.json', manifest)
        self.save('results.json', [])

    def compile(self):
        commands = []
        for degree in (1, 2):
            folder = self.build / 'examples/Adaptation/SWIFT'
            cmake = folder / f'CMakeFiles/SWIFT_ReconstructionP{degree}.dir'
            flags = dict(line.split(' = ', 1) for line in (cmake/'flags.make').read_text().splitlines() if ' = ' in line)
            obj, binary = self.root/f'ComparisonP{degree}.o', self.root/f'ComparisonP{degree}'
            command = ['/usr/bin/nice', '-n', '19', '/opt/local/bin/mpicxx-mpich-clang19', '-I'+str(self.root/'include')]
            for key in ('CXX_DEFINES', 'CXX_INCLUDES', 'CXX_FLAGS'):
                command += shlex.split(flags[key])
            command += ['-c', str(self.root/'examples/Reconstruction.cpp'), '-o', str(obj)]
            commands.append(command)
            subprocess.run(command, cwd=self.repo, env=self.env, check=True)
            link = shlex.split((cmake/'link.txt').read_text())
            link = [str(obj) if arg.endswith('.cpp.o') else arg for arg in link]
            link[link.index('-o')+1] = str(binary)
            commands.append(link)
            subprocess.run(link, cwd=folder, env=self.env, check=True)
        self.save('build-commands.json', commands)
        output = subprocess.run([str(self.root/'ComparisonP1'), '--dump-quality-points'],
            capture_output=True, text=True, check=True, env=self.env)
        rules = {}
        for line in output.stdout.splitlines():
            values = line.split()
            dimension, degree, count = map(int, values[:3])
            coordinates = list(map(float, values[3:]))
            assert len(coordinates) == dimension * count
            rules[f'{dimension}:{degree}'] = [coordinates[i:i+dimension] for i in range(0, len(coordinates), dimension)]
        self.save('lobatto-points.json', rules)

    def generate(self, dimension, degree):
        geometry = 'triangle' if dimension == 2 else 'tetrahedron'
        points = json.loads((self.root/'lobatto-points.json').read_text())[f'{dimension}:{degree}']
        folder = self.root/'sets'/f'{geometry}-p{degree}'
        folder.mkdir()
        output = folder/'covering.json'
        command = [str(self.witness), str(len(points)), '--geometry', geometry, '--output', str(output)]
        self.save(str(folder.relative_to(self.root)/'command.json'), command)
        status = self.execute(command, folder, self.env)
        if status['returncode'] or status['failure'] or not output.exists():
            raise RuntimeError(f'Offline generation failed: {folder}: {status}')
        data = json.loads(output.read_text())
        assert len(data['points']) == len(points) and len(set(map(tuple, data['points']))) == len(points)
        for p in data['points']:
            assert len(p) == dimension and all(math.isfinite(x) and x >= -1e-12 for x in p)
            assert sum(p) <= 1+1e-12
        assert data['covering_radius'] <= data['initial_radius'] + 1e-10
        coordinates = folder/'points.txt'
        coordinates.write_text(f'{dimension} {len(points)}\n' + '\n'.join(
            ' '.join(format(x, '.17g') for x in p) for p in data['points']) + '\n')
        status.update(count=len(points), coordinates=str(coordinates), hash=self.digest(coordinates),
                      initial_radius=data['initial_radius'], covering_radius=data['covering_radius'],
                      search=data['search'])
        self.save(str(folder.relative_to(self.root)/'generation.json'), status)
        reference = folder/'lobatto-input.json'
        reference.write_text(json.dumps(dict(geometry=geometry, points=points)))
        audit_folder = folder/'lobatto-coverage'
        audit_folder.mkdir()
        coverage_command = [str(self.witness), '--evaluate', str(reference), '--output', str(audit_folder/'coverage.json')]
        coverage_status = self.execute(coverage_command, audit_folder, self.env)
        self.save(str(audit_folder.relative_to(self.root)/'command.json'), coverage_command)
        self.save(str(audit_folder.relative_to(self.root)/'status.json'), coverage_status)
        if coverage_status['returncode'] or coverage_status['failure']:
            raise RuntimeError(f'Lobatto coverage evaluation failed: {audit_folder}: {coverage_status}')
        return coordinates

    def run(self):
        manifest = json.loads((self.root/'manifest.json').read_text())
        assert not (self.root/'started.json').exists(), 'Campaign already started; no implicit resume'
        self.save('started.json', dict(pid=os.getpid(), time=time.time()))
        for path, digest in manifest['source_hashes'].items():
            assert self.digest(self.root/path) == digest, f'Frozen source changed: {path}'
        self.witness = Path(manifest['witness_binary'])
        assert self.digest(self.witness) == manifest['witness_hash'], 'Frozen Witness changed'
        manifest['binaries'] = {str(d): self.digest(self.root/f'ComparisonP{d}') for d in (1, 2)}
        self.save('manifest.json', manifest)
        rows, sets = [], {}
        # Compute all reference sets offline, before starting any reconstruction.
        for dimension in (2, 3):
            for degree in (1, 2):
                key = (dimension, degree)
                sets[key] = self.generate(*key)
                print('offline-set-ready', key, flush=True)
        for case in manifest['cases']:
            for degree in (1, 2):
                key = (case['dimension'], degree)
                for mu in (10, 100):
                    for policy in manifest['policies']:
                        name = f"d{case['dimension']}-n{case['n']}-l{case['lobes']}-p{degree}-mu{mu}-{policy}"
                        folder = self.root/name
                        folder.mkdir()
                        command = [str(self.root/f'ComparisonP{degree}'),
                            *[f'--{k}={v}' for k, v in case.items()], '--trace',
                            '--model-fit=1', '--model-distribution-deviatoric=1e-4',
                            '--model-distribution-divergence=1e-2', f'--model-hinge={mu}',
                            '--linear-solver=mumps', '--linear-threads=4']
                        env = dict(self.env)
                        env.pop('SWIFT_COVERING_POINTS', None)
                        if policy == 'covering':
                            env['SWIFT_COVERING_POINTS'] = str(sets[key])
                        self.save(str(folder.relative_to(self.root)/'command.json'), dict(argv=command,
                            environment={k:v for k,v in env.items() if k.startswith(('OMP_', 'OPENBLAS_', 'SWIFT_'))}))
                        load = subprocess.run(['ps', '-axo', 'pid,pcpu,rss,comm'], capture_output=True, text=True, check=True)
                        (folder/'cpu-start.txt').write_text(load.stdout)
                        row = dict(case, degree=degree, mu=mu, policy=policy,
                                   **self.execute(command, folder, env))
                        text = (folder/'run.log').read_text()
                        for line in text.splitlines():
                            if line.startswith(('exit=', 'audit_min_j=', 'solve_setup_seconds=')):
                                row.update(dict(re.findall(r'(\w+)=([^\s]+)', line)))
                        row['accepted'] = [dict(re.findall(r'(\w+)=([^\s]+)', line)) for line in text.splitlines()
                                           if 'swift geometry:' in line and 'phase=accepted' in line]
                        inner = [int(a['inner_last']) for a in row['accepted'] if 'inner_last' in a]
                        row['median_inner'] = statistics.median(inner) if inner else None
                        timing = re.search(r'([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys', text)
                        if timing:
                            row['process_real'], row['cpu_user'], row['cpu_sys'] = map(float, timing.groups())
                        (folder/'cpu-end.txt').write_text(subprocess.run(['ps', '-axo', 'pid,pcpu,rss,comm'],
                            capture_output=True, text=True, check=True).stdout)
                        rows.append(row)
                        self.save('results.json', rows)
                        print(name, row['returncode'], row.get('exit'), row.get('D_inf'), flush=True)
        self.save('completed.json', dict(planned=64, completed=len(rows)))

    def main(self, prepare_only=False, run_only=False):
        signal.signal(signal.SIGTERM, self.stop)
        signal.signal(signal.SIGINT, self.stop)
        try:
            if not run_only:
                self.prepare()
                self.compile()
            if not prepare_only:
                self.run()
        except BaseException as error:
            if self.root.exists():
                self.save('blocked.json', dict(error=repr(error), time=time.time()))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--witness', required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--prepare-only', action='store_true')
    mode.add_argument('--run-only', action='store_true')
    options = parser.parse_args()
    Comparison(options.root, options.witness).main(options.prepare_only, options.run_only)
