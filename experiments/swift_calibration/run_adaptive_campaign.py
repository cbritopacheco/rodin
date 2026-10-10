"""Serial, frozen-binary adaptive affine-hinge calibration."""

import argparse
import csv
import fcntl
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone


class Campaign:
    DISTRIBUTION = (0, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1)
    HINGE = (0.1, 1, 10, 100, 1000)
    JACOBIAN_WEIGHT = (0.1, 1, 10)
    THREADS = 4
    TIMEOUT_SECONDS = 1800
    SAMPLE_SECONDS = 2
    MAX_RSS_BYTES = 8 * 1024**3
    CSV_FIELDS = (
        'attempt', 'dimension', 'n', 'lobes', 'degree', 'dev', 'div', 'mu', 'kj',
        'exit', 'returncode', 'D_inf', 'C', 'energy', 'max_Q', 'min_j', 'outer',
        'inner', 'inner_max', 'inner_median', 'target_hit', 'quality_ok',
        'best_D_inf', 'first_hit_outer', 'first_hit_seconds', 'seconds',
        'peak_rss_bytes', 'cpu_flag', 'log',
    )

    def __init__(self, output):
        self.output = Path(output).resolve()
        self.stop = False
        self.child = None

    @classmethod
    def plan(cls):
        rows = []
        for dimension, resolutions, lobes in (
            (2, (8, 16, 32), (0, 4, 6, 8)),
            (3, (4, 6, 8), (0, 4, 6)),
        ):
            for n, lobe in itertools.product(resolutions, lobes):
                for dev, div, mu, kj in itertools.product(
                    cls.DISTRIBUTION, cls.DISTRIBUTION, cls.HINGE, cls.JACOBIAN_WEIGHT
                ):
                    for degree in (1, 2):
                        rows.append(dict(attempt=len(rows) + 1, dimension=dimension,
                                         n=n, lobes=lobe, degree=degree,
                                         dev=dev, div=div, mu=mu, kj=kj))
        return rows

    @staticmethod
    def digest(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    @staticmethod
    def fields(line):
        return dict(re.findall(r'(\w+)=([^\s]+)', line))

    def save(self, name, value):
        temporary = self.output / (name + '.tmp')
        temporary.write_text(json.dumps(value, indent=2))
        temporary.replace(self.output / name)

    def command(self, row):
        degree = row['degree']
        return [
            str(self.output / f'SWIFT_ReconstructionP{degree}'),
            f"--dimension={row['dimension']}", f"--n={row['n']}",
            f"--lobes={row['lobes']}", '--amp=0.05', '--R0=0.25', '--phase=0',
            '--cx=0.5', '--cy=0.5', '--cz=0.5',
            '--model-fit=1', f"--model-distribution-deviatoric={row['dev']}",
            f"--model-distribution-divergence={row['div']}",
            f"--model-hinge={row['mu']}", f"--model-jacobian-weight={row['kj']}",
            '--model-distortion-weight=1', '--model-jacobian=0.01',
            '--model-distortion=10', '--model-quality-guard=0.1',
            '--model-robust-scale=0', '--globalization-directional-newton=1',
            '--globalization-max-step-over-h=0', '--globalization-armijo=0.0001',
            '--convergence-tolerance-geometric=0', '--convergence-iterations-outer=30',
            '--convergence-iterations-inner=15', '--convergence-tolerance-inner-relative=0.001',
            '--convergence-tolerance-inner-absolute=1e-12',
            '--convergence-tolerance-linear-relative=1e-6', '--convergence-iterations-linear=1000',
            '--convergence-tolerance-energy=1e-8', '--convergence-tolerance-step=0',
            '--convergence-tolerance-step-over-h=0.0005', '--convergence-iterations-stagnation=5',
            '--convergence-iterations-backtracks=32', '--linear-solver=mumps',
            f'--linear-threads={self.THREADS}', '--quadrature-order=0',
            f'--quadrature-surface={8 if degree == 1 else 12}',
            f'--quadrature-volume={2 if degree == 1 else 8}',
            f'--sampling-subdivision={2 if degree == 1 else 16}',
            '--quadrature-validation=32', '--trace',
            '--trace-quality-witness',
            f'--output={self.output / "scratch" / "reconstruction"}',
        ]

    def prepare(self):
        repo = Path(__file__).resolve().parents[2]
        driver = (repo / 'examples/Adaptation/SWIFT/Reconstruction.cpp').read_text()
        sampling = (repo / 'src/Rodin/Adaptation/SWIFT/QualitySamples.h').read_text()
        if ('DirichletBC' in driver or 'QualityCovering::get' not in sampling
                or 'QF::GaussLobatto::get' in sampling):
            raise RuntimeError('Expected free boundaries and canonical closed-form quality coverings')
        self.output.mkdir(parents=True, exist_ok=False)
        (self.output / 'logs').mkdir()
        (self.output / 'scratch').mkdir()
        plan = self.plan()
        assert len(plan) == 30870
        binaries = {}
        for degree in (1, 2):
            name = f'SWIFT_ReconstructionP{degree}'
            source = repo / 'build-p1-3d-clang19/examples/Adaptation/SWIFT' / name
            destination = self.output / name
            shutil.copy2(source, destination)
            binaries[str(degree)] = dict(path=str(destination), sha256=self.digest(destination),
                dependencies=subprocess.check_output(['otool', '-L', str(destination)], text=True))
        snapshot = self.output / 'source'
        snapshot.mkdir()
        sources = {}
        for relative in ('src/Rodin/Adaptation/SWIFT', 'src/Rodin/QF',
                         'examples/Adaptation/SWIFT'):
            destination = snapshot / relative
            shutil.copytree(repo / relative, destination)
            for path in destination.rglob('*'):
                if path.is_file():
                    sources[str(path.relative_to(snapshot))] = self.digest(path)
        shutil.copy2(__file__, self.output / 'run.py')
        protocol = repo / 'experiments/swift_calibration/ADAPTIVE_HINGE_CAMPAIGN.md'
        shutil.copy2(protocol, self.output / protocol.name)
        (self.output / 'source.diff').write_bytes(subprocess.check_output(
            ['git', 'diff', 'HEAD'], cwd=repo))
        self.save('plan.json', plan)
        environment = dict(OMP_NUM_THREADS=str(self.THREADS),
                           OPENBLAS_NUM_THREADS=str(self.THREADS),
                           VECLIB_MAXIMUM_THREADS=str(self.THREADS),
                           OMP_MAX_ACTIVE_LEVELS='1')
        self.save('manifest.json', dict(
            planned=len(plan), planned_2d=17640, planned_3d=13230, repo=str(repo),
            commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
            binaries=binaries, sources=sources, runner_sha256=self.digest(self.output / 'run.py'),
            protocol_sha256=self.digest(protocol), environment=environment,
            order='dimension, n, lobes, dev, div, mu, kj, degree',
            maximum_concurrent_cases=1, warmups=0, repeats=0,
            boundary_conditions='none; free exterior boundary; no gauge',
            quality_witnesses=dict(policy='closed-form reference covering; homothetic on tetrahedra',
                subdivisions={'1': 2, '2': 16},
                simplex_counts={'triangle': {'1': 6, '2': 153},
                                'tetrahedron': {'1': 10, '2': 969}},
                shared_hinge_and_validation=True, supplemental_points=False,
                weights='Unchanged frozen constraint-specific adaptive weights'),
            timeout_seconds=self.TIMEOUT_SECONDS, maximum_rss_bytes=self.MAX_RSS_BYTES,
            cpu_sampling_seconds=self.SAMPLE_SECONDS,
            launch_prefix=['/usr/bin/time', '-l'],
            memory_measurement='Native time peak RSS; process-tree polling for live memory limit',
            outputs='Compressed full logs and per-case results; mesh outputs discarded after parsing',
        ))
        with (self.output / 'results.csv').open('x') as stream:
            csv.DictWriter(stream, fieldnames=self.CSV_FIELDS).writeheader()
        with (self.output / 'launch.log').open('x') as stream:
            process = subprocess.Popen(
                ['/usr/bin/caffeinate', '-i', '/usr/bin/nice', '-n', '10',
                 sys.executable, '-u', str(self.output / 'run.py'),
                 '--output', str(self.output), '--worker'],
                cwd=repo, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        self.save('launch.json', dict(pid=process.pid, planned=len(plan),
            launched=datetime.now(timezone.utc).isoformat()))
        print(json.dumps(dict(pid=process.pid, output=str(self.output), planned=len(plan))))

    def sample(self):
        rows = []
        output = subprocess.check_output(['ps', '-Ao', 'pid,ppid,pcpu,rss,comm'], text=True)
        for line in output.splitlines()[1:]:
            fields = line.split(None, 4)
            if len(fields) == 5:
                rows.append(dict(pid=int(fields[0]), ppid=int(fields[1]),
                                 cpu=float(fields[2]), rss=int(fields[3]) * 1024,
                                 command=fields[4]))
        return rows

    def external_work(self, rows, excluded):
        scientific = ('SWIFT_', 'KelvinBall', 'Reconstruction', 'Comparison',
                      'RodinConvergence', 'clang', 'cc1', 'g++', 'cc1plus')
        return [row for row in rows if row['pid'] not in excluded and row['cpu'] >= 50
                and any(name in Path(row['command']).name for name in scientific)]

    def request_stop(self, *_):
        self.stop = True
        if self.child and self.child.poll() is None:
            os.killpg(self.child.pid, signal.SIGTERM)

    def run(self):
        manifest = json.loads((self.output / 'manifest.json').read_text())
        repo = Path(manifest['repo'])
        with (repo / 'tmp/swift-campaign.lock').open('a+') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            signal.signal(signal.SIGTERM, self.request_stop)
            signal.signal(signal.SIGINT, self.request_stop)
            for binary in manifest['binaries'].values():
                if self.digest(binary['path']) != binary['sha256']:
                    raise RuntimeError('Frozen executable hash mismatch')
            environment = dict(os.environ, **manifest['environment'])
            with (self.output / 'results.csv').open() as stream:
                previous = list(csv.DictReader(stream))
            finished = {int(row['attempt']) for row in previous}
            if len(finished) != len(previous):
                raise RuntimeError('Duplicate completed campaign attempts')
            recorded = set()
            if (self.output / 'results.jsonl').exists():
                with (self.output / 'results.jsonl').open() as stream:
                    recorded = {int(json.loads(line)['attempt']) for line in stream if line.strip()}
            if recorded != finished:
                raise RuntimeError('CSV and JSONL campaign records disagree')
            completed = len(finished)
            hits = sum(row['target_hit'] == '1' for row in previous)
            for row in json.loads((self.output / 'plan.json').read_text()):
                if row['attempt'] in finished:
                    continue
                if self.stop:
                    break
                command = self.command(row)
                self.save('progress.json', dict(completed=completed, planned=manifest['planned'],
                    target_hits=hits, status='running', current=row, command=command))
                with (self.output / 'commands.jsonl').open('a') as stream:
                    stream.write(json.dumps(dict(**row, argv=command)) + '\n')
                    stream.flush()
                    os.fsync(stream.fileno())
                log = self.output / 'logs' / f"{row['attempt']:05d}.log"
                started = time.monotonic()
                cpu_samples = []
                peak_rss = 0
                resource_exit = None
                with log.open('w') as stream:
                    self.child = subprocess.Popen(['/usr/bin/time', '-l', *command], cwd=self.output / 'scratch',
                        env=environment, stdout=stream, stderr=subprocess.STDOUT,
                        start_new_session=True)
                    while self.child.poll() is None:
                        try:
                            self.child.wait(timeout=self.SAMPLE_SECONDS)
                        except subprocess.TimeoutExpired:
                            pass
                        processes = self.sample()
                        own = {self.child.pid}
                        while True:
                            descendants = {p['pid'] for p in processes if p['ppid'] in own}
                            if descendants <= own:
                                break
                            own |= descendants
                        peak_rss = max(peak_rss, sum(p['rss'] for p in processes if p['pid'] in own))
                        external = self.external_work(processes, own | {os.getpid(), os.getppid()})
                        if external:
                            cpu_samples.append(dict(seconds=time.monotonic() - started, processes=external))
                        if self.stop or time.monotonic() - started > self.TIMEOUT_SECONDS or peak_rss > self.MAX_RSS_BYTES:
                            resource_exit = ('interrupted' if self.stop else
                                'memory-limit' if peak_rss > self.MAX_RSS_BYTES else 'timeout')
                            if self.child.poll() is None:
                                os.killpg(self.child.pid, signal.SIGTERM)
                                try:
                                    self.child.wait(timeout=10)
                                except subprocess.TimeoutExpired:
                                    os.killpg(self.child.pid, signal.SIGKILL)
                            break
                    code = self.child.wait()
                result = dict(row, returncode=code, seconds=time.monotonic() - started,
                              peak_rss_bytes=peak_rss, cpu_flag=bool(cpu_samples))
                summaries, states = [], []
                with log.open() as stream:
                    for line in stream:
                        memory = re.search(r'(\d+)\s+maximum resident set size', line)
                        if memory:
                            result['peak_rss_bytes'] = max(result['peak_rss_bytes'], int(memory.group(1)))
                        if line.startswith('exit='):
                            summaries.append(self.fields(line))
                        elif 'swift geometry:' in line:
                            states.append(self.fields(line))
                if summaries:
                    result.update(summaries[-1])
                if resource_exit or not summaries:
                    result['exit'] = resource_exit or 'process-failure'
                    result['target_hit'] = '0'
                accepted = [s for s in states if s['phase'] in ('initial', 'accepted')]
                if accepted:
                    result['best_D_inf'] = min(float(s['geom_sup']) for s in accepted)
                    first_hit = next((s for s in accepted if
                        float(s['geom_sup']) <= float(s['geom_sup_target']) and
                        float(s['min_j']) > 0.01 and float(s['max_qrel']) < 10), None)
                    if first_hit:
                        result['first_hit_outer'] = first_hit['outer']
                        result['first_hit_seconds'] = first_hit['seconds']
                    inner = [int(s['inner_last']) for s in accepted if s['phase'] == 'accepted']
                    result['inner_median'] = statistics.median(inner) if inner else 0
                compressed = log.with_suffix('.log.gz')
                with log.open('rb') as source, gzip.open(compressed, 'wb') as destination:
                    shutil.copyfileobj(source, destination)
                log.unlink()
                result['log'] = str(compressed)
                with (self.output / 'results.csv').open('a') as stream:
                    csv.DictWriter(stream, fieldnames=self.CSV_FIELDS,
                                   extrasaction='ignore').writerow(result)
                    stream.flush()
                    os.fsync(stream.fileno())
                with (self.output / 'results.jsonl').open('a') as stream:
                    stream.write(json.dumps(dict(result, accepted=accepted, cpu_samples=cpu_samples)) + '\n')
                    stream.flush()
                    os.fsync(stream.fileno())
                if resource_exit or not summaries or code or result.get('target_hit') != '1':
                    with (self.output / 'failures.jsonl').open('a') as stream:
                        stream.write(json.dumps(result) + '\n')
                for path in (self.output / 'scratch').iterdir():
                    if path.is_file():
                        path.unlink()
                completed += 1
                hits += int(result.get('target_hit', '0') == '1')
                print(f"{completed}/{manifest['planned']} d={row['dimension']} n={row['n']} "
                      f"lobes={row['lobes']} P{row['degree']} exit={result.get('exit')} "
                      f"C={result.get('C', 'NA')}", flush=True)
                self.save('progress.json', dict(completed=completed, planned=manifest['planned'],
                    target_hits=hits, status='between-cases', last=result))
                if self.stop:
                    break
            name = 'stopped.json' if self.stop else 'completed.json'
            self.save(name, dict(completed=completed, planned=manifest['planned'],
                                target_hits=hits, time=datetime.now(timezone.utc).isoformat()))
            self.save('progress.json', dict(completed=completed, planned=manifest['planned'],
                target_hits=hits, status='stopped' if self.stop else 'complete'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    parser.add_argument('--launch', action='store_true')
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    campaign = Campaign(args.output)
    if args.worker:
        campaign.run()
    elif args.launch:
        campaign.prepare()
    else:
        plan = campaign.plan()
        print(json.dumps(dict(planned=len(plan), first=plan[0], last=plan[-1])))
