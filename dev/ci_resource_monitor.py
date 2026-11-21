#!/usr/bin/env python3
"""Sample Linux runner resources and publish explicitly scoped CI statistics."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time

SAMPLE_INTERVAL = 15  # Seconds; sampled peaks can miss shorter excursions.


def snapshot(root=Path('/')):
    result = {'time_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'monotonic': time.monotonic()}
    for name in ('memory.current', 'memory.peak', 'memory.max', 'memory.events'):
        try:
            result[name] = (root / 'sys/fs/cgroup' / name).read_text().strip()
        except OSError:
            result[name] = None
    try:
        result['host_memory'] = {line.split(':', 1)[0]: line.split(':', 1)[1].strip()
                                 for line in (root / 'proc/meminfo').read_text().splitlines()
                                 if line.startswith(('MemTotal:', 'MemAvailable:', 'SwapFree:'))}
    except OSError:
        result['host_memory'] = None
    try:
        rows = (root / 'proc/stat').read_text().splitlines()
        # Guest time is already included in user/nice; do not count it twice.
        ticks = [int(value) for value in rows[0].split()[1:9]]
        result['cpu'] = {'total': sum(ticks), 'idle': ticks[3] + ticks[4],
                         'cores': sum(line.startswith('cpu') and line[3:4].isdigit() for line in rows)}
    except (OSError, ValueError, IndexError):
        result['cpu'] = None
    return result


def number(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def events(sample):
    try:
        return {key: int(value) for key, value in
                (line.split() for line in sample['memory.events'].splitlines())}
    except (AttributeError, KeyError, TypeError, ValueError):
        return {}


def statistics(samples):
    """Describe host samples; never sum process RSS or invent missing counters."""
    result = {'samples': len(samples), 'interval_seconds': SAMPLE_INTERVAL,
              'elapsed_seconds': None, 'cpu_percent': None, 'cpu_cores': None,
              'host_total_bytes': None, 'peak_host_used_bytes': None, 'minimum_available_bytes': None,
              'peak_cgroup_bytes': None, 'oom_events': None, 'oom_kills': None}
    if not samples:
        return result
    first, last = samples[0], samples[-1]
    result['elapsed_seconds'] = max(0, last['monotonic'] - first['monotonic'])
    cpus = [sample['cpu'] for sample in samples if sample.get('cpu')]
    if len(cpus) >= 2:
        total = cpus[-1]['total'] - cpus[0]['total']
        idle = cpus[-1]['idle'] - cpus[0]['idle']
        if total > 0 and 0 <= idle <= total:
            result['cpu_percent'] = 100 * (total - idle) / total
            result['cpu_cores'] = cpus[-1]['cores']
    used, available, peaks = [], [], []
    for sample in samples:
        memory = sample.get('host_memory') or {}
        try:
            total = int(memory['MemTotal'].split()[0]) * 1024
            free = int(memory['MemAvailable'].split()[0]) * 1024
            if 0 <= free <= total:
                result['host_total_bytes'] = total
                used.append(total - free)
                available.append(free)
        except (KeyError, ValueError, IndexError):
            pass
        peak = number(sample.get('memory.peak'))
        if peak is not None:
            peaks.append(peak)
    result['peak_host_used_bytes'] = max(used, default=None)
    result['minimum_available_bytes'] = min(available, default=None)
    result['peak_cgroup_bytes'] = max(peaks, default=None)
    before, after = events(first), events(last)
    for key, target in [('oom', 'oom_events'), ('oom_kill', 'oom_kills')]:
        if key in before and key in after and after[key] >= before[key]:
            result[target] = after[key] - before[key]
    return result


def record(directory):
    with (directory / 'samples.jsonl').open('a') as stream:
        stream.write(json.dumps(snapshot()) + '\n')
        stream.flush()


def monitor(directory):
    # A sentinel avoids killing a possibly recycled PID. The runner also
    # cleans up this child on job termination through its inherited job token.
    while not (directory / 'stop').exists():
        record(directory)
        deadline = time.monotonic() + SAMPLE_INTERVAL
        while time.monotonic() < deadline and not (directory / 'stop').exists():
            time.sleep(0.2)


def start(directory):
    directory.mkdir(parents=True, exist_ok=True)
    # Refuse to mix measurements if this step is accidentally launched twice.
    with (directory / 'samples.jsonl').open('x') as stream:
        stream.write(json.dumps(snapshot()) + '\n')
    subprocess.Popen([sys.executable, str(Path(__file__).resolve()), 'monitor', str(directory)],
                     stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                     stderr=subprocess.DEVNULL, start_new_session=True)
    print(f'Resource sampling started: {directory}; interval {SAMPLE_INTERVAL}s.')


def finish(directory):
    from ci_summary import resources
    if directory.exists():
        (directory / 'stop').touch()
        record(directory)  # Include CPU counters and OOM events at the end.
    samples = []
    path = directory / 'samples.jsonl'
    if path.exists():
        for line in path.read_text().splitlines():
            try:
                samples.append(json.loads(line))
            except json.JSONDecodeError:
                pass  # A truncated final record is diagnostic, not a test failure.
    # A simultaneous sampler write can arrive after the final record.
    samples.sort(key=lambda sample: sample['monotonic'])
    result = statistics(samples)
    if directory.exists():
        (directory / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print('CI resources: ' + json.dumps(result, sort_keys=True), flush=True)
    resources(result)


def main():
    operation = sys.argv[1]
    if operation in {'start', 'monitor', 'finish'}:
        directory = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(os.environ['RUNNER_TEMP']) / 'rodin-resources'
        {'start': start, 'monitor': monitor, 'finish': finish}[operation](directory)
    else:
        # Existing benchmark console diagnostics retain their original interface.
        output = Path(operation)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('a') as stream:
            while True:
                line = 'CI memory: ' + json.dumps(snapshot(), sort_keys=True)
                print(line, flush=True)
                stream.write(line + '\n')
                stream.flush()
                time.sleep(60)


if __name__ == '__main__':
    main()
