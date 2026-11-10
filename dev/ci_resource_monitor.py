#!/usr/bin/env python3
"""Emit periodic Linux memory readings without inspecting process arguments."""

import json
from pathlib import Path
import sys
import time


def snapshot(root=Path('/')):
    result = {'time_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
    for name in ('memory.current', 'memory.max', 'memory.events'):
        path = root / 'sys/fs/cgroup' / name
        try:
            result[name] = path.read_text().strip()
        except OSError:
            result[name] = None
    try:
        result['host_memory'] = {line.split(':', 1)[0]: line.split(':', 1)[1].strip()
                                 for line in (root / 'proc/meminfo').read_text().splitlines()
                                 if line.startswith(('MemTotal:', 'MemAvailable:', 'SwapFree:'))}
    except OSError:
        result['host_memory'] = None
    return result


def main():
    output = Path(sys.argv[1])
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
