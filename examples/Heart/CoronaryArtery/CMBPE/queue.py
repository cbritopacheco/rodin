#!/usr/bin/env python3
"""Runs a list of jobs concurrently within a core budget.

Each line of the job file is tab-separated:

    name  np  workdir  after  command

`np` cores are reserved while the job runs; `after` is "-" or a comma-separated
list of job names that must have finished successfully first. A job runs
`bash -c command` inside `workdir`, with stdout and stderr in workdir/log.txt,
and leaves workdir/.done on success. Jobs whose .done exists are skipped, so
the queue can be interrupted and relaunched.

    python3 queue.py jobs.tsv --cores 60
"""
import argparse
import os
import subprocess
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs")
    ap.add_argument("--cores", type=int, default=os.cpu_count())
    ap.add_argument("--poll", type=float, default=5.0)
    a = ap.parse_args()

    jobs = []
    with open(a.jobs) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            name, np, workdir, after, cmd = line.split("\t", 4)
            deps = [] if after.strip() in ("", "-") else after.split(",")
            jobs.append(dict(name=name, np=int(np), dir=workdir, deps=deps, cmd=cmd))

    names = {j["name"] for j in jobs}
    for j in jobs:
        for d in j["deps"]:
            if d not in names:
                sys.exit("job %s depends on unknown job %s" % (j["name"], d))
        if j["np"] > a.cores:
            sys.exit("job %s needs %d cores > budget %d" % (j["name"], j["np"], a.cores))

    state = {}
    for j in jobs:
        state[j["name"]] = "done" if os.path.exists(os.path.join(j["dir"], ".done")) else "pending"
    skipped = sum(1 for s in state.values() if s == "done")
    if skipped:
        print("[queue] %d job(s) already done" % skipped, flush=True)

    running = {}
    free = a.cores
    t0 = time.time()

    def stamp():
        return "%7.1f min" % ((time.time() - t0) / 60.0)

    while True:
        # reap
        for name, (proc, job) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            free += job["np"]
            del running[name]
            if rc == 0:
                open(os.path.join(job["dir"], ".done"), "w").close()
                state[name] = "done"
                print("[queue %s] done   %s" % (stamp(), name), flush=True)
            else:
                state[name] = "failed"
                print("[queue %s] FAILED %s (exit %d, see %s/log.txt)"
                      % (stamp(), name, rc, job["dir"]), flush=True)

        # propagate failures to dependants
        for j in jobs:
            if state[j["name"]] == "pending" and any(state[d] == "failed" for d in j["deps"]):
                state[j["name"]] = "failed"
                print("[queue %s] skip   %s (a dependency failed)" % (stamp(), j["name"]), flush=True)

        # launch, in file order, every ready job that fits
        for j in jobs:
            if state[j["name"]] != "pending" or j["np"] > free:
                continue
            if not all(state[d] == "done" for d in j["deps"]):
                continue
            os.makedirs(j["dir"], exist_ok=True)
            log = open(os.path.join(j["dir"], "log.txt"), "w")
            proc = subprocess.Popen(["bash", "-c", j["cmd"]], cwd=j["dir"], stdout=log,
                                    stderr=subprocess.STDOUT)
            running[j["name"]] = (proc, j)
            state[j["name"]] = "running"
            free -= j["np"]
            print("[queue %s] start  %s (%d cores, %d free)" % (stamp(), j["name"], j["np"], free),
                  flush=True)

        pending = [n for n, s in state.items() if s in ("pending", "running")]
        if not pending:
            break
        time.sleep(a.poll)

    failed = [n for n, s in state.items() if s == "failed"]
    print("[queue %s] finished: %d done, %d failed" % (stamp(), sum(1 for s in state.values()
                                                                    if s == "done"), len(failed)))
    for n in failed:
        print("  failed:", n)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
