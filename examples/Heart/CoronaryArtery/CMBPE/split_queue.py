#!/usr/bin/env python3
"""Splits a running campaign queue in two without touching the running jobs.

Run it AFTER killing only the scheduler (queue.py); the mpirun jobs it
started keep running. Writes, next to jobs.tsv:
  jobs_main.tsv  every job except the selected ones; each job still running is
                 replaced by a placeholder that holds its cores until the run
                 ends, checks that its last step was reached and writes .done;
  jobs_sel.tsv   the selected jobs (default: names starting with pra2400/ or dt/).
    python3 split_queue.py ~/Results/cmbpe/coronary [--select pra2400/,dt/]
"""
import argparse, os, re, shlex

ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--select", default="pra2400/,dt/")
a = ap.parse_args()
out = os.path.abspath(os.path.expanduser(a.out))
sel = tuple(s for s in a.select.split(",") if s)

jobs = []
for line in open(os.path.join(out, "jobs.tsv")):
    if line.strip():
        jobs.append(line.rstrip("\n").split("\t", 4))

# directories in which a CoronaryArtery process is running now
busy = set()
for pid in filter(str.isdigit, os.listdir("/proc")):
    try:
        cmd = open("/proc/%s/cmdline" % pid, "rb").read()
        if b"CoronaryArtery" in cmd and b"split_queue" not in cmd:
            busy.add(os.path.realpath(os.readlink("/proc/%s/cwd" % pid)))
    except OSError:
        pass

def finished(d, cmd):
    """The run reached its last step (it ended while no scheduler was watching)."""
    m = re.search(r"-coronary_nsteps (\d+)", cmd)
    log = os.path.join(d, "log.txt")
    return bool(m) and os.path.isfile(log) and \
        ("Step %s " % m.group(1)) in open(log, errors="ignore").read()

main, chosen = [], []
for name, np, d, after, cmd in jobs:
    done = os.path.exists(os.path.join(d, ".done"))
    if not done and os.path.realpath(d) not in busy and finished(d, cmd):
        open(os.path.join(d, ".done"), "w").close()
        done = True
        print("marked done  %s (finished while unscheduled)" % name)
    if name.startswith(sel):
        if os.path.realpath(d) in busy:
            raise SystemExit("%s is running: stop it first or leave it in the main queue" % name)
        chosen.append((name, np, d, after, cmd))
    elif os.path.realpath(d) in busy and not done:
        m = re.search(r"-coronary_nsteps (\d+)", cmd)
        last = m.group(1) if m else ""
        q = shlex.quote(os.path.realpath(d))
        wait = ("busy() { for p in $(pgrep -f CoronaryArtery); do "
                "[ \"$(readlink /proc/$p/cwd 2>/dev/null)\" = %s ] && return 0; done; return 1; }; "
                "while busy; do sleep 60; done; "
                "grep -q 'Step %s ' %s/log.txt && touch %s/.done" % (q, last, q, q))
        main.append((name, np, os.path.join(d, ".wait"), after, wait))
        print("placeholder  %s (running)" % name)
    else:
        main.append((name, np, d, after, cmd))

names_main = {j[0] for j in main}
for j in chosen:
    for dep in (j[3].split(",") if j[3] not in ("", "-") else []):
        if dep in names_main:
            raise SystemExit("%s depends on %s, which stays in the main queue" % (j[0], dep))

for fn, js in (("jobs_main.tsv", main), ("jobs_sel.tsv", chosen)):
    with open(os.path.join(out, fn), "w") as f:
        for j in js:
            f.write("\t".join(j) + "\n")
    print("%-13s %d jobs" % (fn, len(js)))
for j in chosen:
    print("  selected", j[0])
