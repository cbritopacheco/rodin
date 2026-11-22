#!/usr/bin/env python3
"""Prepares the remaining coronary jobs to run from a frozen copy of the code.

  python3 resume_jobs.py OUT FROZEN [--drop tree/] [--long 10]

- marks .done every job whose log reached its last step (the scheduler was
  stopped before it could do so);
- rewrites the executable, mesh and script paths of every job to FROZEN
  (FROZEN/CoronaryArtery, FROZEN/meshes/, FROZEN/matched_secant.py), so that
  switching branches or rebuilding the repository no longer affects the runs;
- orders the remaining jobs by priority for the paper and drops the given
  prefixes; with --long N, adds healthy rmu and rest runs of N cycles
  (periodicity check, note S4);
- writes OUT/jobs_resume.tsv and prints what will run.
"""
import argparse, os, re

ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("frozen")
ap.add_argument("--drop", default="")
ap.add_argument("--long", type=int, default=0)
a = ap.parse_args()
out = os.path.abspath(os.path.expanduser(a.out))
F = os.path.abspath(os.path.expanduser(a.frozen))

jobs = [l.rstrip("\n").split("\t", 4) for l in open(os.path.join(out, "jobs.tsv")) if l.strip()]

def rewrite(cmd):
    cmd = re.sub(r"\S+/CFD/CoronaryArtery\b", F + "/CoronaryArtery", cmd)
    cmd = re.sub(r"\S+/CMBPE/meshes/", F + "/meshes/", cmd)
    cmd = re.sub(r"\S+/CMBPE/matched_secant\.py", F + "/matched_secant.py", cmd)
    return cmd

def finished(d, cmd):
    m = re.search(r"-coronary_nsteps (\d+)", cmd)
    log = os.path.join(d, "log.txt")
    return bool(m) and os.path.isfile(log) and ("Step %s " % m.group(1)) in open(log, errors="ignore").read()

drop = tuple(p for p in a.drop.split(",") if p)
prio = ["pra700/", "pra2400/", "mesh/", "dt/", "tree/"]
def rank(name):
    if name.startswith("long"):
        return 0.5          # longest runs first
    if name.endswith("/match"):
        return 2.5          # matched constants right after the base runs
    return next((i for i, p in enumerate(prio) if name.startswith(p)), len(prio))

# runs still alive must not be scheduled again
busy = set()
for pid in filter(str.isdigit, os.listdir("/proc")):
    try:
        if b"CoronaryArtery" in open("/proc/%s/cmdline" % pid, "rb").read():
            busy.add(os.path.realpath(os.readlink("/proc/%s/cwd" % pid)))
    except OSError:
        pass
alive = [n for n, _, d, _, _ in jobs if os.path.realpath(d) in busy]
if alive:
    raise SystemExit("still running, wait for them or stop them first: " + ", ".join(alive))

keep, done = [], []
for name, np, d, after, cmd in jobs:
    if os.path.exists(os.path.join(d, ".done")) or (not name.endswith("/match") and finished(d, cmd)):
        open(os.path.join(d, ".done"), "w").close() if os.path.isdir(d) else None
        done.append(name)
    if name.startswith(drop):
        continue
    keep.append([name, np, d, after, rewrite(cmd)])

names = {j[0] for j in keep}
for j in keep:   # a dependency that was dropped would block the job
    j[3] = ",".join(x for x in j[3].split(",") if x in names) if j[3] not in ("", "-") else "-"
    j[3] = j[3] or "-"

if a.long:
    base = {j[0]: j for j in keep}
    for c in ("rmu", "rest"):
        src = base.get("pra700/healthy/" + c)
        if src:
            m = re.search(r"-coronary_nsteps (\d+)", src[4])
            n = int(round(int(m.group(1)) * a.long / 3.0))
            d = os.path.join(out, "long%d" % a.long, "healthy", c)
            keep.append(["long%d/healthy/%s" % (a.long, c), src[1], d, "-",
                         src[4].replace(m.group(0), "-coronary_nsteps %d" % n)])

keep.sort(key=lambda j: rank(j[0]))
with open(os.path.join(out, "jobs_resume.tsv"), "w") as f:
    for j in keep:
        f.write("\t".join(j) + "\n")

left = [j for j in keep if not os.path.exists(os.path.join(j[2], ".done"))]
print("already done (%d): %s" % (len(done), ", ".join(sorted(done))))
print("to run (%d jobs, %d cores at most at once):" % (len(left), sum(int(j[1]) for j in left)))
for j in left:
    print("  %-34s np=%-3s after=%s" % (j[0], j[1], j[3]))
stale = [j[0] for j in keep if re.search(r"/(rodin[^/]*)/(build|examples)", j[4])]
if stale:
    print("WARNING: paths still point to a repository in:", ", ".join(stale))
