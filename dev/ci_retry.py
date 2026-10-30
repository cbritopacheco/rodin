#!/usr/bin/env python3
"""Retry one interrupted CI job once, after its trusted push run completes."""

import json
import os
import re
import subprocess
from pathlib import Path
from urllib.parse import quote

from ci_summary import append, cell


INTERRUPTION = re.compile(r"the runner has received a shutdown signal|the hosted runner.*lost communication", re.I)
FAILURE = re.compile(
    r"error:|fatal error|\[\s*FAILED\s*\]|\*\*\*Failed|\*\*\*Timeout|"
    r"AddressSanitizer|UndefinedBehaviorSanitizer|runtime error:|"
    r"out of memory|oom.kill|killed process|segmentation fault|"
    r"assertion.*failed|tests failed|benchmark.*failed|checksum mismatch|"
    r"Process completed with exit code (?!143\b)[1-9][0-9]*\b|Traceback \(most recent call last\)",
    re.I,
)


class Recovery:
    """Validate origin, inspect attempt history, then request one bounded rerun.

    Job names identify lanes across attempts. Only the latest occurrence can
    qualify, and any earlier occurrence exhausts its automatic retry budget.
    Logs and API failures are never interpreted as infrastructure evidence.
    """

    def __init__(self, repository, event, api=None):
        self.repository = repository
        self.event = event
        self.api = api or self.github
        self.prefix = f"repos/{repository}"
        self.rows = []
        self.run = None
        self.message = "Inspection did not complete; no retry is confirmed."

    @staticmethod
    def github(endpoint, method="GET", raw=False):
        command = ["gh", "api", "--method", method, endpoint]
        if raw:
            command.append("--allow-escape-sequences")
        if not raw:
            command += ["--paginate", "--slurp"] if "/jobs?" in endpoint else []
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return result.stdout if raw else json.loads(result.stdout)

    def current(self, run):
        if (run.get("name") not in {"Tests", "Benchmarks"}
                or run.get("event") != "push"
                or run.get("status") != "completed"
                or run.get("conclusion") != "failure"
                or run.get("head_branch") not in {"develop", "master"}
                or (run.get("head_repository") or {}).get("full_name") != self.repository):
            return False
        ref = quote(run["head_branch"], safe="")
        branch = self.api(f"{self.prefix}/git/ref/heads/{ref}")
        return branch["object"]["sha"] == run["head_sha"]

    @staticmethod
    def interrupted(job, logs):
        # SIGTERM (128 + 15) accompanies runner shutdown; the explicit runner
        # marker is still mandatory and every other failure veto remains.
        return (job.get("conclusion") in {"failure", "cancelled"}
                and bool(INTERRUPTION.search(logs))
                and not FAILURE.search(logs))

    def execute(self):
        run_id = int(self.event["workflow_run"]["id"])
        # Read current server state: the triggering event can already be stale.
        run = self.api(f"{self.prefix}/actions/runs/{run_id}")
        self.run = run
        if not self.current(run):
            return "No retry: run is successful, cancelled, untrusted, active or superseded."
        pages = self.api(f"{self.prefix}/actions/runs/{run_id}/jobs?filter=all&per_page=100")
        jobs = [job for page in pages for job in page["jobs"]]
        histories = {}
        for job in jobs:
            histories.setdefault(job["name"], []).append(job)
        candidates = []
        inspection_errors = []
        for name, history in sorted(histories.items()):
            job = max(history, key=lambda item: item["id"])
            if job.get("conclusion") not in {"failure", "cancelled", "timed_out"}:
                if len(history) > 1 and job.get("conclusion") == "success":
                    self.rows.append((job, "Passed after rerun", "No further retry required"))
                continue
            if len(history) != 1:
                self.rows.append((job, "Retry exhausted", "An earlier attempt already exists"))
                continue
            if job.get("conclusion") == "timed_out":
                self.rows.append((job, "No retry", "Job timeout requires investigation"))
                continue
            try:
                logs = self.api(f"{self.prefix}/actions/jobs/{job['id']}/logs", raw=True)
            except Exception as error:
                self.rows.append((job, "Inspection failed", "Logs unavailable; no retry requested"))
                inspection_errors.append(error)
                continue
            failure = FAILURE.search(logs)
            if failure:
                self.rows.append((job, "No retry", "Failure evidence: " + failure.group()))
            elif not self.interrupted(job, logs):
                self.rows.append((job, "No retry", "Insufficient runner-interruption evidence"))
            else:
                self.rows.append((job, "Eligible; deferred", "One job can be retried per completed attempt"))
                candidates.append((job, len(self.rows) - 1))
        if candidates:
            job, row = candidates[0]
            # Recheck after reading potentially lengthy logs; never revive an
            # old commit after a new push or an intervening manual rerun.
            latest = self.api(f"{self.prefix}/actions/runs/{run_id}")
            if latest["run_attempt"] != run["run_attempt"] or not self.current(latest):
                self.rows = [(job, "No retry", "Run changed during inspection")
                             if state == "Eligible; deferred" else (job, state, reason)
                             for job, state, reason in self.rows]
                return "No retry: run changed during inspection."
            self.rows[row] = (job, "Request failed", "No retry confirmed unless GitHub accepts the request")
            self.api(f"{self.prefix}/actions/jobs/{job['id']}/rerun", method="POST", raw=True)
            self.rows[row] = (job, "Retry requested", "Runner interruption; single retry budget consumed")
            return (f"Requested the single retry for job {job['id']} ({job['name']}) at "
                    f"{run['head_sha']}. GitHub also reruns its dependents. "
                    "All numerical and completeness gates remain mandatory.")
        # An unavailable log never qualifies its lane for recovery, but must
        # not discard independently verified interruptions in other lanes. If
        # nothing can be retried, retain the inspection failure as an error.
        if inspection_errors:
            raise inspection_errors[0]
        return "No eligible interrupted job remains; failures require investigation."

    def summary(self):
        server = os.environ.get("GITHUB_SERVER_URL", "https://github.com")
        base = f"{server}/{self.repository}/actions/runs/{int(self.event['workflow_run']['id'])}"
        rows = ["## CI failure inspection and retry decisions", "", cell(self.message), "",
                f"[Source workflow, attempts and job logs]({base})", ""]
        if self.run:
            rows += [f"Commit: {cell(self.run['head_sha'])}; inspected attempt: {self.run['run_attempt']}.", ""]
        if self.rows:
            rows += ["| Job / logs | Failed step | Decision | Reason |", "| --- | --- | --- | --- |"]
            # Keep the summary bounded; full job logs remain on the source run.
            for job, state, reason in self.rows[:200]:
                failed = [step['name'] for step in job.get('steps', [])
                          if step.get('conclusion') in {'failure', 'cancelled', 'timed_out'}]
                rows.append(f"| [{cell(job['name'])}]({base}/job/{job['id']}) | "
                            f"{cell(', '.join(failed) or 'See job logs')} | {cell(state)} | {cell(reason)} |")
            if len(self.rows) > 200:
                rows += ["", "Additional failures are visible in the source workflow."]
        rows += ["", "A requested retry is not a successful result. Numerical, coverage and pooling gates remain mandatory."]
        return "\n".join(rows)


def main():
    event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
    recovery = Recovery(os.environ["GITHUB_REPOSITORY"], event)
    try:
        recovery.message = recovery.execute()
        print(recovery.message)
        if any(state == 'Retry requested' for _, state, _ in recovery.rows):
            # requested workflow_run events do not fire for reruns. Refresh the
            # API-backed badges explicitly after accepting an automatic retry.
            result = subprocess.run(['gh', 'workflow', 'run', 'CIBadges.yml', '--repo',
                                     os.environ['GITHUB_REPOSITORY'], '-f',
                                     'branch=' + recovery.run['head_branch']],
                                    capture_output=True, text=True)
            if result.returncode:
                print('Badge refresh request unavailable; the next run event will refresh status.')
    finally:
        # Keep inspection failures visible even if fetching logs or rerunning
        # fails. Never turn reporting failures into successful CI recovery.
        append(recovery.summary())


if __name__ == "__main__":
    main()
