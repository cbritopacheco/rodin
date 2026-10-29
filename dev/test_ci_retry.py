"""Recovery policy tests with an offline GitHub API fixture."""

import copy
import unittest
from pathlib import Path

import ci_retry


class RecoveryTest(unittest.TestCase):
    def setUp(self):
        self.run = {"id": 8, "name": "Tests", "event": "push", "status": "completed",
                    "conclusion": "failure", "head_branch": "develop", "head_sha": "abc",
                    "head_repository": {"full_name": "owner/repo"}, "run_attempt": 1}
        self.jobs = [{"id": 12, "name": "MPI / Slow / NP3", "conclusion": "failure"}]
        self.logs = "The runner has received a shutdown signal. The operation was canceled."
        self.calls = []
        self.recovery = ci_retry.Recovery("owner/repo", {"workflow_run": {"id": 8}}, self.api)

    def api(self, endpoint, method="GET", raw=False):
        self.calls.append((endpoint, method))
        if method == "POST":
            return ""
        if "/git/ref/" in endpoint:
            return {"object": {"sha": "abc"}}
        if "/jobs?" in endpoint:
            return [{"jobs": self.jobs[:1]}, {"jobs": self.jobs[1:]}]
        if endpoint.endswith("/logs"):
            return self.logs
        return copy.deepcopy(self.run)

    def posted(self):
        return [endpoint for endpoint, method in self.calls if method == "POST"]

    def test_interrupted_lane_is_retried_with_its_exact_id(self):
        self.assertIn("single retry", self.recovery.execute())
        self.assertEqual(self.posted(), ["repos/owner/repo/actions/jobs/12/rerun"])

    def test_numerical_and_resource_failures_are_not_retried(self):
        for error in ("[ FAILED ] Stokes", "runtime error: invalid pointer", "error: compilation",
                      "AddressSanitizer", "Out of memory", "Killed process 123",
                      "checksum mismatch", "***Timeout", "Assertion x failed"):
            with self.subTest(error=error):
                self.calls.clear()
                self.logs = error + "\nThe runner has received a shutdown signal"
                self.recovery.execute()
                self.assertFalse(self.posted())

    def test_shutdown_evidence_is_required(self):
        self.logs = "The operation was canceled."
        self.recovery.execute()
        self.assertFalse(self.posted())

    def test_repeat_failure_exhausts_retry_budget(self):
        self.jobs.append(dict(self.jobs[0], id=13))
        self.recovery.execute()
        self.assertFalse(self.posted())

    def test_only_one_job_is_retried_per_completion(self):
        self.jobs.append(dict(self.jobs[0], id=14, name="Other lane"))
        self.recovery.execute()
        self.assertEqual(len(self.posted()), 1)

    def test_untrusted_cancelled_active_and_superseded_runs_are_rejected(self):
        baseline = copy.deepcopy(self.run)
        for field, value in (("event", "pull_request"), ("head_branch", "feature"),
                             ("head_repository", {"full_name": "fork/repo"}),
                             ("conclusion", "cancelled"), ("status", "in_progress"),
                             ("head_sha", "old"), ("name", "Unknown")):
            with self.subTest(field=field):
                self.calls.clear()
                self.run = dict(baseline, **{field: value})
                self.recovery.execute()
                self.assertFalse(self.posted())

    def test_manual_rerun_during_inspection_prevents_race(self):
        original = self.api
        def racing(endpoint, **kwargs):
            result = original(endpoint, **kwargs)
            if endpoint.endswith("/logs"):
                self.run["run_attempt"] = 2
            return result
        self.recovery.api = racing
        self.recovery.execute()
        self.assertFalse(self.posted())

    def test_missing_logs_fail_closed(self):
        original = self.api
        def missing(endpoint, **kwargs):
            if endpoint.endswith("/logs"):
                raise RuntimeError("Logs unavailable")
            return original(endpoint, **kwargs)
        self.recovery.api = missing
        with self.assertRaisesRegex(RuntimeError, "Logs unavailable"):
            self.recovery.execute()
        self.assertFalse(self.posted())

    def test_unavailable_logs_do_not_block_other_verified_interruptions(self):
        for missing_name, eligible_name in [("A missing", "Z eligible"), ("Z missing", "A eligible")]:
            with self.subTest(missing=missing_name):
                self.calls.clear()
                self.jobs = [dict(self.jobs[0], id=12, name=missing_name),
                             dict(self.jobs[0], id=14, name=eligible_name)]
                original = self.api
                def missing(endpoint, **kwargs):
                    if endpoint.endswith("/jobs/12/logs"):
                        raise RuntimeError("Logs unavailable")
                    return original(endpoint, **kwargs)
                recovery = ci_retry.Recovery("owner/repo", {"workflow_run": {"id": 8}}, missing)
                self.assertIn("single retry", recovery.execute())
                self.assertEqual(self.posted(), ["repos/owner/repo/actions/jobs/14/rerun"])
                decisions = {job['id']: state for job, state, _ in recovery.rows}
                self.assertEqual(decisions[12], "Inspection failed")
                self.assertEqual(decisions[14], "Retry requested")
                self.assertIn("Logs unavailable; no retry requested", recovery.summary())

    def test_all_missing_logs_are_reported_before_inspection_fails(self):
        self.jobs.append(dict(self.jobs[0], id=14, name="Other lane"))
        original = self.api
        def missing(endpoint, **kwargs):
            if endpoint.endswith("/logs"):
                raise RuntimeError("Logs unavailable")
            return original(endpoint, **kwargs)
        self.recovery.api = missing
        with self.assertRaisesRegex(RuntimeError, "Logs unavailable"):
            self.recovery.execute()
        self.assertFalse(self.posted())
        self.assertEqual({job['id'] for job, _, _ in self.recovery.rows}, {12, 14})

    def test_summary_reports_every_failure_and_exact_log_links(self):
        self.jobs += [dict(self.jobs[0], id=14, name="Other | <lane>")]
        self.recovery.message = self.recovery.execute()
        summary = self.recovery.summary()
        self.assertIn("Retry requested", summary)
        self.assertIn("Eligible; deferred", summary)
        self.assertIn("Other &#124; &lt;lane&gt;", summary)
        self.assertIn("/actions/runs/8/job/12", summary)
        self.assertIn("abc", summary)
        self.assertIn("inspected attempt: 1", summary)

    def test_summary_explains_rejected_and_exhausted_retries(self):
        self.logs = "runtime error: invalid pointer"
        self.recovery.message = self.recovery.execute()
        self.assertIn("Failure evidence: runtime error:", self.recovery.summary())
        self.jobs.append(dict(self.jobs[0], id=13))
        self.recovery.rows.clear()
        self.recovery.message = self.recovery.execute()
        self.assertIn("Retry exhausted", self.recovery.summary())
        self.assertFalse(self.posted())

    def test_workflow_uses_trusted_code_and_narrow_permissions(self):
        import yaml
        root = Path(__file__).resolve().parents[1]
        workflow = yaml.safe_load((root / ".github/workflows/RetryInterruptedJobs.yml").read_text())
        self.assertEqual(workflow["permissions"], {"contents": "read", "checks": "read", "actions": "write"})
        steps = workflow["jobs"]["Recover"]["steps"]
        self.assertEqual(steps[0]["with"]["ref"], "${{ github.sha }}")
        self.assertNotIn("secrets.", str(steps))
        publication = yaml.safe_load((root / ".github/actions/publish-pages/action.yml").read_text())
        upload, deploy = publication["runs"]["steps"][1:3]
        self.assertEqual(upload["with"]["name"], deploy["with"]["artifact_name"])
        self.assertIn("github.run_attempt", upload["with"]["name"])
        for name in ("BenchmarkGroup", "BenchmarkConfiguration", "Benchmarks", "Documentation"):
            data = yaml.safe_load((root / f".github/workflows/{name}.yml").read_text())
            for job in data["jobs"].values():
                for step in job.get("steps", []):
                    if step.get("uses") == "actions/upload-artifact@v4":
                        self.assertTrue(step["with"]["overwrite"])


if __name__ == "__main__":
    unittest.main()
