"""Checks the offline sample parser, not benchmark performance or PDE accuracy."""

import io
import unittest

from summarize_xctrace import XctraceProfile


class XctraceProfileTest(unittest.TestCase):
  def setUp(self):
    self.xml = """<trace-query-result><node><schema name="time-profile"/>
      <row><weight id="w">1000000</weight></row>
      <row><thread><process id="p" fmt="RodinWorkload (123)"/></thread>
        <process ref="p"/><weight ref="w"/><tagged-backtrace id="s">
          <frame id="f" name="kernel"/><frame ref="f"/>
          <frame id="b" name="Rodin::Bench::run"/>
        </tagged-backtrace></row>
      <row><process ref="p"/><weight ref="w"/><tagged-backtrace ref="s"/></row>
      <row><process ref="p"/><weight id="v">2000000</weight>
        <tagged-backtrace><frame name="fixture"/><frame ref="b"/></tagged-backtrace>
      </row>
      <row><process id="other" fmt="Other (456)"/><weight ref="v"/>
        <tagged-backtrace ref="s"/></row>
      </node></trace-query-result>"""

  def test_weight_defined_without_stack_and_process_reference(self):
    report = XctraceProfile(io.StringIO(self.xml),
      required=("Rodin::Bench::run", "kernel"), process="RodinWorkload",
      categories=(("kernel", "kernel"), ("wrapper", "::run"))).summarize()
    self.assertEqual(report["samples_with_weight"], 5)
    self.assertEqual(report["selected_samples"], 2)
    self.assertEqual(report["total_sampled_ms"], 7)
    self.assertEqual(report["selected_sampled_ms"], 2)
    self.assertEqual(report["category_patterns"], {"kernel": "kernel", "wrapper": "::run"})
    self.assertEqual(report["leaf"][0]["sampled_ms"], 2)
    self.assertEqual(report["inclusive"][0]["percent_selected"], 100)
    self.assertEqual(sum(x["percent_selected"] for x in report["categories"]), 200)

  def test_fixture_is_excluded_by_kernel_selector(self):
    report = XctraceProfile(io.StringIO(self.xml),
      required=("kernel",)).summarize()
    self.assertEqual(report["selected_sampled_ms"], 4)
    self.assertNotIn("fixture", [x["name"] for x in report["inclusive"]])

  def test_original_kernel_is_excluded_explicitly(self):
    row = """<row><weight ref="w"/><tagged-backtrace><frame ref="f"/>
      <frame ref="b"/><frame name="ReferenceIntegral::integrate"/>
      </tagged-backtrace></row>"""
    xml = self.xml.replace("</node>", row + "</node>")
    report = XctraceProfile(io.StringIO(xml), required=("kernel",),
      excluded=("ReferenceIntegral",)).summarize()
    self.assertEqual(report["selected_sampled_ms"], 4)
    self.assertNotIn("ReferenceIntegral::integrate",
      [x["name"] for x in report["inclusive"]])

  def test_unknown_reference_is_rejected(self):
    with self.assertRaises(KeyError):
      XctraceProfile(io.StringIO(self.xml.replace('ref="w"', 'ref="missing"', 1)),
        required=("kernel",)).summarize()

  def test_empty_selection_is_rejected(self):
    with self.assertRaises(ValueError):
      XctraceProfile(io.StringIO(self.xml), required=("absent",)).summarize()

  def test_duplicate_category_names_are_rejected(self):
    with self.assertRaises(ValueError):
      XctraceProfile(io.StringIO(self.xml),
        categories=(("same", "kernel"), ("same", "::run")))

  def test_zero_weight_is_rejected(self):
    with self.assertRaises(ValueError):
      XctraceProfile(io.StringIO(self.xml.replace("1000000", "0").replace("2000000", "0")),
        required=("kernel",)).summarize()

  def test_other_schema_is_rejected(self):
    with self.assertRaises(ValueError):
      XctraceProfile(io.StringIO(self.xml.replace('name="time-profile"',
        'name="other"')), required=("kernel",)).summarize()


if __name__ == "__main__":
  unittest.main()
