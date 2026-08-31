import importlib.util
import sys
import types
import unittest
from pathlib import Path


class FakeInterval:
    def __init__(self, lower_bound, upper_bound):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

    def __contains__(self, position):
        return self.lower_bound <= position <= self.upper_bound


# The coverage logic can be tested without installing the native pysam dependency.
sys.modules.setdefault("pysam", types.ModuleType("pysam"))
interval_module = types.ModuleType("interval")
interval_module.Interval = FakeInterval
sys.modules.setdefault("interval", interval_module)
spec = importlib.util.spec_from_file_location(
    "calculate_cnv", Path(__file__).parents[1] / "atacamp" / "cnv.py"
)
calculate_cnv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(calculate_cnv)


class FakeAlignment:
    def __init__(self, length=100):
        self.length = length
        self.calls = []

    def get_reference_length(self, chromosome):
        return self.length

    def count_coverage(self, chromosome, start, end):
        self.calls.append((chromosome, start, end))
        # One base of coverage per genomic position across the four bases.
        return ([1] * (end - start), [], [], [])


class CoverageWindowTests(unittest.TestCase):
    def setUp(self):
        self.finder = calculate_cnv.find_amp.__new__(calculate_cnv.find_amp)
        self.finder.interval = 10

    def test_custom_left_window_is_honored(self):
        alignment = FakeAlignment()
        self.assertEqual(self.finder.interval_cov(alignment, "chr1", 50, 20), 1)
        self.assertEqual(alignment.calls, [("chr1", 30, 50)])

    def test_windows_are_clipped_at_contig_boundaries(self):
        alignment = FakeAlignment()
        self.assertEqual(self.finder.interval_cov(alignment, "chr1", 5), 1)
        self.assertEqual(self.finder.interval_cov1(alignment, "chr1", 95), 1)
        self.assertEqual(alignment.calls, [("chr1", 0, 5), ("chr1", 95, 100)])

    def test_find_amplicon_extends_left_for_in_range_breakpoint(self):
        self.finder.cov = 0.5
        alignment = FakeAlignment(length=100)
        left, right = self.finder.find_amplicon(alignment, "chr1", 50)
        self.assertEqual((left, right), (0, 100))


if __name__ == "__main__":
    unittest.main()
