import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ST"))
from pipeline_utils import input_path, output_directory, select_library


class PipelineHelpersTest(unittest.TestCase):
    def test_input_validation_and_nested_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "input.h5ad"
            data.write_text("fixture")
            resolved = input_path(data)
            output_directory(root / "results" / "run1")
            self.assertEqual(resolved, data.resolve())
            self.assertEqual(data.read_text(), "fixture")
            with self.assertRaises(FileNotFoundError):
                input_path(root / "missing.h5ad")
            with self.assertRaises(FileNotFoundError):
                input_path(data, directory=True)

    def test_real_library_id_is_inferred(self):
        data = SimpleNamespace(uns={"spatial": {"sample_001": {}}})
        self.assertEqual(select_library(data), "sample_001")
        with self.assertRaises(ValueError):
            select_library(data, "spatial")

    def test_multilibrary_requires_explicit_selection(self):
        data = SimpleNamespace(uns={"spatial": {"one": {}, "two": {}}})
        with self.assertRaises(ValueError):
            select_library(data)
        self.assertEqual(select_library(data, "two"), "two")
        with self.assertRaises(ValueError):
            select_library(SimpleNamespace(uns={}))


if __name__ == "__main__":
    unittest.main()
