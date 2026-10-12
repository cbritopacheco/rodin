"""Regression checks for the repository's two-space record indentation."""

import os
from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parent.parent


class FormatConfigTests(unittest.TestCase):
    def check_format(self, expected):
        binary = os.environ.get("CLANG_FORMAT_BINARY", "clang-format")
        source = expected.replace("  ", "    ")
        for text in (source, expected):
            with self.subTest(idempotent=text == expected):
                result = subprocess.run(
                    [binary, "--style=file", "--assume-filename=src/FormatFixture.cpp"],
                    input=text, capture_output=True, text=True, cwd=ROOT,
                    check=True)
                self.assertEqual(result.stdout, expected)

    def test_plain_struct_and_method_body(self):
        self.check_format("""namespace Example
{
  struct Table
  {
    int value;
    int get() const
    {
      return value;
    }
  };
}
""")

    def test_nested_struct_and_access_labels(self):
        self.check_format("""namespace Example
{
  class Store
  {
    public:
    struct Table
    {
      int value;
    };

    protected:
    int capacity;

    private:
    int next;
  };
}
""")


if __name__ == "__main__":
    unittest.main()
