# script that tests all the python scripts in the test/python directory


import unittest
import subprocess
from pathlib import Path


this_file_path = Path(__file__).parent


class testALL(unittest.TestCase):
    def test_1(self):
        cmd = ["python3", f"{this_file_path}/main.py"]
        output = subprocess.run(cmd)
        self.assertEqual(output.returncode, 0)

    def test_2(self):
        cmd = ["python3", f"{this_file_path}/rrt.py"]
        output = subprocess.run(cmd)
        self.assertEqual(output.returncode, 0)

    def test_3(self):
        cmd = ["python3", f"{this_file_path}/rrt_free.py"]
        output = subprocess.run(cmd)
        self.assertEqual(output.returncode, 0)


if __name__ == "__main__":
    unittest.main()
