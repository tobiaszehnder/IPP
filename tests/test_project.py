import os
import shutil
import subprocess
import tempfile


def test_cpp_backend_project_script_runs():
    # Create a temporary output directory
    tmpdir = tempfile.mkdtemp()

    try:
        # Run the script using uv run
        _ = subprocess.run(
            [
                "uv",
                "run",
                "python",
                "src/ipp/project.py",
                "-o",
                tmpdir,
                "-n",
                "10",
                "data/dummy_regions.bed",
                "mm39",
                "hg38",
                "data/dummy.pwaln.bin",
            ],
            capture_output=True,
            text=True,
            check=True,
        )

        # Check files were created in tmpdir
        assert os.path.isdir(tmpdir)
        assert len(os.listdir(tmpdir)) > 0

    finally:
        # Clean up the temporary directory
        shutil.rmtree(tmpdir)
