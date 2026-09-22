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
                "1",
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


def test_project_script_creates_nested_out_dir():
    # Parent dir does not exist yet - out_dir is nested inside it
    tmpdir = tempfile.mkdtemp()
    shutil.rmtree(tmpdir)
    nested_out_dir = os.path.join(tmpdir, "experiment1", "run3")

    try:
        # Run the script using uv run
        _ = subprocess.run(
            [
                "uv",
                "run",
                "python",
                "src/ipp/project.py",
                "-o",
                nested_out_dir,
                "-n",
                "1",
                "data/dummy_regions.bed",
                "mm39",
                "hg38",
                "data/dummy.pwaln.bin",
            ],
            capture_output=True,
            text=True,
            check=True,
        )

        # Check the nested output directory was created with files in it
        assert os.path.isdir(nested_out_dir)
        assert len(os.listdir(nested_out_dir)) > 0

    finally:
        # Clean up the temporary directory
        shutil.rmtree(tmpdir)
