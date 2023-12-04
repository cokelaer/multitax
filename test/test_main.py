import easydev
import os
import tempfile
import subprocess
import sys
from . import test_dir
from click.testing import CliRunner
from sequana_pipelines.multitax.main import main


sharedir = f"{test_dir}/data"
krakendb = f"{test_dir}/data/krakendb"


def test_standalone_subprocess():
    directory = tempfile.TemporaryDirectory()
    cmd = f"""sequana_multitax --input-directory {sharedir}
            --working-directory {directory.name} --force --databases {krakendb} """
    subprocess.call(cmd.split())


def test_standalone_script():
    directory = tempfile.TemporaryDirectory()

    runner = CliRunner()
    results = runner.invoke(main, ["--input-directory", sharedir, "--working-directory", directory.name, "--force",
"--databases", krakendb])
    assert results.exit_code == 0

    # 2 databases
    runner = CliRunner()
    results = runner.invoke(main, ["--input-directory", sharedir, "--working-directory", directory.name, "--force",
"--databases", krakendb, krakendb])
    assert results.exit_code == 0


def test_version():
    cmd = "sequana_multitax --version"
    subprocess.call(cmd.split())


