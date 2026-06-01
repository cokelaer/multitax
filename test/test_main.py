import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from sequana_pipelines.multitax.main import main

from . import test_dir

sharedir = f"{test_dir}/data"
krakendb = f"{test_dir}/data/krakendb"
paireddir = f"{test_dir}/data/paired"
simpledir = f"{test_dir}/data/simple"


def test_standalone_subprocess():
    directory = tempfile.TemporaryDirectory()
    cmd = f"""sequana_multitax --input-directory {sharedir}
            --working-directory {directory.name} --force --databases {krakendb} """
    subprocess.call(cmd.split())


def test_standalone_script():
    directory = tempfile.TemporaryDirectory()

    runner = CliRunner()
    results = runner.invoke(
        main, ["--input-directory", sharedir, "--working-directory", directory.name, "--force", "--databases", krakendb]
    )
    assert results.exit_code == 0

    # 2 databases
    runner = CliRunner()
    results = runner.invoke(
        main,
        [
            "--input-directory",
            sharedir,
            "--working-directory",
            directory.name,
            "--force",
            "--databases",
            krakendb,
            krakendb,
        ],
    )
    assert results.exit_code == 0


def test_version():
    cmd = "sequana_multitax --version"
    subprocess.call(cmd.split())


def test_paired_reads_single_db(tmp_path):
    runner = CliRunner()
    results = runner.invoke(
        main, ["--input-directory", paireddir, "--working-directory", str(tmp_path), "--force", "--databases", krakendb]
    )
    assert results.exit_code == 0


def test_paired_reads_multiple_dbs(tmp_path):
    runner = CliRunner()
    results = runner.invoke(
        main,
        [
            "--input-directory",
            paireddir,
            "--working-directory",
            str(tmp_path),
            "--force",
            "--databases",
            krakendb,
            krakendb,
        ],
    )
    assert results.exit_code == 0


def test_store_unclassified(tmp_path):
    runner = CliRunner()
    results = runner.invoke(
        main,
        [
            "--input-directory",
            sharedir,
            "--working-directory",
            str(tmp_path),
            "--force",
            "--databases",
            krakendb,
            "--store-unclassified",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((Path(str(tmp_path)) / "config.yaml").read_text())
    assert config["sequana_taxonomy"]["store_unclassified"] is True


def test_kraken_confidence(tmp_path):
    runner = CliRunner()
    results = runner.invoke(
        main,
        [
            "--input-directory",
            sharedir,
            "--working-directory",
            str(tmp_path),
            "--force",
            "--databases",
            krakendb,
            "--kraken-confidence",
            "0.5",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((Path(str(tmp_path)) / "config.yaml").read_text())
    assert config["sequana_taxonomy"]["confidence"] == 0.5


def test_blast_unclassified_config(tmp_path):
    runner = CliRunner()
    results = runner.invoke(
        main,
        [
            "--input-directory",
            sharedir,
            "--working-directory",
            str(tmp_path),
            "--force",
            "--databases",
            krakendb,
            "--do-blast-unclassified",
        ],
    )
    assert results.exit_code == 0
    config = yaml.safe_load((Path(str(tmp_path)) / "config.yaml").read_text())
    assert config["blast"]["do"] is True
    assert config["sequana_taxonomy"]["store_unclassified"] is True


@pytest.mark.skipif(
    shutil.which("kraken2") is None or shutil.which("sequana_taxonomy") is None,
    reason="kraken2 / sequana_taxonomy not available; functional run skipped",
)
def test_functional(tmp_path):
    # End-to-end run on valid (non-empty) data. The deliberately-empty
    # `empty_R1_.fastq.gz` fixture is excluded as it has no reads to classify.
    runner = CliRunner()
    results = runner.invoke(
        main,
        [
            "--input-directory",
            sharedir,
            "--exclude-pattern",
            "empty",
            "--working-directory",
            str(tmp_path),
            "--force",
            "--databases",
            krakendb,
        ],
    )
    assert results.exit_code == 0

    completed = subprocess.run(["bash", "multitax.sh"], cwd=str(tmp_path))
    assert completed.returncode == 0

    assert (tmp_path / "data" / "summary.html").exists()
    assert (tmp_path / "data" / "kraken" / "kraken.csv").exists()
