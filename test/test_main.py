import easydev
import os
import tempfile
import subprocess
import sys
from . import test_dir

sharedir = f"{test_dir}/data"
krakendb = f"{test_dir}/data/krakendb"


def test_standalone_subprocess():
    directory = tempfile.TemporaryDirectory()
    cmd = f"""sequana_pipelines_multitax --input-directory {sharedir}
            --working-directory {directory.name} --force --databases {krakendb} """
    subprocess.call(cmd.split())


def test_standalone_script():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.multitax.main as m
    from sequana import sequana_config_path
    sys.argv = ["test", "--input-directory", sharedir, "--working-directory", directory.name, "--force", "--databases", krakendb]
    m.main()



def test_version():
    cmd = "sequana_pipelines_multitax --version"
    subprocess.call(cmd.split())


