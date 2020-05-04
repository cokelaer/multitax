import easydev
import os
import tempfile
import subprocess
import sys
from sequana.pipelines_common import get_pipeline_location as getpath

sharedir = getpath('multitax')


def test_standalone_subprocess():
    directory = tempfile.TemporaryDirectory()
    cmd = """sequana_pipelines_multitax --input-directory {} 
            --working-directory {} --force""".format(sharedir, directory.name)
    subprocess.call(cmd.split())


def test_standalone_script():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.multitax.main as m
    sys.argv = ["test", "--input-directory", sharedir, 
            "--working-directory", directory.name, "--force"]
    m.main()

def test_full():

    with tempfile.TemporaryDirectory() as directory:
        wk = directory
        print(sharedir)
        cmd = 'sequana_pipelines_multitax --input-directory {} --input-pattern E*fastq.gz'
        cmd += ' --working-directory {}  --force --databases toydb '
        cmd = cmd.format(sharedir , wk)
        subprocess.call(cmd.split())

        stat = subprocess.call("sh multitax.sh".split(), cwd=wk)
        print(os.listdir(wk))

        assert os.path.exists(wk + "/multiqc/multiqc_report.html")


def test_version():
    cmd = "sequana_pipelines_multitax --version"
    subprocess.call(cmd.split())

