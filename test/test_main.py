import easydev
import os
import tempfile
import subprocess
import sys
#from sequana.pipelines_common import get_pipeline_location as getpath
#sharedir = getpath('multitax')
from . import test_dir

sharedir = f"{test_dir}/data" 

from sequana.scripts import taxonomy
try:
    taxonomy.main(["tst", '--download', 'toydb'])
except:
    pass



def test_standalone_subprocess():
    directory = tempfile.TemporaryDirectory()
    cmd = """sequana_pipelines_multitax --input-directory {}  
            --working-directory {} --force --databases toydb """.format(sharedir, directory.name)
    subprocess.call(cmd.split())


def test_standalone_script():
    directory = tempfile.TemporaryDirectory()
    import sequana_pipelines.multitax.main as m
    from sequana import sequana_config_path
    sys.argv = ["test", "--input-directory", sharedir, 
            "--working-directory", directory.name, "--force", "--databases",
            sequana_config_path + "/kraken_toydb"]
    m.main()



def test_version():
    cmd = "sequana_pipelines_multitax --version"
    subprocess.call(cmd.split())


