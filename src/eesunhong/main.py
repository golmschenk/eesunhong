import subprocess
import sys
from pathlib import Path


def main():
    module_parent_directory = Path(__file__).parent.parent
    executable_directory = module_parent_directory
    if executable_directory.name == 'src' and executable_directory.joinpath('../build').exists():
        executable_directory = executable_directory.joinpath('../build')
    executable_path = executable_directory.joinpath('eesunhong_main')
    subprocess.run([executable_path, *sys.argv])
