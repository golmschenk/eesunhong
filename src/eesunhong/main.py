import subprocess
import sys
from pathlib import Path
import platform
import os


def main():
    module_parent_directory = Path(__file__).parent.parent
    executable_directory = module_parent_directory
    if executable_directory.name == 'src' and executable_directory.joinpath('../build').exists():
        executable_directory = executable_directory.joinpath('../build')
    executable_path = executable_directory.joinpath('eesunhong_main')
    platform_system_name = platform.system()
    subprocess_run_kwargs = {}
    if platform_system_name == 'Windows':
        # On Windows, the path to the DLLs packaged with the wheel needs to be explicitly added.
        environment = os.environ.copy()
        environment["PATH"] = f'{Path(__file__).parent.parent.joinpath("eesunhong.libs")};{environment["PATH"]}'
        subprocess_run_kwargs['env'] = environment
    subprocess.run([executable_path, *sys.argv], **subprocess_run_kwargs)
