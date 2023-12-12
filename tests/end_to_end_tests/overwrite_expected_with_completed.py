import shutil
from pathlib import Path


def overwrite_expected_with_completed():
    root_source_directory = Path('tests/end_to_end_tests')
    for directory in [path for path in root_source_directory.glob('**/') if path.is_dir()]:
        run_directory = directory.joinpath('run_directory')
        if not run_directory.exists():
            continue
        expected_directory = run_directory.parent.parent.joinpath(directory.name).joinpath('expected_resulting_run_directory')
        for file_path in run_directory.iterdir():
            shutil.copy(file_path, expected_directory.joinpath(file_path.name))
            print(f'Copied {file_path.name}.')


if __name__ == '__main__':
    overwrite_expected_with_completed()
