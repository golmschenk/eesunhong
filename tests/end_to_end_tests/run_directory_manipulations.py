import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


def clean_up_run(run_directory):
    run_directory.joinpath('fit.lc_run_1').unlink(missing_ok=True)
    run_directory.joinpath('resid.run_1').unlink(missing_ok=True)
    run_directory.joinpath('run_1.dat').unlink(missing_ok=True)
    run_directory.joinpath('run_1.out').unlink(missing_ok=True)
    run_directory.joinpath('run_2.in').unlink(missing_ok=True)


def verify_directories_match(run_directory, expected_resulting_run_directory):
    expected_path_list = list(expected_resulting_run_directory.glob('*'))
    run_path_list = list(run_directory.glob('*'))
    verify_file_lists_match(run_directory, run_path_list, expected_path_list)
    for expected_path in expected_path_list:
        run_path = run_directory.joinpath(expected_path.name)
        assert run_path.exists()
        verify_run_files_match(run_path, expected_path)


def verify_file_lists_match(run_directory, run_path_list, expected_path_list):
    if set([path.name for path in run_path_list]) != set([path.name for path in expected_path_list]):
        lines_to_display = 100
        with open(run_directory.joinpath('run_1.out')) as run_1_out_file:
            run_1_out_content = run_1_out_file.readlines()
        print('============================', file=sys.stderr)
        print('Ending content of run_1.out ', file=sys.stderr)
        print(''.join(run_1_out_content[-lines_to_display:]), file=sys.stderr)
        print('============================', file=sys.stderr)
    assert set([path.name for path in run_path_list]) == set([path.name for path in expected_path_list])


def verify_run_files_match(run_path, expected_run_path):
    with expected_run_path.open() as expected_file, run_path.open() as run_file:
        line_number = 1
        expected_list = re.split(r'(\s+)', expected_file.read())
        actual_list = re.split(r'(\s+)', run_file.read())
        for (expected_item, actual_item) in zip(expected_list, actual_list):
            expected_white_space_match = re.fullmatch(r'\s+', expected_item)
            if expected_white_space_match is not None:
                actual_white_space_match = re.fullmatch(r'\s+', actual_item)
                assert actual_white_space_match is not None, f'''
                    When comparing the expected {expected_run_path} and the actual {run_path}
                    on line {line_number}, expected a segment of white space and found {actual_item}.  
                '''
                if '\n' in expected_item:
                    line_number += 1
                continue
            try:
                expected_number = float(expected_item)
                actual_number = float(actual_item)
                relative_tolerance = 0.01
                assert actual_number == pytest.approx(expected_number, rel=relative_tolerance), f'''
                    When comparing the expected {expected_run_path} and the actual {run_path}
                    on line {line_number}, the number {expected_number} was expected and the actual was {actual_number}
                    which does not match to a relative tolerance of {relative_tolerance}.  
                '''
                continue
            except ValueError:
                assert actual_item == expected_item, f'''
                    When comparing the expected {expected_run_path} and the actual {run_path}
                    on line {line_number}, the string {expected_item} was expected and the actual was {actual_item}.  
                '''
                continue


def verify_directories_file_list_does_not_match(run_directory, expected_resulting_run_directory):
    expected_path_list = list(expected_resulting_run_directory.glob('*'))
    run_path_list = list(run_directory.glob('*'))
    assert set([path.name for path in run_path_list]) != set([path.name for path in expected_path_list])


def run_trial_with_output_file_cleaning(current_file_path):
    run_directory = Path(current_file_path).parent.joinpath('run_directory')
    expected_directory = Path(current_file_path).parent.joinpath('expected_resulting_run_directory')
    clean_up_run(run_directory)
    run_trial(current_file_path, run_directory, expected_directory)
    clean_up_run(run_directory)


def run_trial_from_run_directory_template(current_file_path):
    template_run_directory = Path(current_file_path).parent.joinpath('template_run_directory')
    run_directory = Path(current_file_path).parent.joinpath('run_directory')
    expected_directory = Path(current_file_path).parent.joinpath('expected_resulting_run_directory')
    if run_directory.exists():
        shutil.rmtree(run_directory)
    shutil.copytree(template_run_directory, run_directory)
    run_trial(current_file_path, run_directory, expected_directory)
    shutil.rmtree(run_directory)


def run_trial(current_file_path, run_directory, expected_directory):
    executable_name = 'eesunhong_main'
    if shutil.which(executable_name) is not None:
        executable = executable_name
    else:
        executable = Path(current_file_path).parent.parent.parent.parent.joinpath(f'build/{executable_name}')
    verify_directories_file_list_does_not_match(run_directory, expected_directory)
    run_in_path = run_directory.joinpath('run_1.in')
    run_out_path = run_directory.joinpath('run_1.out')
    with run_in_path.open() as input_file, run_out_path.open('w') as output_file:
        subprocess.run([executable], cwd=run_directory, stdin=input_file,
                       stdout=output_file)
    verify_directories_match(run_directory, expected_directory)
