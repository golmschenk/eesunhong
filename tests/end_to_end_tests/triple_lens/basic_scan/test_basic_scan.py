from pathlib import Path

from tests.end_to_end_tests.run_directory_manipulations import run_trial_with_output_file_cleaning
from tests.end_to_end_tests.triple_lens.executable_path import triple_lens_executable_path


def test_basic_scan():
    executable_path = Path(__file__).parent.parent.parent.parent.parent.joinpath(triple_lens_executable_path)
    run_trial_with_output_file_cleaning(__file__, executable_path=executable_path)
