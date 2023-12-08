from pathlib import Path

import pytest

from tests.end_to_end_tests.run_directory_manipulations import run_trial_from_run_directory_template


def test_basic_oseek_step1():
    executable_path_from_root = Path('minuit_torbtpar_rvgMagsclast1/minuit_torbtpar_rvgMagsclast1.xO')
    executable_path = Path(__file__).parent.parent.parent.parent.parent.joinpath(executable_path_from_root)
    run_trial_from_run_directory_template(__file__, executable_path=executable_path)
