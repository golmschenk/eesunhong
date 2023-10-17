from pathlib import Path

import pytest

from tests.end_to_end_tests.run_directory_manipulations import run_trial_with_output_file_cleaning


@pytest.mark.slow
def test_advanced_dseek():
    executable_path_from_root = Path('reference/minuit_torbtpar_rvgMagsclast1/minuit_torbtpar_rvgMagsclast1.xO')
    executable_path = Path(__file__).parent.parent.parent.parent.parent.joinpath(executable_path_from_root)
    run_trial_with_output_file_cleaning(__file__, executable_path=executable_path)
