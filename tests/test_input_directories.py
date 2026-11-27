from unittest.mock import patch

from exotic.inputs import imaging_files, save_directory


def test_imaging_files_reports_invalid_directory_before_retry(tmp_path):
    invalid_dir = tmp_path / 'misspelled imaging'
    image = tmp_path / 'observation.fits'
    image.touch()

    with patch('exotic.inputs.log_info') as log_info, \
            patch('exotic.inputs.user_input', return_value=str(tmp_path)) as user_input:
        assert imaging_files(str(invalid_dir)) == [str(image)]

    assert repr(str(invalid_dir)) in log_info.call_args.args[0]
    assert log_info.call_args.kwargs['error'] is True
    user_input.assert_called_once()


def test_save_directory_reports_invalid_directory_before_retry(tmp_path):
    invalid_dir = tmp_path / 'misspelled results'

    with patch('exotic.inputs.log_info') as log_info, \
            patch('exotic.inputs.user_input', return_value=str(tmp_path)) as user_input:
        assert save_directory(str(invalid_dir)) == str(tmp_path)

    assert repr(str(invalid_dir)) in log_info.call_args.args[0]
    assert log_info.call_args.kwargs['error'] is True
    user_input.assert_called_once()
