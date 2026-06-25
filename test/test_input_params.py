import unittest
from typing import Literal

import pytest
from pydantic import ValidationError

from tiebenn.config.input_params import InputParams


class TestInputParams(unittest.TestCase):

    def test_input_params_default_values(self):
        config = InputParams()

        assert config.client == Literal["fdsn"]
        assert config.sds_dir == None
        assert config.picker == Literal["seisbench_phasenet"]
        assert config.ph_assoc == Literal["pyocto"]
        assert config.vel_mode == Literal["auto"]
        assert config.velmod is None
        assert config.max_epic_dist == 150.0
        assert config.min_detections == 3
        assert config.secs_before == 0
        assert config.plots is False
        assert config.denoise is True
        assert config.mult_windows is False
        assert config.nll3d is False

    # set values and values types
    def test_input_params_min_detection_is_smaller_than_3(self):
        config = InputParams(min_detections=1)
        assert config.min_detections == 3
        assert isinstance(config.min_detections, int)

    def test_input_params_min_detection_is_larger_than_3(self):
        config = InputParams(min_detections=10)
        assert config.min_detections == 10
        assert isinstance(config.min_detections, int)

    def test_input_params_nll3d_is_set_to_true_returns_false(self):
        config = InputParams(nll3d=True)
        assert config.nll3d is False
        assert isinstance(config.nll3d, bool)

    def test_input_params_client_returns_lowercase_value(self):
        config = InputParams(client="FDSN")
        assert config.client == "fdsn"
        assert isinstance(config.client, str)

    def test_input_params_ph_assoc_returns_lowercase_value(self):
        config = InputParams(ph_assoc="P")
        assert config.ph_assoc == "p"
        assert isinstance(config.ph_assoc, str)

    def test_input_params_picker_returns_lowercase_value(self):
        config = InputParams(picker="PN")
        assert config.picker == "pn"
        assert isinstance(config.picker, str)

    def test_input_params_secs_before_returns_int(self):
        config = InputParams(secs_before=1.0)
        assert config.secs_before == 1
        assert isinstance(config.secs_before, int)

    def test_input_params_sds_dir_value_str(self):
        config = InputParams(sds_dir="sds_dir")
        assert isinstance(config.sds_dir, str)

    def test_input_params_sds_dir_value_none(self):
        config = InputParams(sds_dir=None)
        assert config.sds_dir is None

    def test_input_params_denoise_value_bool(self):
        config = InputParams()
        assert isinstance(config.denoise, bool)

    def test_input_params_plots_value_bool(self):
        config = InputParams()
        assert config.plots is False
        assert isinstance(config.plots, bool)

    def test_input_params_max_epic_dist_value_bool(self):
        config = InputParams(max_epic_dist=50)
        assert isinstance(config.max_epic_dist, float)
        assert config.max_epic_dist == pytest.approx(50)

    # invalid values
    def test_input_params_client_invalid(self):
        with pytest.raises(ValidationError):
            InputParams(client="test")

    def test_input_params_picker_invalid(self):
        with pytest.raises(ValidationError):
            InputParams(picker="test")

    def test_input_params_ph_assoc_invalid(self):
        with pytest.raises(ValidationError):
            InputParams(ph_assoc="test")

    def test_input_params_vel_mode_invalid(self):
        with pytest.raises(ValidationError):
            InputParams(vel_mode="test")

    # valid values
    def test_input_params_valid_client(self):
        for c in ["sds", "fdsn"]:
            p = InputParams(client=c)
            assert p.client == c

    def test_input_params_valid_picker(self):
        for c in ['eqtransformer', 'eqt', 'sb_eqt', 'seisbench_eqt', 'seisbench_eqtransformer', 'sb_eqtransformer',
                  'phasenet', 'pn', 'sb_pn', 'seisbench_pn', 'sb_phasenet', 'seisbench_phasenet']:
            p = InputParams(picker=c)
            assert p.picker == c

    def test_input_params_valid_vel_mode(self):
        for c in ['manual', 'man', 'm', 'automatic', 'auto', 'a']:
            p = InputParams(vel_mode=c)
            assert p.vel_mode == c

    def test_input_params_valid_ph_assoc(self):
        for c in ['gamma', 'g', 'pyocto', 'p']:
            p = InputParams(ph_assoc=c)
            assert p.ph_assoc == c

    def test_input_params_valid_velmod_int(self):
        for c in [0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 16, 18, 19, 20]:
            p = InputParams(velmod=c)
            assert p.velmod == c
            assert isinstance(p.velmod, int)

    def test_input_params_valid_velmod_None(self):
        p = InputParams(velmod=None)
        assert p.velmod is None

    def test_input_params_valid_event_file_none(self):
        p = InputParams(event_file=None)
        assert p.event_file is None

    def test_input_params_valid_event_file_str(self):
        p = InputParams(event_file="dummy_file")
        assert isinstance(p.event_file, str)


if __name__ == "__main__":
    unittest.main()
