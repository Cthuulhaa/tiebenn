
import unittest
import argparse

import pytest

from tiebenn.clients.cli_client import CliClient
from tiebenn.config.input_params import InputParams



def get_namespace_args():
    args = argparse.Namespace(
        client="fdsn",
        sds_dir=None,
        velmod=2,
        nll3d=False,
        secs_before=0,
        event_file="data/dummy_event.csv",
        max_epic_dist=50,
        ph_assoc="g",
        picker="SeisBench_PhaseNet",
        denoise=True,
        vel_mode="auto",
        plots=True,
        min_detections=3,
        mult_windows=False
    )

    return args


class TestCliClient(unittest.TestCase):

    def test_get_event_time_returns_value(self):
        x = "2024 12 24 01:55:04 51.337 12.548"
        y = "2024-12-24 01:55:04 51.337 12.548"
        expected_output = "2024-12-24 01:55:04"

        t = CliClient.get_event_time(x, y)

        self.assertIsInstance(t, str)
        self.assertEqual(t, expected_output)

    def test_get_event_time_raises_exception(self):
        x = ""
        y = ""
        expected_output = ('Accepted Datetime formats for the events are: "dd-Mon-yyyy hh:mm:ss",'
                           ' "yyyy-mm-ddThh:mm:ss", or "yyyy-mm-dd hh:mm:ss"')

        with pytest.raises(ValueError) as exc_info:
            CliClient.get_event_time(x, y)
            assert expected_output in str(exc_info.value)

    def test_get_lon_lat_returns_values(self):
        x = "2024 12 24 01:55:04 51.337 12.548"
        expected_output = (12.548, 51.337)

        ev_lon, ev_lat = CliClient.get_lon_lat(x)

        self.assertIsInstance(ev_lon, float)
        self.assertIsInstance(ev_lat, float)
        self.assertEqual(ev_lon, expected_output[0])
        self.assertEqual(ev_lat, expected_output[1])

    def test_check_client_raises_exception(self):
        client = "sds"
        sds_dir = None
        expected_output = 'Parameter sds_dir must be set for SDS client'

        with pytest.raises(ValueError) as e:
            CliClient.check_client(client, sds_dir)
            assert expected_output in str(e.value)

    def test_check_velocity_model_vel_model_manual_velmod_none_raises_exception(self):
        vel_mode = "m"
        velmod = None
        expected_output = 'Velocity model must be selected in manual mode'

        with pytest.raises(ValueError) as e:
            CliClient.check_velocity_model(vel_mode, velmod)
            assert expected_output in str(e.value)

    def test_check_velocity_model_raises_file_not_found(self):
        vel_mode = "m"
        velmod = 222
        expected_output = f"Velocity model v{velmod} does not exist"

        with pytest.raises(FileNotFoundError) as e:
            CliClient.check_velocity_model(vel_mode, velmod)
            assert expected_output in str(e.value)

    def test_get_params_returns_params(self):

        params = get_namespace_args()
        cc = CliClient(params)
        p = cc.get_input_params

        self.assertIsInstance(p, InputParams)


if __name__ == "__main__":
    unittest.main()
