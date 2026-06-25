import argparse
import os
from datetime import datetime

import pytest
from tiebenn.clients.cli_client import CliClient
from tiebenn.clients.dict_client import DictClient


def get_input_message():
    return {"orig": datetime(2024, 12, 24, 1, 55, 4),
            "lat": 51.337,
            "lon": 12.548,
            "qvalue": 0.424,
            "locname": "Dummy Location",
            "detid": 1234567,
            "provider": "dummy_provider",
            "clusterlen": 1,
            "changecnt": 0,
            "evid": "qwert",
            "timelag": 512}


def get_input_params():
    return {"max_epic_dist": 150,
            "picker": "seisbench_phasenet",
            "client": "fdsn",
            "min_detections": 3,
            "plots": False,
            "vel_mode": "auto",
            "ph_assoc": "pyocto",
            "denoise": True,
            "mult_windows": False}


def get_namespace_args():
    # pf = Path(__file__).parent / "data" / "dummy_event.csv"
    args = argparse.Namespace(
        event_file=os.path.join(os.getcwd(), "data/dummy_event.csv"),  # pf.name,  # "data/dummy_event.csv",
        client="fdsn",
        sds_dir=None,
        ph_assoc="p",
        picker='seisbench_phasenet',
        denoise=True,
        vel_mode="auto",
        velmod=None,
        max_epic_dist=150,
        min_detections=3,
        secs_before=0,
        plots=False,
        mult_windows=False,
        nll3d=False
    )
    return args


def set_env():
    current_dir = os.getcwd()
    os.environ["PATH"] += os.pathsep + current_dir


def test_integration_cli_client():
    params = get_namespace_args()
    cc = CliClient(args=params)
    cc.run()
    assert 1 == 1


def test_integration_dict_client():
    expected_output = {'datetime': '24-12-2024_01:55:4.530991',
                       'latitude': 51.328221,
                       'longitude': 12.54744,
                       'depth': 22.558594,
                       'RMS': 0.292245}

    params = get_input_params()
    msg = get_input_message()
    dc = DictClient(input_msg_dict=msg, input_param_dict=params)
    result = dc.run()

    assert type(result) == dict
    assert len(result) == 5
    assert float(result['depth']) == pytest.approx(expected_output['depth'], abs=1.0)
