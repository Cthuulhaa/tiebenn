import unittest

from tiebenn.clients.dict_client import DictClient
from tiebenn.config.input_params import InputParams
from datetime import datetime


class TestDictClient(unittest.TestCase):

    """
    def test_get_input_returns_default_values_for_non_defined_keys(self):
        input_dict = {"client": "sds"}
        input_msg_dict = {}

        dc = DictClient(input_msg_dict=input_msg_dict, input_dict=input_dict)
        output = dc.get_input_params

        self.assertIsInstance(output, InputParams)
        self.assertEqual(len(output.__dict__), 14)
    """

    def test_get_event_time_str_returns_str(self):
        input_msg_dict = {"orig": datetime(2026, 3, 2, 13, 18, 10)}
        expected_output = "2026-03-02 13:18:10"

        dc = DictClient(input_msg_dict)
        output = dc.get_event_time_str()

        self.assertIsInstance(output, str)
        self.assertEqual(output, expected_output)
