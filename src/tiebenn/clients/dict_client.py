from datetime import datetime
from typing import Tuple

from tiebenn.clients.base_client import BaseClient
from tiebenn.config.input_params import InputParams
from tiebenn.constants.default_config_params import DEFAULT_CLIENT, DEFAULT_SDS_DIR, DEFAULT_PICKER
from tiebenn.constants.default_config_params import DEFAULT_EPIC_DIST, MIN_DETECTIONS, DEFAULT_SECS_BEFORE
from tiebenn.constants.default_config_params import DEFAULT_PH_ASSOC, DEFAULT_VEL_MODE, DEFAULT_VELMOD
from tiebenn.constants.default_config_params import DEFAULT_PLOTS, DEFAULT_DENOISE, DEFAULT_MULT_WINDOWS, DEFAULT_NLL3D
from tiebenn.logger.logger_settings import logger
from tiebenn.tools.algorithm import process_event


class DictClient(BaseClient):

    def __init__(self, input_msg_dict: dict, input_param_dict: dict = {}):
        self.input_msg_dict = input_msg_dict
        self.input_param_dict = input_param_dict

    @property
    def get_input_params(self) -> InputParams:
        return InputParams(event_file=None,
                           client=self.input_param_dict.get("client", DEFAULT_CLIENT),
                           sds_dir=self.input_param_dict.get("sds_dir", DEFAULT_SDS_DIR),
                           picker=self.input_param_dict.get("picker", DEFAULT_PICKER),
                           ph_assoc=self.input_param_dict.get("ph_assoc", DEFAULT_PH_ASSOC),
                           vel_mode=self.input_param_dict.get("vel_mode", DEFAULT_VEL_MODE),
                           velmod=self.input_param_dict.get("velmod", DEFAULT_VELMOD),
                           max_epic_dist=self.input_param_dict.get("max_epic_dist", DEFAULT_EPIC_DIST),
                           min_detections=self.input_param_dict.get("min_detections", MIN_DETECTIONS),
                           secs_before=self.input_param_dict.get("secs_before", DEFAULT_SECS_BEFORE),
                           plots=self.input_param_dict.get("plots", DEFAULT_PLOTS),
                           denoise=self.input_param_dict.get("denoise", DEFAULT_DENOISE),
                           mult_windows=self.input_param_dict.get("mult_windows", DEFAULT_MULT_WINDOWS),
                           nll3d=self.input_param_dict.get("nll3d", DEFAULT_NLL3D))

    def get_event_time_str(self) -> str:
        t_original: datetime = self.input_msg_dict.get("orig")
        return t_original.strftime("%Y-%m-%d %H:%M:%S")

    def get_event_coordinates(self) -> Tuple[float, float]:
        lat = self.input_msg_dict.get("lat")
        lon = self.input_msg_dict.get("lon")
        return lat, lon

    def run(self) -> dict:
        logger.info("DictClient is used for input data.")
        ev_time = self.get_event_time_str()
        lat, lon = self.get_event_coordinates()
        params = self.get_input_params

        results = process_event(params=params, ev_lat=lat, ev_lon=lon, ev_time=ev_time)
        return results
