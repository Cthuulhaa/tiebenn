from importlib.resources import files
from typing import Tuple

from tiebenn.clients.base_client import BaseClient
from tiebenn.config.input_params import InputParams
from tiebenn.constants.valid_config_params import VALID_VEL_MODE_MANUAL, VALID_VEL_MODE_AUTO
from tiebenn.logger.logger_settings import logger
from tiebenn.tools.algorithm import process_event
from tiebenn.tools.nicetools import strmonth2num


class CliClient(BaseClient):

    def __init__(self, args):
        self.args = args

    @property
    def get_input_params(self) -> InputParams:

        return InputParams(event_file=self.args.event_file,
                           max_epic_dist=self.args.max_epic_dist,
                           picker=self.args.picker.lower(),
                           client=self.args.client.lower(),
                           sds_dir=self.args.sds_dir,
                           min_detections=self.args.min_detections,
                           vel_mode=self.args.vel_mode.lower(),
                           velmod=self.args.velmod,
                           ph_assoc=self.args.ph_assoc.lower(),
                           secs_before=self.args.secs_before,
                           plots=self.args.plots,
                           denoise=self.args.denoise,
                           mult_windows=self.args.mult_windows,
                           nll3d=self.args.nll3d
                           )

    @staticmethod
    def get_lon_lat(x: str) -> Tuple[float, float]:
        ev_lon = float(x.split()[5])
        ev_lat = float(x.split()[4])
        return ev_lon, ev_lat

    @staticmethod
    def get_event_time(x: str, y: str) -> str:
        try:
            ev_year = x.split()[2]
            ev_month = strmonth2num(x.split()[1])
            ev_day = x.split()[0]
            ev_time = f"{ev_year}-{ev_month}-{ev_day} {x.split()[3]}"
            return ev_time
        except:
            try:
                ev_time = f"{y.split()[0]} {y.split()[1]}"
                return ev_time
            except:
                raise ValueError(
                    'Accepted Datetime formats for the events are: "dd-Mon-yyyy hh:mm:ss", "yyyy-mm-ddThh:mm:ss", or "yyyy-mm-dd hh:mm:ss"')

    # todo: move to input params (???)
    @staticmethod
    def check_client(client: str, sds_dir: str):
        if client == 'sds' and not sds_dir:
            raise ValueError('Parameter sds_dir must be set for SDS client')
        elif client == 'sds':
            logger.info(f'SDS directory set to: {sds_dir}')

    # todo: move to input params (???)
    @staticmethod
    def check_velocity_model(vel_mode: str, velmod: int):
        if vel_mode in VALID_VEL_MODE_MANUAL:
            if velmod is None:
                raise ValueError('Velocity model must be selected in manual mode')

            if velmod not in [6, 7, 12, 13, 17]:
                model_path = files('tiebenn.data.velocity_models').joinpath(f"v{velmod}")
                if not model_path.is_file():
                    raise FileNotFoundError(f"Velocity model v{velmod} does not exist")
                logger.info(f"Velocity model {velmod} selected.")
            elif velmod in [6, 7]:
                logger.info('Crust1.0 model will be used for seismic location')
            else:
                raise RuntimeError('Depth estimation with 3D grids currently disabled')

    # todo: move to input params (???)
    def check_vel_mode_valid(self, params: InputParams):
        if params.vel_mode.lower() in VALID_VEL_MODE_MANUAL:
            if params.velmod is None:
                raise ValueError('Velocity model must be selected in manual mode')
            CliClient.check_velocity_model(params.vel_mode, params.velmod)
        elif params.vel_mode.lower() not in VALID_VEL_MODE_AUTO:
            raise ValueError('vel_mode must be manual (m) or automatic (a).')

    def run(self) -> None:
        params = self.get_input_params
        self.check_client(params.client, params.sds_dir)
        self.check_vel_mode_valid(params)

        f = open(params.event_file, 'r')
        for x in f:
            x = x.replace('T', ' ')
            y = x
            x = x.replace('-', ' ').replace('_', ' ')

            ev_time = CliClient.get_event_time(x, y)
            lon, lat = CliClient.get_lon_lat(x)

            _ = process_event(params=params, ev_lat=lat, ev_lon=lon, ev_time=ev_time)
        return None
