import argparse

from tiebenn.clients.cli_client import CliClient
from tiebenn.logger.logger_settings import logger
from tiebenn.tools.nicetools import str2bool


def none_or_str(value):
    if value.lower() == "none":
        return None
    return value


def get_params_from_terminal() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--event_file', type=str,
                        help='The complete path and filename containing the event/s for hypocenter estimation')
    parser.add_argument('--max_epic_dist', type=float,
                        help='The stations used for arrival-time picking will be withing this epicentral distance (in km)')
    parser.add_argument('--picker', default='PhaseNet', type=str,
                        help='Picker used to phase-picking: SeisBench (eqt or pn)')
    parser.add_argument('--nll3d', default=False, type=str2bool,
                        help='Use 3D velocity model for depth estimation in NonLinLoc')
    parser.add_argument('--client', default='FDSN', type=str, help='FDSN or SDS')
    parser.add_argument('--sds_dir', default='/', type=none_or_str, help='full path to SeisComp3 directory')
    parser.add_argument('--min_detections', default=3, type=int,
                        help='Minimal detections required to hypocenter estimation')
    parser.add_argument('--plots', default=False, type=str2bool, help='If True, it will plot several figures')
    parser.add_argument('--vel_mode', default='auto', type=str,
                        help='Velocity model mode: automatic or manual selection. If manual, velmod must be specified')
    parser.add_argument('--velmod', type=int,
                        help='Velocity model for phase association and hypocenter estimation with NonLinLoc')
    parser.add_argument('--denoise', default=False, type=str2bool,
                        help='True to use DeepDenoiser on collected waveforms up to a distance of 100 km')
    parser.add_argument('--ph_assoc', default='no', type=str,
                        help='Phase associator. Options are GaMMa and PyOcto (not case sensitive)')
    parser.add_argument('--mult_windows', default=False, type=str2bool,
                        help='Picks in windows with variable start time')
    parser.add_argument('--secs_before', default=0, type=int,
                        help='Seconds before origin time in case mult_window=False')

    return parser.parse_args()


def main_cli():
    params = get_params_from_terminal()
    client_terminal = CliClient(params)
    logger.info("CliClient class is used: input parameters are given from the terminal.")
    client_terminal.run()
