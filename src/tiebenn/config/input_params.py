from typing import Literal

from pydantic import BaseModel, field_validator, Field
from pydantic_core.core_schema import ValidationInfo

from tiebenn.constants.default_config_params import (DEFAULT_CLIENT, DEFAULT_SDS_DIR, DEFAULT_PICKER, DEFAULT_PH_ASSOC,
                                                     DEFAULT_VEL_MODE, DEFAULT_VELMOD, DEFAULT_EPIC_DIST,
                                                     DEFAULT_SECS_BEFORE, DEFAULT_EVENT_FILE,
                                                     DEFAULT_PLOTS, DEFAULT_DENOISE, DEFAULT_MULT_WINDOWS,
                                                     DEFAULT_NLL3D, MIN_DETECTIONS)
from tiebenn.logger.logger_settings import logger


class InputParams(BaseModel):
    event_file: str | None = DEFAULT_EVENT_FILE

    client: Literal['sds', 'fdsn'] = DEFAULT_CLIENT
    sds_dir: None | str = DEFAULT_SDS_DIR

    picker: Literal[
        'eqtransformer', 'eqt', 'sb_eqt', 'seisbench_eqt', 'seisbench_eqtransformer', 'sb_eqtransformer', 'phasenet', 'pn', 'sb_pn', 'seisbench_pn', 'sb_phasenet', 'seisbench_phasenet'] = DEFAULT_PICKER
    ph_assoc: Literal['gamma', 'g', 'pyocto', 'p'] = DEFAULT_PH_ASSOC

    vel_mode: Literal['manual', 'man', 'm', 'automatic', 'auto', 'a'] = DEFAULT_VEL_MODE
    velmod: Literal[0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 16, 18, 19, 20, None] = DEFAULT_VELMOD

    max_epic_dist: int | float = DEFAULT_EPIC_DIST
    min_detections: int = Field(default=MIN_DETECTIONS)
    secs_before: int | float | list[int] | list[float] = DEFAULT_SECS_BEFORE

    plots: bool = DEFAULT_PLOTS
    denoise: bool = DEFAULT_DENOISE
    mult_windows: bool = DEFAULT_MULT_WINDOWS
    nll3d: bool = DEFAULT_NLL3D

    @field_validator("picker", mode="before")
    @classmethod
    def normalize_picker(cls, value: str) -> str:
        return value.lower()

    @field_validator("client", mode="before")
    @classmethod
    def normalize_client(cls, value: str) -> str:
        return value.lower()

    @field_validator("ph_assoc", mode="before")
    @classmethod
    def normalize_ph_assoc(cls, value: str) -> str:
        return value.lower()

    @field_validator("nll3d", mode="before")
    @classmethod
    def check_nll3d_valid(cls, value: bool) -> bool:
        if value:
            logger.warning('3D grid depth estimation is temporarily disabled. Set nll3d=False')
            return False
        return value

    @field_validator("min_detections", mode="before")
    @classmethod
    def check_min_detections_valid(cls, value: int) -> int:
        if value < MIN_DETECTIONS:
            logger.warning(f'Minimum detections cannot be lower than {MIN_DETECTIONS}.')
            logger.warning(f'min_detections is set to {MIN_DETECTIONS}.')
            return MIN_DETECTIONS
        return value

    @field_validator("max_epic_dist", mode="before")
    @classmethod
    def check_max_epic_dist_valid(cls, value: float | int) -> float:
        return float(value)

    @field_validator("secs_before", mode="before")
    @classmethod
    def check_secs_before_valid(cls, value: int | float) -> int:
        return int(value)

    @field_validator("sds_dir", mode="before")
    @classmethod
    def set_sds_dir(cls, value: str | None, info: ValidationInfo) -> str | None:
        client = info.data.get("client")

        if client == "fdsn":
            return None
        if (client == "sds") and (value is None):
            return "/"

        return value
