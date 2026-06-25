import logging
import os
import sys
from pathlib import Path

# for the container the logging path is:
# "/app/logs/app.log"

def get_log_path():
    log_dir = Path(os.getenv("LOG_DIR", "./logs"))
    try:
        log_dir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        log_dir = Path("./logs")
        log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir / "app.log"

# Configure logger
logger = logging.getLogger("tiebenn")  # give it a name
logger.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s | %(name)s | %(filename)s:%(lineno)d | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

# Create console handler
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(formatter)
console_handler.setLevel(logging.INFO)

# file logger
log_path = get_log_path()

# Create file handler
file_handler = logging.FileHandler(log_path, mode="w")  # it will overwrite the file every time the code is started!
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)

if not logger.handlers:
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
