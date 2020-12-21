from .version import __version__
# from .config import settings

import os
import yaml
config = yaml.safe_load(open(os.path.dirname(__file__) + "/config.yaml"))

