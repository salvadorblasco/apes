from dataclasses import dataclass
import datetime

import consts

@dataclass
class Project:
    default_dir: str = '.'
    author: str = 'Anonymous'
    title: str = 'No title'
    comments: str = ''
    created: str = ''
    last_modified: str = ''
    temperature: float = 298.15
    filename: str = ''
    weighting: int = consts.WEIGHT_AUTO

    def __post_init__(self):
        self.created = str(datetime.datetime.now())
        self.last_modified = str(datetime.datetime.now())
        self.set_weight_auto()

    def set_weight_auto(self):
        self.weighting = consts.WEIGHT_AUTO

    def set_weight_unit(self):
        self.weighting = consts.WEIGHT_UNIT
