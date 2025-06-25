from dataclasses import dataclass, field
import datetime

import consts

@dataclass
class Project:
    """Manage options for project."""

    default_dir: str = '.'
    author: str = 'Anonymous'
    title: str = 'No title'
    comments: str = ''
    created: str = field(default_factory=lambda: str(datetime.datetime.now()))
    last_modified: str = field(default_factory=lambda: str(datetime.datetime.now()))
    temperature: float = 298.15     # in kelvin
    filename: str = ''
    weighting: int = consts.WEIGHT_AUTO

    def __post_init__(self):
        if not self.created:
            self.created = str(datetime.datetime.now())
        if not self.last_modified:
            self.last_modified = str(datetime.datetime.now())
        self.set_weight_auto()

    def set_weight_auto(self) -> None:
        """Set weighting to automatic mode."""
        self.weighting = consts.WEIGHT_AUTO

    def set_weight_unit(self) -> None:
        """Set weighting to unit mode."""
        self.weighting = consts.WEIGHT_UNIT
