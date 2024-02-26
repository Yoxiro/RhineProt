# Core Library
import json
import math
from typing import Any, Dict, List

from pkg_resources import resource_filename

with open(resource_filename(__name__, "data/hydrophobicity.json"), "r") as f:
    _Hydrophobicity: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "data/hydrophilicity.json"), "r") as f:
    _hydrophilicity: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "data/residuemass.json"), "r") as f:
    _residuemass: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "data/pK1.json"), "r") as f:
    _pK1: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "data/pK2.json"), "r") as f:
    _pK2: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "data/pI.json"), "r") as f:
    _pI: Dict[str, float] = json.load(f)




