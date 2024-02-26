# Core Library
import json
import math
from typing import Any, Dict, List

import pandas

# First Party Library
from RhineAMP.Descriptors.AAC import AAC

from pkg_resources import resource_filename

with open(resource_filename(__name__, "Data/hydrophobicity.json"), "r") as f:
    _Hydrophobicity: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "Data/hydrophilicity.json"), "r") as f:
    _hydrophilicity: Dict[str, float] = json.load(f)

with open(resource_filename(__name__, "Data/residuemass.json"), "r") as f:
    _residuemass: Dict[str, float] = json.load(f)

# with open(resource_filename(__name__, "Data/pK1.json"), "r") as f:
#     _pK1: Dict[str, float] = json.load(f)
#
# with open(resource_filename(__name__, "Data/pK2.json"), "r") as f:
#     _pK2: Dict[str, float] = json.load(f)
#
# with open(resource_filename(__name__, "Data/pI.json"), "r") as f:
#     _pI: Dict[str, float] = json.load(f)

# print(_Hydrophobicity)
def _normalized_aac(Protein_Sequence:str)-> pandas.Series:
    aac_data = AAC(Protein_Sequence)
    mean = round(aac_data.mean(),8)
    std = round(aac_data.std(ddof=0),8)
    # For the situation that the frequencies of amino acids equals
    if std == 0:
        return pandas.Series([1]*aac_data.size,
                             index=aac_data.index,
                             dtype="float")
    else:
        return (aac_data-mean)/std