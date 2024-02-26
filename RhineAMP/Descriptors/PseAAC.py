# Core Library
import json
import math
from typing import Any, Dict, List

import pandas

# First Party Library
from RhineAMP.Descriptors.AAC import AAC

from pkg_resources import resource_filename


def _normalize_property(Properties: pandas.Series) -> pandas.Series:
    """
    normalize the given properties of amino acids
    :param Properties:
    :return:
    """
    mean = round(Properties.mean(), 8)
    std = round(Properties.std(ddof=0), 8)
    return (Properties - mean) / std


with open(resource_filename(__name__, "Data/hydrophobicity.json"), "r") as f:
    _hydrophobicity = _normalize_property(pandas.Series(json.load(f)))

with open(resource_filename(__name__, "Data/hydrophilicity.json"), "r") as f:
    _hydrophilicity = _normalize_property(pandas.Series(json.load(f)))

with open(resource_filename(__name__, "Data/residuemass.json"), "r") as f:
    _residuemass = _normalize_property(pandas.Series(json.load(f)))


# with open(resource_filename(__name__, "Data/pK1.json"), "r") as f:
#     _pK1: Dict[str, float] = json.load(f)
#
# with open(resource_filename(__name__, "Data/pK2.json"), "r") as f:
#     _pK2: Dict[str, float] = json.load(f)
#
# with open(resource_filename(__name__, "Data/pI.json"), "r") as f:
#     _pI: Dict[str, float] = json.load(f)

# print(_hydrophobicity)
# print(_hydrophobicity.mean())
# print(_hydrophobicity.std(ddof=0))

def _normalized_aac(Protein_Sequence: str) -> pandas.Series:
    aac_data = AAC(Protein_Sequence)
    mean = round(aac_data.mean(), 8)
    std = round(aac_data.std(ddof=0), 8)
    # For the situation that the frequencies of amino acids equals
    if std == 0:
        return pandas.Series([1] * aac_data.size,
                             index=aac_data.index,
                             dtype="float")
    else:
        return (aac_data - mean) / std


def PseAAC(Protein_Sequence: str,
           Lambda: int,
           weight: float = 0.05) -> pandas.Series:
    """
    calculate the PseAAC of the given Protein Sequence.
    The used amino acid properties in this function is hydrophobicity, hydrophilicity
    and residue mass
    :param Protein_Sequence:
        str,uppercase
    :param Lambda:
        integer
    :param weight:
        float
        the weight factor for the sequence order effect.
        Default value = 0.05
    :return:
        pandas.Series
    """

    # param validation
    Protein_Length = len(Protein_Sequence)
    if Protein_Length <= Lambda:
        raise Exception("Protein_Length should be larger than Lambda, Lambda larger instead")


def _get_correlation_function(Ri:str,Rj:str):
    """
    get the correlation function between Ri and Rj.
    :param Ri:
    :param Rj:
    :return:
    """
    return (math.pow(_hydrophobicity[Ri]-_hydrophobicity[Rj],2)+
            math.pow(_hydrophilicity[Ri]-_hydrophilicity[Rj],2)+
            math.pow(_residuemass[Ri]-_residuemass[Rj],2))/3
def _get_correlation_factor(Protein_Sequence:str,tier:int):
    Protein_Length = len(Protein_Sequence)
    cor_factor = 0
    for i in range(0,len(Protein_Sequence)-tier):
        """
        ABCDEFG LEN=7 TIER=1 
        then i can be 7-1-1
        """
        j = i + tier
        Ri = Protein_Sequence[i]
        Rj = Protein_Sequence[j]
        cor_factor += _get_correlation_function(Ri,Rj)
    return cor_factor/float((Protein_Length-tier))


