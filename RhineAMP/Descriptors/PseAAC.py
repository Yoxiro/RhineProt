# Core Library
import json
import math

import pandas
from pkg_resources import resource_filename

# First Party Library
from RhineAMP.Descriptors.AAC import AAC


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
    hydrophobicity = _normalize_property(pandas.Series(json.load(f)))

with open(resource_filename(__name__, "Data/hydrophilicity.json"), "r") as f:
    hydrophilicity = _normalize_property(pandas.Series(json.load(f)))

with open(resource_filename(__name__, "Data/residuemass.json"), "r") as f:
    residuemass = _normalize_property(pandas.Series(json.load(f)))

with open(resource_filename(__name__, "Data/pK1.json"), "r") as f:
    pK1 = _normalize_property(pandas.Series(json.load(f)))

with open(resource_filename(__name__, "Data/pK2.json"), "r") as f:
    pK2 = _normalize_property(pandas.Series(json.load(f)))

with open(resource_filename(__name__, "Data/pI.json"), "r") as f:
    pI = _normalize_property(pandas.Series(json.load(f)))


# print(_hydrophobicity)
# print(_hydrophobicity.mean())
# print(_hydrophobicity.std(ddof=0))

def _normalized_aac(Protein_Sequence: str) -> pandas.Series:
    """
    Get the AAC of the given protein sequence
    The normalization method here is to be discussed
    :param Protein_Sequence:
    :return:
    """
    aac_data = AAC(Protein_Sequence)
    return aac_data


def PseAAC(Protein_Sequence: str,
           Lambda: int,
           weight: float = 0.05,
           properties=None) -> pandas.Series:
    """
    calculate the PseAAC of the given Protein Sequence.
    The used amino acid properties in this function is hydrophobicity, hydrophilicity
    and residue mass
    :param properties:
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

    PseAAC_series = pandas.Series([0] * (20 + Lambda), dtype=float,
                                  index=["PseAAC" + str(i + 1) for i in range(20 + Lambda)])
    aac = _normalized_aac(Protein_Sequence)

    sum_of_theta = 0
    for i in range(1, Lambda + 1):
        sum_of_theta += _get_correlation_factor(Protein_Sequence, i, properties)
    numerator = 1 + weight * sum_of_theta
    PseAAC_series[0:20] = aac / numerator * 100
    for i in range(20, Lambda + 20):
        PseAAC_series[i] = (weight * _get_correlation_factor(Protein_Sequence, i - 19, properties)) \
                           / numerator * 100
    return PseAAC_series


def PseAAC_Amphiphilic(Protein_Sequence: str,
                       Lambda: int = 30,
                       weight: float = 0.05,
                       properties=None) -> pandas.Series:
    """
    calculate the PseAAC of the given Protein Sequence.
    The used amino acid properties in this function is hydrophobicity, hydrophilicity
    and residue mass
    :param properties:
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
    if properties is None:
        properties = [0, 1]
    properties_num = len(properties)
    Protein_Length = len(Protein_Sequence)
    if Protein_Length <= Lambda:
        raise Exception("Protein_Length should be larger than Lambda, Lambda larger instead")

    PseAAC_series = pandas.Series([0] * (20 + Lambda * properties_num), dtype=float,
                                  index=["PseAAC" + str(i + 1) for i in range(20 + Lambda*properties_num)])
    aac = _normalized_aac(Protein_Sequence)

    tau = _get_correlation_factor_Amphiphilic(Protein_Sequence=Protein_Sequence,
                                              Lambda=Lambda,
                                              properties=properties)
    tau_sum = tau.sum()
    # print(tau_sum)
    # print(aac)
    numerator = 1 + weight * tau_sum
    # print(numerator)
    PseAAC_series[0:20] = aac/numerator*100
    PseAAC_series[20:20 + Lambda * properties_num] = weight*tau/numerator*100
    # print(tau)
    return PseAAC_series


def _get_correlation_function(Ri: str, Rj: str,
                              properties=None):
    """
    get the correlation function between Ri and Rj.
    :param Ri:
    :param Rj:
    :return:
    """
    if properties is None:
        properties = [0, 1, 2]
    property_num = len(properties)
    properties_all = [hydrophilicity, hydrophobicity, residuemass, pI, pK1, pK2]
    res = 0.0
    for i in properties:
        res += math.pow(properties_all[i][Ri] - properties_all[i][Rj], 2)
    return res / property_num


def _get_correlation_factor_Amphiphilic(Protein_Sequence: str,
                                                properties,
                                                Lambda)->pandas.Series:
    if properties is None:
        properties = [0, 1]
    properties_all = [hydrophilicity, hydrophobicity, residuemass, pI, pK1, pK2]
    res_pandas = pandas.Series([0]*Lambda*len(properties),dtype=float)
    Protein_Length = len(Protein_Sequence)
    for i in range(Lambda):
        for j in range(len(properties)):
            res_pandas[i*len(properties)+j] = _sum_of_H(Protein_Sequence,tier=i+1,property=properties_all[properties[j]])/(Protein_Length-i-1)
    return res_pandas
def _sum_of_H(Protein_Sequence:str,tier:int,property:pandas.Series)->float:
    """
    求sumH
    :param Protein_Sequence:
    :param tier:
    :param property:
    :return:
    """
    Protein_Length = len(Protein_Sequence)
    res = 0.0
    for i in range(0,Protein_Length-tier):
        """
        ABCDEFG TIER=1 LENGTH=7
        """
        j = i + tier
        Ri = Protein_Sequence[i]
        Rj = Protein_Sequence[j]
        res = res + property[Ri]*property[Rj]
    return res


def _get_correlation_factor(Protein_Sequence: str, tier: int, properties):
    Protein_Length = len(Protein_Sequence)
    cor_factor = 0
    for i in range(0, len(Protein_Sequence) - tier):
        """
        ABCDEFG LEN=7 TIER=1 
        then i can be 7-1-1
        """
        j = i + tier
        Ri = Protein_Sequence[i]
        Rj = Protein_Sequence[j]
        cor_factor += _get_correlation_function(Ri, Rj, properties)
    return cor_factor / float((Protein_Length - tier))
