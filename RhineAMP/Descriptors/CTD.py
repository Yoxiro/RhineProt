import pandas
from typing import Dict

def _Classification_Validation(Classification:Dict[str,str]):
    sum_of_values = sum(
        list(
            map(
                lambda x:len(x),
                Classification.values()
            )
        )
    )
    if 20 != sum_of_values:
        raise Exception("Invalid Classification")
def CTD_C(Protein_Sequence:str,Classification:Dict[str,str] = None)->pandas.Series:
    """

    :param Protein_Sequence:
    :param Classification: the classification of proteins
    :return: the C part of CTD of the given protein sequence
    """
    Protein_Sequence = Protein_Sequence.strip()

    _Classification_Validation(Classification)
    values = list(Classification.values())

    t1 = list(values[0])
    t2 = list(values[1])
    t3 = list(values[2])
    char2class = dict()
    for char in t1:
        char2class[char] = 0
    for char in t2:
        char2class[char] = 1
    for char in t3:
        char2class[char] = 2

    group = {"CTD_TYPE_1":0,
             "CTD_TYPE_2":0,
             "CTD_TYPE_3":0}

    ctdc = pandas.Series(group,name="CTD_C",dtype=float)
    for index in range(len(Protein_Sequence)):
        char = Protein_Sequence[index]

        ctdc.iloc[char2class[char]] += 1.0
    ctdc = ctdc/len(Protein_Sequence)
    return ctdc

def _CTD_C(Protein_Sequence:str,Classification:Dict[str,str]):
    pass
if __name__ == "__main__":
    Protein_Sequence = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTESA"

    a = {
        '1': 'RKEDQN',
        "2": 'GASTPHY',
        "3": 'CLVIMFW'
    }
    x = CTD_C(Protein_Sequence,a)
    print(x)