import pandas
from typing import Dict

_combined_dict = {
    'hydrophobicity_PRAM900101': {
        "group1": "RKEDQN",
        "group2": "GASTPHY",
        "group3": "CLVIMFW"
    },
    'hydrophobicity_ARGP820101': {
        "group1": "QSTNGDE",
        "group2": "RAHCKMV",
        "group3": "LYPFIW"
    },
    'hydrophobicity_ZIMJ680101': {
        "group1": "QNGSWTDERA",
        "group2": "HMCKV",
        "group3": "LPFYI"
    },
    'hydrophobicity_PONP930101': {
        "group1": "KPDESNQT",
        "group2": "GRHA",
        "group3": "YMFWLCVI"
    },
    'hydrophobicity_CASG920101': {
        "group1": "KDEQPSRNTG",
        "group2": "AHYMLV",
        "group3": "FIWC"
    },
    'hydrophobicity_ENGD860101': {
        "group1": "RDKENQHYP",
        "group2": "SGTAW",
        "group3": "CVLIMF"
    },
    'hydrophobicity_FASG890101': {
        "group1": "KERSQD",
        "group2": "NTPG",
        "group3": "AYHWVMFLIC"
    }
}
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
    if Classification is None:
        Classification = _combined_dict["hydrophobicity_PRAM900101"]
    else:
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

    group = {"CTD_TYPE_1": 0,
             "CTD_TYPE_2": 0,
             "CTD_TYPE_3": 0}

    ctdc = pandas.Series(group, name="CTD_C", dtype=float)
    for index in range(len(Protein_Sequence)):
        char = Protein_Sequence[index]

        ctdc.iloc[char2class[char]] += 1.0
    ctdc = ctdc / len(Protein_Sequence)
    return ctdc
def CTD_T(Protein_Sequence:str,Classification:Dict[str,str]=None):

    Protein_Sequence = Protein_Sequence.strip()
    if Classification is None:
        Classification = _combined_dict["hydrophobicity_PRAM900101"]
    else:
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

    ctdt_temp = pandas.DataFrame([[0]*3]*3, dtype=float)

    for index in range(len(Protein_Sequence)-1):
        char1 = Protein_Sequence[index]
        char2 = Protein_Sequence[index+1]
        # If the former char is the same as thr latter one, then there is no transition
        if char2 != char1:
            ctdt_temp.iloc[char2class[char1],char2class[char2]] += 1

    FirstandSecond = ctdt_temp.iloc[0,1] +ctdt_temp.iloc[1,0]
    FirstandThird = ctdt_temp.iloc[0,2] +ctdt_temp.iloc[2,0]
    SecondandThird = ctdt_temp.iloc[1,2] +ctdt_temp.iloc[2,1]

    ctdt = pandas.Series([FirstandSecond,
                          FirstandThird,
                          SecondandThird],
                         dtype=float,
                         index=["Type1_2","Type1_3","Type2_3"],
                         name = "CTD_T")
    ctdt = ctdt / (len(Protein_Sequence)-1)
    return ctdt

def _CTD_T(Protein_Sequence:str,Classification:Dict[str,str]):

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

    ctdt_temp = pandas.DataFrame([[0] * 3] * 3, dtype=float)

    for index in range(len(Protein_Sequence) - 1):
        char1 = Protein_Sequence[index]
        char2 = Protein_Sequence[index + 1]
        # If the former char is the same as thr latter one, then there is no transition
        if char2 != char1:
            ctdt_temp.iloc[char2class[char1], char2class[char2]] += 1

    FirstandSecond = ctdt_temp.iloc[0, 1] + ctdt_temp.iloc[1, 0]
    FirstandThird = ctdt_temp.iloc[0, 2] + ctdt_temp.iloc[2, 0]
    SecondandThird = ctdt_temp.iloc[1, 2] + ctdt_temp.iloc[2, 1]

    ctdt = pandas.Series([FirstandSecond,
                          FirstandThird,
                          SecondandThird],
                         dtype=float,
                         index=["Type1_2", "Type1_3", "Type2_3"],
                         name="CTD_T")
    ctdt = ctdt / (len(Protein_Sequence) - 1)
    return ctdt

if __name__ == "__main__":
    Protein_Sequence = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTESA"

    a = {
        '1': 'RKEDQN',
        "2": 'GASTPHY',
        "3": 'CLVIMFW'
    }
    x = CTD_T(Protein_Sequence,a)
    print(x)