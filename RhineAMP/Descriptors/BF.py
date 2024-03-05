import pandas

_AAletter = list("ACDEFGHIKLMNPQRSTVWY")


def BF(Protein_Sequence: str) -> pandas.Series:
    A2Index = {_AAletter[index]: index for index in range(20)}
    Serieses = [_generateBF(A2Index[char]) for char in Protein_Sequence]
    res = pandas.concat(Serieses)
    return res


def _generateBF(index: int):
    res = pandas.Series([0] * 20)
    res[index] = 1
    return res


if __name__ == "__main__":
    LL37 = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTESA"
    BF(LL37)
