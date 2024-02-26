import pandas


def AAC(AminoAcid_Sequence: str) -> pandas.Series:
    """
    calculate the AAC of the input AminoAcid_Sequence
    :param AminoAcid_Sequence: str,uppercase
    :return:
    """
    AminoAcid_Sequence = AminoAcid_Sequence.strip()
    AminoAcid_Length = len(AminoAcid_Sequence)

    index_string = ["A", "I", "L", "V", "M", "F", "W", "Y",
                    "N", "C", "Q", "S", "T",
                    "R", "H", "K",
                    "D", "E",
                    "G", "P"]
    aac_series = pandas.Series([0 for _ in index_string], index=index_string, dtype=float)
    for char in AminoAcid_Sequence:
        if char not in index_string:
            raise Exception("Invalid amino acid symbols")
        aac_series[char] += 1
    aac_series = aac_series / AminoAcid_Length
    return aac_series


def DPC(AminoAcid_Sequence: str) -> pandas.Series:
    """

    :param AminoAcid_Sequence: str,uppercase
    :return:
    """
    AminoAcid_Sequence = AminoAcid_Sequence.strip()
    AminoAcid_Length = len(AminoAcid_Sequence)
    index_string = ["A", "I", "L", "V", "M", "F", "W", "Y",
                    "N", "C", "Q", "S", "T",
                    "R", "H", "K",
                    "D", "E",
                    "G", "P"]
    index_string_after = [i + j for i in index_string for j in index_string]
    aac_series = pandas.Series([0 for _ in index_string_after], index=index_string_after, dtype=float)
    for char_index in range(len(AminoAcid_Sequence) - 1):
        chars = AminoAcid_Sequence[char_index:char_index+2]
        if chars not in index_string_after:
            raise Exception("Invalid amino acid symbols")
        aac_series[chars] += 1
    aac_series = aac_series / (AminoAcid_Length - 1)
    return aac_series


if __name__ == "__main__":
    LL_37 = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES"
    aac_of_LL37 = DPC(LL_37)
    print(aac_of_LL37)
    print(aac_of_LL37["LL"])
