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