import pandas
def AAC(Protein_Sequence: str,header = None) -> pandas.Series:
    """
    calculate the AAC of the input AminoAcid_Sequence
    :param Protein_Sequence: str,uppercase
    :return:
    """
    if header is None:
        header = "0"
    Protein_Sequence = Protein_Sequence.strip()
    Protein_Length = len(Protein_Sequence)

    index_string = ["A", "I", "L", "V", "M", "F", "W", "Y",
                    "N", "C", "Q", "S", "T",
                    "R", "H", "K",
                    "D", "E",
                    "G", "P"]
    aac_series = pandas.Series([0 for _ in index_string],
                               index=index_string,
                               dtype=float,
                               name=header)
    for char in Protein_Sequence:
        if char not in index_string:
            raise Exception("Invalid amino acid symbols")
        aac_series[char] += 1
    aac_series = aac_series / Protein_Length
    return aac_series