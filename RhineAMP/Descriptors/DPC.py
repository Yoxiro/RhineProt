import pandas


def DPC(Protein_Sequence: str, header=None) -> pandas.Series:
    """

    :param header:
    :param Protein_Sequence: str,uppercase
    :return:
    """
    if header is None:
        header = "0"
    Protein_Sequence = Protein_Sequence.strip()
    AminoAcid_Length = len(Protein_Sequence)
    index_string = ["A", "I", "L", "V", "M", "F", "W", "Y",
                    "N", "C", "Q", "S", "T",
                    "R", "H", "K",
                    "D", "E",
                    "G", "P"]
    index_string_after = [i + j for i in index_string for j in index_string]
    aac_series = pandas.Series([0 for _ in index_string_after],
                               index=index_string_after,
                               dtype=float,
                               name=header)
    for char_index in range(len(Protein_Sequence) - 1):
        chars = Protein_Sequence[char_index:char_index + 2]
        if chars not in index_string_after:
            raise Exception("Invalid amino acid symbols")
        aac_series[chars] += 1
    aac_series = aac_series / (AminoAcid_Length - 1)
    return aac_series


def DPC_2D_array(Protein_Sequence: str) -> pandas.DataFrame:
    # Protein_Sequence = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTESA"
    index_string = ["A", "I", "L", "V", "M", "F", "W", "Y",
                    "N", "C", "Q", "S", "T",
                    "R", "H", "K",
                    "D", "E",
                    "G", "P"]
    index_string_after = [i + j for i in index_string for j in index_string]
    dpc_array = pandas.DataFrame([[0.0] * 20] * 20,
                                 index=index_string,
                                 columns=index_string,
                                 dtype=float)
    for char_index in range(len(Protein_Sequence) - 1):
        chars = Protein_Sequence[char_index:char_index + 2]
        if chars not in index_string_after:
            raise Exception("Invalid amino acid symbols")
        dpc_array.loc[Protein_Sequence[char_index], Protein_Sequence[char_index + 1]] += 1
    dpc_array = dpc_array / (len(Protein_Sequence) - 1)
    return dpc_array
