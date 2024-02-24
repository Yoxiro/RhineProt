import pandas

def AAC(AminoAcid_Sequence:str) -> pandas.Series:
    """
    calculate the AAC of the input AminoAcid_Sequence
    :param AminoAcid_Sequence: str,uppercase
    :return:
    """
    AminoAcid_Sequence = AminoAcid_Sequence.strip()
    AminoAcid_Length = len(AminoAcid_Sequence)

    index_string = ["A","I","L","V","M","F","W","Y",
                    "N","C","Q","S","T",
                    "R","H","K",
                    "D","E",
                    "G","P"]
    aac_dict = {item:0 for item in index_string}
    for char in AminoAcid_Sequence:
        if char not in index_string:
            raise Exception("Invalid amino acid symbols")
        aac_dict[char] +=1
    aac_series = pandas.Series(aac_dict,dtype=float)
    aac_series = aac_series/AminoAcid_Length
    return aac_series

if __name__ == "__main__":
    LL_37 = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES"
    aac_of_LL37=AAC(LL_37)
    print(aac_of_LL37)