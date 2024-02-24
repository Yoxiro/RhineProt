import Descriptor

class AminoAcidDescriptor:
    def __init__(self,AminoAcid_Sequence:str):
        self.AminoAcid_Sequence = AminoAcid_Sequence.strip()
        self.AminoAcid_Length = len(AminoAcid_Sequence)
        self.AAC = Descriptor.AAC(self.AminoAcid_Sequence)
