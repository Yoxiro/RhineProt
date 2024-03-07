import RhineAMP.FastaRead
import multiprocessing
import pandas


def _FEGS_buffer(args: list):
    header = args[0]
    sequence = args[1]
    return RhineAMP.FEGS(sequence, header)


def _AAC_Buffer(args: list):
    header = args[0]
    sequence = args[1]
    return RhineAMP.AAC(sequence, header)


class Sequences:

    def __init__(self):
        self._df_to_csv = None
        self._Protein_Sequences = None

    def LoadFromFasta(self, filepath):
        self._Protein_Sequences = RhineAMP.FastaRead(filepath)

    def AAC(self,
            core: int = 1,
            save_result: bool = True,
            filepath=None):
        """

        :param enable_multiprocessing: whether able Multiprocessing
        :param core: the core to be used if Multiprocessing is enabled
        :param save_result: whether save the result
        :param filepath: the filepath
        :return: pandas.DataFrame
        """
        pool = multiprocessing.Pool(processes=core)
        result = pool.map(_AAC_Buffer, self._Protein_Sequences.values.tolist())
        df = pandas.concat(result, axis=1)
        if save_result:
            self._to_csv(filepath, df)
        return df

    def FEGS(self, core: int = 1, save_result: bool = True, filepath=None):
        """
        To apply FEGS to the Protein Sequences input
        :param core: the core of CPU to be used
        :param save_result: whether to save the file or not
        :param filepath: the filepath
        :return: pandas.DataFrame
        """
        if core is None:
            core = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=core)
        result = pool.map(_FEGS_buffer, self._Protein_Sequences.values.tolist())
        df = pandas.concat(result, axis=1)
        if save_result:
            self._to_csv(filepath, df)
        return df

    def _to_csv(self, filepath="Saved_File.csv", df=None):
        if df is None:
            if self._df_to_csv is None:
                print("couldn't find the file to be saved")
                return 0
            df = self._df_to_csv
        if filepath is None:
            filepath = "Saved_File.csv"
        df.to_csv(filepath)
