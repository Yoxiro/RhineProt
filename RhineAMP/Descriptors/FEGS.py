import math
import numpy
import pandas

from pkg_resources import resource_filename

from RhineAMP.Descriptors.AAC import AAC
from RhineAMP.Descriptors.DPC import DPC_2D_array, DPC

"""
import os 
os.chdir(r"F:\code\RhineAMP\RhineAMP\Descriptors")
"""


def _single() -> numpy.ndarray:
    """

    :return:
    """
    return numpy.asarray([[math.cos(2 * math.pi / 20 * i),
                           math.sin(2 * math.pi / 20 * i),
                           1] for i in range(1, 21)],
                         dtype=float).transpose()


def _pairs(single: numpy.ndarray) -> numpy.ndarray:
    """

    :param single:
    :return:
    """
    after = numpy.repeat(single[:, :, numpy.newaxis], repeats=20, axis=2)
    after_T = after.transpose((0, 2, 1))
    pairs = (after * 3 + after_T) / 4
    return pairs


def _pairs_sum(Protein_Sequence: str, Squence: pandas.Series) -> numpy.ndarray:
    """

    :param Protein_Sequence:
    :param Squence: AAINDEX_Seq.iloc[:,i]
    :return:
    """
    # to generate pairs
    single = _single()
    pairs = _pairs(single=single)
    # to generate dpc in 3 20 20
    dpc = DPC_2D_array(Protein_Sequence)
    dpc = dpc.reindex(index=Squence)
    dpc = dpc[Squence]
    dpc_array = numpy.array(dpc)
    dpc_array = numpy.repeat(dpc_array[numpy.newaxis, :, :], repeats=3, axis=0)
    # to generate sum pairs
    sum_pairs = dpc_array * pairs
    sum_pairs = sum_pairs.sum(axis=1).sum(axis=1)
    return sum_pairs


def _psai(Protein_Sequence: str, Sequence: pandas.Series) -> numpy.ndarray:
    """

    :param Protein_Sequence:
    :param Sequence:
    :return:
    """
    single = _single()
    psai = numpy.array([[0] * (len(Protein_Sequence) + 1)] * 3,
                       dtype=float)
    pairs_sum = _pairs_sum(Protein_Sequence, Sequence)
    for i in range(len(Protein_Sequence)):
        position = numpy.array(Protein_Sequence[i] == Sequence)
        position = numpy.repeat(position[numpy.newaxis, :], repeats=3, axis=0)
        phi = position * single
        phi = phi.sum(axis=1)
        psai[:, i + 1] = psai[:, i] + phi + pairs_sum
    return psai


def _euclidean(psai: numpy.ndarray) -> numpy.ndarray:
    """

    :param psai:
    :return:
    """
    repeat_times = psai.shape[1]
    psai_2d = numpy.repeat(psai[:, :, numpy.newaxis], repeats=repeat_times, axis=2)
    psai_2d_T = psai_2d.transpose((0, 2, 1))
    psai_minus = psai_2d - psai_2d_T
    psai_sqr = psai_minus * psai_minus
    psai_sqr_add = psai_sqr.sum(axis=0)
    psai_sqrt = numpy.sqrt(psai_sqr_add)[1:, 1:]
    return psai_sqrt


def _generating_M(Protein_Sequence: str, Sequence: pandas.Series):
    """

    :param Protein_Sequence:
    :param Sequence:
    :return:
    """
    psai = _psai(Protein_Sequence, Sequence)
    euclidean = _euclidean(psai)
    sum_of_euclidean = _sum_of_euclidean(euclidean)
    M = euclidean / sum_of_euclidean
    return M


def eigenvalue(M: numpy.ndarray):
    """

    :param M:
    :return:
    """
    eigenvalues = numpy.linalg.eigvals(M)
    leading_eigenvalue = max(eigenvalues, key=abs)
    return leading_eigenvalue


def FEGS(Protein_Sequence: str) -> pandas.Series:
    """

    :param Protein_Sequence:
    :return:
    """
    AAINDEX_Seq = _sort()
    eigenvalue_list = []
    for i in range(AAINDEX_Seq.shape[1]):
        Sequence = AAINDEX_Seq.iloc[:, i]
        M = _generating_M(Protein_Sequence, Sequence)
        eigenvalues = eigenvalue(M)
        eigenvalue_list.append(eigenvalues)
    eigenvalue_Series = pandas.Series(eigenvalue_list,
                                      index=["FEGS" + str(i + 1) for i in range(len(eigenvalue_list))])
    eigenvalue_Series = eigenvalue_Series / len(Protein_Sequence)
    dpc = DPC(Protein_Sequence)
    aac = AAC(Protein_Sequence)
    fegs = pandas.concat([eigenvalue_Series, dpc, aac])
    return fegs


def _sum_of_euclidean(euclidean: numpy.ndarray) -> numpy.ndarray:
    """

    :param euclidean:
    :return:
    """
    sum_of_euclidean = euclidean.copy()
    index = sum_of_euclidean.shape[1]
    euclidean_offset = euclidean.diagonal(offset=1)
    for i in range(0, index):
        euclidean_offset_cumsum = euclidean_offset[i:].cumsum()
        sum_of_euclidean[i, i + 1:] = euclidean_offset_cumsum
    sum_of_euclidean = numpy.triu(sum_of_euclidean)
    sum_of_euclidean += sum_of_euclidean.T
    sum_of_euclidean[numpy.diag_indices(index)] = 1
    return sum_of_euclidean


def _sort(file=None) -> pandas.DataFrame:
    """

    :param file:
    :return:
    """
    if file is None:
        file = resource_filename(__name__, "Data/AAINDEX.csv")
    AAINDEX = pandas.read_csv(file, index_col=0)
    AAINDEX_Ranked = AAINDEX.rank()
    AAINDEX_Ranked.dtype = int
    AAINDEX_Seq = pandas.DataFrame(index=range(1, 21),
                                   columns=AAINDEX_Ranked.columns)
    for coli in range(0, AAINDEX_Ranked.shape[1]):
        AAINDEX_Series = AAINDEX_Ranked.iloc[:, coli]
        AAINDEX_Series_After = AAINDEX_Series.sort_values()
        AAINDEX_Seq.iloc[:, coli] = AAINDEX_Series_After.index
    return AAINDEX_Seq


if __name__ == "__main__":
    Protein_Sequence = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTESA"
    _sort()
