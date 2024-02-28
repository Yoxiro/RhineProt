# 步骤
# 1. 生成螺旋
# 2. 氨基酸特征值排序
# 3. 蛋白质覆盖
import math
import numpy
import pandas

from pkg_resources import resource_filename

from RhineAMP.Descriptors.AAC import AAC
from RhineAMP.Descriptors.DPC import DPC_2D_array


def _single() -> numpy.ndarray:
    """

    :return:
    """
    return numpy.asarray([[math.cos(2 * math.pi / 20 * i),
                           math.sin(2 * math.pi / 20 * i),
                           1] for i in range(1, 21)],
                         dtype=float).transpose()


def _pairs(single: numpy.ndarray) -> numpy.ndarray:
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
    psai_after = psai[:, 1:] - psai[:, :-1]
    psai_after = psai_after * psai_after
    psai_after = psai_after.sum(axis=0)
    psai_sqrt = numpy.sqrt(psai_after)[1:]
    return psai_sqrt


def _sort(file=None) -> pandas.DataFrame:
    """"""
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
    _sort()
