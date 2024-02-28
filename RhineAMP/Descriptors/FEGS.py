# 步骤
# 1. 生成螺旋
# 2. 氨基酸特征值排序
# 3. 蛋白质覆盖
import math
import numpy
import pandas

from pkg_resources import resource_filename


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


def _psai(Protein_Sequence: str) -> numpy.ndarray:
    single = _single()
    pairs = _pairs(single=single)


def _sort(file=None) -> pandas.core.frame.DataFrame:
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
