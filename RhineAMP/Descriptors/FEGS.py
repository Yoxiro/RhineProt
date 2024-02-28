#步骤
#1. 生成螺旋
#2. 氨基酸特征值排序
#3. 蛋白质覆盖
import math
import numpy

def _single()->numpy.ndarray:
    """

    :return:
    """
    return numpy.asarray([[math.cos(2*math.pi/20*i),
                               math.sin(2*math.pi/20*i),
                               1] for i in range(1,21)],
                       dtype=float).transpose()
def _pairs(single:numpy.ndarray) ->numpy.ndarray:
    after = numpy.repeat(single[:, :,numpy.newaxis], repeats=20, axis=2)
    after_T=after.transpose((0,2,1))
    pairs = (after*3+after_T)/4
    return pairs

def _psai(Protein_Sequence:str) -> numpy.ndarray:
    single = _single()
    pairs = _pairs(single=single)
    
def deasd
if __name__ == "__main__":
    L = [math.cos(2*math.pi/20*i) for i in range(1,21)]
    X = [[0 for i in range(20)] for j in range(20)]
    for i in range(1,21):
        for j in range(1,21):
            X[i-1][j-1] = L[i-1]*0.75+L[j-1]*0.25
    X=numpy.array(X)
