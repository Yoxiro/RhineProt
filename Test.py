import pandas
import numpy as np

from RhineAMP.MachineLearning.RF import RF
from RhineAMP.MachineLearning.MLP import MLP
from RhineAMP.MachineLearning.KNN import KNN
from RhineAMP.MachineLearning.SVM import SVM
from RhineAMP.MachineLearning.LR import LR

positive = pandas.read_csv("Result/FEGS_Positive.csv",
                           index_col=0).T
negative = pandas.read_csv("Result/FEGS/Negative_Sample_1.csv",
                           index_col=0).T
print(positive.shape)
X = pandas.concat([positive, negative], axis=0)
X = np.array(X)
X.tofile("ex.csv", sep=",")
y = np.array([0] * positive.shape[0] + [1] * negative.shape[0])
# RF(X, y, fold=10)
MLP(X, y, fold=10)
# KNN(X, y, fold=10)
# LR(X, y, fold=10)
# SVM(X, y, fold=10)
