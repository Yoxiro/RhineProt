import pandas
import numpy as np

from RhineAMP.MachineLearning.RF import RF

positive = pandas.read_csv("Result/FEGS_Positive.csv",
                           index_col=0).T
print("Po shape")
print(positive.shape)

negative = pandas.read_csv("Result/FEGS/Negative_Sample_1.csv",
                           index_col=0).T
print("Po shape")
print(positive.shape)
print(positive)
print(negative)
X = pandas.concat([positive, negative], axis=0)
print(X.shape)
X = np.array(X)
X.tofile("ex.csv", sep=",")
y = np.array([0] * positive.shape[0] + [1] * negative.shape[0])
RF(X, y)
