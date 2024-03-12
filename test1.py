import pandas
df = pandas.read_csv("Result/FEGS_Negative.csv",index_col=0)
print(df.shape)
for i in range(10):
    path = "Result/FEGS/Negative_sample_"+str(i)+".csv"
    res = df.sample(n = 1253,axis=1)
    print(res.shape)
    res.to_csv(path)