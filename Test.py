import pandas as pd

# 创建一个示例DataFrame
df = pd.DataFrame({
    'A': [1, 2, 3],
    'B': [4, 5, 6],
    'C': [7, 8, 9]
})

# 使用transpose()方法进行转置
df_transposed = df.transpose()
print(df)
print(df_transposed)