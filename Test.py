import pandas as pd  
  
# 创建一个示例DataFrame  
data = {  
    'A': [4, 2, 9, 6],  
    'B': [1, 5, 3, 7],  
    'C': [8, 7, 2, 1]  
}  
  
df = pd.DataFrame(data)  
  
# 使用rank()方法为每个元素在其列中分配排名  
df_ranked = df.rank(ascending=True)  
  
print(type(df))