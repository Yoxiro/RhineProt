# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 14:59:15 2024

@author: Youzijie
"""
import os
from typing import Dict,List
os.chdir(r"F:\code\RhineAMP\RhineAMP\Descriptors")
with open("aaindex1.txt","r") as f:
    lines = f.readlines()
new_lines = []
for item in lines:
    new_lines.append(item.strip())
lines = new_lines

line = []
line.append([lines[1].split(sep=" ")[-1]])
for i in range(1,len(lines)):
    if lines[i] =="//":
        line[-1].append(lines[i-2])
        line[-1].append(lines[i-1])
        try:
            line.append([lines[i+1].split(sep=" ")[-1]])
        except:
            print("done")
new_line = []
for i in line:
    i1sp = i[1].split()
    i2sp = i[2].split()
    if "NA" in i1sp or "NA" in i2sp:
        print(i)
        continue
    new_line.append([i[0]])
    tov1 = [eval(k) for k in i1sp]
    tov2 = [eval(k) for k in i2sp]
    tov_sum = tov1+tov2
    new_line[-1].append(tov_sum)
line_dict ={item[0]:item[1] for item in new_line}

from typing import Dict, List, Set  
  
def remove_keys_with_duplicate_lists(dict_of_lists: Dict[str, List[float]]) -> Dict[str, List[float]]:  
    # 创建一个新的字典来存储没有重复元素的列表的项  
    new_dict = {}  
  
    # 遍历原始字典的键和值  
    for key, value_list in dict_of_lists.items():  
        # 如果列表中没有重复元素，则将其添加到新字典中  
        if len(value_list) == len(set(value_list)):  
            new_dict[key] = value_list  
  
    return new_dict

result = remove_keys_with_duplicate_lists(line_dict)

import pandas

index_1 = "A/L,R/K,N/M,D/F,C/P,Q/S,E/T,G/W,H/Y,I/V"
index_2 = ["A","R","N","D","C","Q","E","G","H","I",
           "L","K","M","F","P","S","T","W","Y","V"]
re_p = pandas.DataFrame(result,index=index_2)
re_p.to_csv("AAINDEX.csv")
