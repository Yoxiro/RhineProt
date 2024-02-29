import pandas
def FastaRead(Filepath:str)->pandas.DataFrame:
    """
    Read Protein_Sequences from fasta file
    :param Filepath:
    :return:
    """
    fa_dict = {}
    with open(Filepath,"r") as fa:
        for line in fa:
            # 去除末尾换行符
            line = line.replace('\n', '')
            if line.startswith('>'):
                # 去除 > 号
                seq_name = line[1:]
                fa_dict[seq_name] = ''
            else:
                # 去除末尾换行符并连接多行序列
                fa_dict[seq_name] += line.replace('\n', '')
    sequence_dataframe = pandas.DataFrame(fa_dict)
    return sequence_dataframe