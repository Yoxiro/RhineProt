import pandas
def FastaRead(Filepath:str)->pandas.DataFrame:
    """
    Read Protein_Sequences from fasta file and Return the DataFrame with "header","sequence","sequence_length"
    :param Filepath:
    :return:
    """
    records = []
    with open(Filepath, 'r') as file:
        header = None
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if header is not None:
                    records.append({'header': header.strip(">"),
                                    'sequence': sequence,
                                    'sequence_length':len(sequence)})
                header = line.strip()
                sequence = ''
            else:
                sequence += line.strip()
                # 添加最后一个记录
        if header is not None:
            records.append({'header': header.strip(">"),
                            'sequence': sequence,
                            'sequence_length':len(sequence)})

            # 将记录转换为DataFrame
    df = pandas.DataFrame(records)
    return df
