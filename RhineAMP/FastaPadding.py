import pandas


def FastaPadding(fa:pandas.DataFrame):
    """
    pad the fa with "-"
    :param fa:
    :return:
    """
    max_length = max(fa["sequence_length"])
    pad = ["-"*(max_length - fa.iloc[i,2]) for i in range(fa.shape[0])]
    pad_Series = pandas.Series(pad)
    fa.iloc[:,1] = fa.iloc[:,1]+pad_Series
    return fa



