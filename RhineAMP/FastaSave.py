import pandas


def FastaSave(fa: pandas.DataFrame, filename: str):
    length = fa.shape[0]
    with open(filename, "w", encoding="utf-8") as f:
        for i in range(length):
            f.write(">" + fa.iloc[i, 0] + "\n")
            f.write(fa.iloc[i, 1] + "\n")
