import pandas as pd

def printline(fout, dict, fieldnames):
    l = len(fieldnames) - 1
    for i, key in enumerate(fieldnames):
        fout.write(dict[key])
        if i < l:
            fout.write('\t')
        else:
            fout.write('\n')

def getDataCSV(file_in, file_out, sep='\t', skipheader=True, comment='#'):
    if skipheader and file_out:
        with open(file_in) as fin:
            with open(file_out, 'w') as fout:
                fout.write(fin.readline())
    return pd.read_csv(file_in, sep=sep, comment=comment)