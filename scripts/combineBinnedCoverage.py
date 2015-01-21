import argparse as ap
import iogen, iocsv
import pandas as pd
import os, sys

OUT_SUFFIX = 'wig.combined.xls'
KEYS = ['contig', 'position']

def parse_args():
    stats_suffix = 'wig'
    os.environ['COLUMNS'] = '110'
    parser = ap.ArgumentParser(description='Combine .wig files into one spreadsheet')
    parser.add_argument('files_in', metavar='input_file,input_file,...',
        nargs='?', type=iogen.parameterGroup, help='Path(s) to *.' + stats_suffix)
    parser.add_argument('-o', metavar='output_file',
        dest='file_out', type=iogen.parameterGroup,
        help='Path to output Excel file (will be set to *.' + OUT_SUFFIX + ' by default)')
    parser.add_argument('-e', action='store_true',
        help='Open output Excel file on exit')
    return parser.parse_args()

def processData(data, dataCombined):
    
    bounds = data[KEYS[1]]
    bounds = bounds[bounds.str.contains('^variableStep')]
    bounds.replace('.*chr([0-9XY]*) .*', r'\1', regex=True, inplace=True)
    bounds.drop_duplicates(inplace=True)
    
    index = bounds.index.tolist() + [data.index[len(data) - 1]]
    data.dropna(inplace=True)
    bins = pd.cut(data.index, index)
    groups = data[KEYS[1]].groupby(bins)
    data[KEYS[0]] = groups.transform(lambda row: bounds[row.index.min() - 1])
    
    data = data.loc[:, KEYS + [data.columns[1]]]
    data[KEYS[1]] = data[KEYS[1]].astype(int, copy=False)
    
    if dataCombined is None:
        dataCombined = data
    else:
        dataCombined = dataCombined.merge(data, on=KEYS, how='outer')
    return dataCombined

def main():
    args = parse_args()
    print(args)
    inp = open(args.files_in[0])
    dataCombined = None
    for line in inp:
        lineContains = line.split('\t')
        file_in = lineContains[1]
        print("file name attempting to load:")
        print(file_in)
        sampleName = os.path.split(os.path.splitext(file_in)[0])[1]
        columns = [KEYS[1], sampleName]
        data = pd.read_csv(file_in, sep='\t', names=columns, skiprows=1)
        dataCombined = processData(data, dataCombined)
    file_out = os.path.join(args.file_out[0] + '.' + OUT_SUFFIX)
    dataCombined.to_excel(file_out, sheet_name="Coverage", index=False)
    sys.stdout.write(file_out)
    if args.e:
        iogen.openFile(file_out)

if __name__ == '__main__':
    main()