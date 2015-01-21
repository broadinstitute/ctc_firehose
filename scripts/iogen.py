import subprocess as sp
import os, sys

def openFile(filepath):
    if sys.platform == 'win32':
        os.startfile(filepath)
    elif sys.platform == 'linux2':
        sp.call(('xdg-open', filepath))
    elif sys.platform == 'darwin':
        sp.call(('open', filepath))

def getOutFilePaths(files_in, files_out, out_suffix):
    if not files_out:
        files_out = []
        for file_in in files_in:
            with open(file_in) as fin:
                files_out.append(os.path.splitext(file_in)[0] + '.' + out_suffix)
    return files_out

class Inputs(object):
    def __init__(self, args, out_suffix):
        self.files_in = args.files_in
        self.files_out = getOutFilePaths(
            args.files_in, args.files_out, out_suffix)
        self.t = args.t

def parameterGroup(s):
    try:
        return s.split(',')
    except:
        raise TypeError