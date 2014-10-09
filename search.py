#!/usr/bin/env python

import argparse
import os
import shutil
import subprocess
import sweetconfig

def dbpath2seq(dbpath):
    return os.path.join(dbpath, 'sequence.fasta')

class Blast:
    def __init__(self):
        self.toolname = 'blast'

    @staticmethod
    def writeraw(dbpath, proteinspath, outpath):
        subprocess.call(['tblastn', '-db', os.path.join(dbpath, os.path.basename(dbpath)), '-query', proteinspath, '-out', outpath])

    @staticmethod
    def writeparseable(dbpath, proteinspath, outpath):
        subprocess.call(['tblastn', '-db', os.path.join(dbpath, os.path.basename(dbpath)), '-query', proteinspath, '-outfmt', '5', '-out', outpath])


class Exonerate:
    def __init__(self):
        self.toolname = 'exonerate'

    @staticmethod
    def writeraw(dbpath, proteinspath, outpath):
        with open(outpath, 'w') as outfile:
            subprocess.call(['exonerate', '--model', 'protein2genome', '--query',
                             proteinspath, '--target', dbpath2seq(dbpath)], stdout=outfile)

    @staticmethod
    def writeparseable(dbpath, proteinspath, outpath):
        with open(outpath, 'w') as outfile:
            subprocess.call(['exonerate', '--model', 'protein2genome', '--query', proteinspath,
                             '--showalignment', 'no', '--target', dbpath2seq(dbpath)], stdout=outfile)


class Genewise:
    def __init__(self):
        self.toolname = 'genewise'

    @staticmethod
    def writeraw(dbpath, proteinspath, outpath):
        with open(outpath, 'w') as outfile:
            subprocess.call(['genewise', proteinspath, dbpath2seq(dbpath)], stdout=outfile)


class Genblasta:
    def __init__(self):
        self.toolname = 'genblasta'

    @staticmethod
    def _prepareenv():
        newenv = os.environ.copy()
        newenv['PATH'] = sweetconfig.genblasta_cfg['oldblastbin'] + ':' + newenv['PATH']
        return newenv

    @staticmethod
    def writeraw(dbpath, proteinspath, outpath):
        seqpath_abs = os.path.abspath(dbpath2seq(dbpath))
        outpath_abs = os.path.abspath(outpath)

        auxdirpath = os.path.join(os.path.dirname(outpath_abs), '.aux.' + os.path.basename(outpath_abs))
        if not os.path.exists(auxdirpath):
            os.mkdir(auxdirpath)

        proteinspath_aux = os.path.abspath(os.path.join(auxdirpath, 'query.fasta'))
        shutil.copyfile(proteinspath, proteinspath_aux)

        subprocess.call(['genblasta', '-q', proteinspath_aux, '-t', seqpath_abs, '-o', outpath_abs],
                        cwd = auxdirpath, env=Genblasta._prepareenv())


def runraw(args):
    filename = ('' if args.runname is None else os.path.basename(args.runname) + '.') + args.toolname + '.out'
    print('running ' + args.toolname + ' (raw output)...')
    tools[args.toolname].writeraw(args.dbname, args.proteins, filename)
    print(filename + ' was written')


def runparseable(args):
    filename = ('' if args.runname is None else os.path.basename(args.runname) + '.') + args.toolname + '.parseable.out'
    print('running ' + args.toolname + ' (parseable output)...')
    tools[args.toolname].writeparseable(args.dbname, args.proteins, filename)
    print(filename + ' was written')


tools = dict()
for x in [Blast(), Genewise(), Genblasta(), Exonerate()]:
    tools[x.toolname] = x

parser = argparse.ArgumentParser(description="Sweet bioinformatics tool.")
parser.add_argument('operation', choices=['raw', 'parseable'], help='Run mode')
parser.add_argument('toolname', choices=tools.keys(), help='Name of a tool to be used')
parser.add_argument('proteins', help='Input proteins file name')
parser.add_argument('dbname', help='Sequence database name (created with initialize.py)')
parser.add_argument('--runname', dest='runname', help='Run name (saving output to runname.%someextension%)')

args = parser.parse_args()

if args.operation == 'raw':
    runraw(args)
elif args.operation == 'parseable':
    runparseable(args)