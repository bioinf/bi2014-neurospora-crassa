#!/usr/bin/env python

import argparse
import os
from os import path
import shutil
import subprocess
import sweetconfig
import tempfile
import genblastA_to_gff3
from sweet_utils import dbpath2seq


class Blast:
    def __init__(self):
        self.toolname = 'blast'

    @staticmethod
    def writeraw(dbpath, proteinspath, outpath):
        subprocess.check_call(['tblastn', '-db', path.join(dbpath, path.basename(dbpath)), '-query', proteinspath, '-out', outpath])

    @staticmethod
    def writeparseable(dbpath, proteinspath, outpath):
        subprocess.check_call(['tblastn', '-db', path.join(dbpath, path.basename(dbpath)), '-query', proteinspath, '-outfmt', '5', '-out', outpath])


class Exonerate:
    def __init__(self):
        self.toolname = 'exonerate'

    @staticmethod
    def writeraw(dbpath, proteinspath, outpath):
        with open(outpath, 'w') as outfile:
            subprocess.check_call(
                    ['exonerate', '--model', 'protein2genome', '--query',
                    proteinspath, '--softmasktarget', 'yes', 
                    '--target', dbpath2seq(dbpath),
                    '--ryo', '<< %pi >>\n'], stdout=outfile)

    @staticmethod
    def writeparseable(dbpath, proteinspath, outpath):
        Exonerate.writeraw(dbpath, proteinspath, outpath)

class Genewise:
    def __init__(self):
        self.toolname = 'genewise'

    @staticmethod
    def writeraw(dbpath, proteinspath, outpath):
        pass

class Genblasta:
    def __init__(self):
        self.toolname = 'genblasta'

    @staticmethod
    def _prepareenv():
        newenv = os.environ.copy()
        newenv['PATH'] = sweetconfig.genblasta_cfg['oldblastbin'] + ':' + newenv['PATH']
        return newenv

    @staticmethod
    def _prepare_aux_dir(proteinspath, outpath):
        outpath_abs = path.abspath(outpath)

        auxdirpath = path.join(path.dirname(outpath_abs), '.aux.' + path.basename(outpath_abs))
        if not path.exists(auxdirpath):
            os.mkdir(auxdirpath)

        proteinspath_aux = path.abspath(path.join(auxdirpath, 'query.fasta'))
        shutil.copyfile(proteinspath, proteinspath_aux)

        return auxdirpath, proteinspath_aux


    @staticmethod
    def writeraw(dbpath, proteinspath, outpath):
        auxdir, proeins_aux = Genblasta._prepare_aux_dir(proteinspath, outpath)
        subprocess.check_call(['genblasta', '-q', os.path.abspath(proeins_aux), '-t', os.path.abspath(dbpath2seq(dbpath))
                               , '-o', os.path.abspath(outpath)], cwd=auxdir, env=Genblasta._prepareenv())

    @staticmethod
    def writeparseable(dbpath, proteinspath, outpath):
        Genblasta.writeraw(dbpath, proteinspath, outpath)


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

if __name__ == '__main__':
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
