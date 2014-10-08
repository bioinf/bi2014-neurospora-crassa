#!/usr/bin/env python

import tempfile
import argparse
import io
import os
import subprocess
import configparser
from Bio.Blast import NCBIXML

if not os.path.exists('sweet.ini'):
    initial_cfg = configparser.ConfigParser()
    initial_cfg['general'] = {}
    initial_cfg['genblasta'] = {}
    with open('sweet.ini', 'w') as configfile:
        initial_cfg.write(configfile)

cfg = configparser.ConfigParser()
cfg.read('sweet.ini')

general_cfg = cfg['general']
genblasta_cfg = cfg['genblasta']


def dbname2path(dbname):
    return dbname + '/' + dbname + '.fasta'

class Blast:
    def __init__(self):
        self.toolname = 'blast'

    @staticmethod
    def writeraw(dbname, proteinsfile, output):
        subprocess.call(['tblastn', '-db', dbname + '/' + dbname, '-query', proteinsfile], stdout=output)

    @staticmethod
    def writeparsable(dbname, proteinsfile, output):
        subprocess.call(['tblastn', '-db', dbname + '/' + dbname, '-query', proteinsfile, '-outfmt', '5'], stdout=output)


class Exonerate:
    def __init__(self):
        self.toolname = 'exonerate'

    @staticmethod
    def writeraw(dbname, proteinsfile, output):
        subprocess.call(['exonerate', '--model', 'protein2genome', '--query',
                         proteinsfile, '--target', dbname2path(dbname)], stdout=output)

    @staticmethod
    def writeparsable(dbname, proteinsfile, output):
        subprocess.call(['exonerate', '--model', 'protein2genome', '--showvulgar', 'no',
                         '--showalignment', 'no', '--showquerygff', 'yes', '--query',
                         proteinsfile, '--target', dbname2path(dbname)], stdout=output)


class Genewise:
    def __init__(self):
        self.toolname = 'genewise'

    @staticmethod
    def writeraw(dbname, proteinsfile, output):
        subprocess.call(['genewise', proteinsfile, dbname2path(dbname)], stdout=output)


class Genblasta:
    def __init__(self):
        self.toolname = 'genblasta'

    @staticmethod
    def _prepareenv():
        newenv = os.environ.clone()
        newenv['PATH'] = genblasta_cfg['genblasta']['oldblastbin'] + ':' + newenv['PATH']
        return newenv

    @staticmethod
    def writeraw(dbname, proteinsfile, output):
        subprocess.call(['genblasta', '-q', proteinsfile, '-t', dbname2path(dbname)], stdout=output,
                        env=Genblasta._prepareenv())



def runraw(args):
    filename = ('' if args.runname is None else args.runname + '.') + args.toolname + '.out'
    print('running ' + args.toolname + ' (raw output)...')
    with open(filename, 'w') as output:
        tools[args.toolname].writeraw(args.dbname, args.proteins, output)
    print(filename + ' was written')


def runbasicdata(args):
    filename = ('' if args.runname is None else args.runname  + '.') + args.toolname + '.parsable.out'
    print('running ' + args.toolname + ' (parsable output)...')
    with open(filename, 'w') as output:
        tools[args.toolname].writeparsable(args.dbname, args.proteins, output)
    print(filename + ' was written')


tools = dict()
for x in [Blast(), Genewise(), Genblasta(), Exonerate()]:
    tools[x.toolname] = x

parser = argparse.ArgumentParser(description="Sweet bioinformatics tool.")
parser.add_argument('operation', choices=['raw', 'parsable'], help='Run mode')
parser.add_argument('toolname', choices=tools.keys(), help='Name of a tool to be used')
parser.add_argument('proteins', help='Input proteins file name')
parser.add_argument('dbname', help='Sequence database name (created with initialize.py)')
parser.add_argument('--runname', dest='runname', help='Run name (saving output to runname.%someextension%)')

args = parser.parse_args()

if args.operation == 'raw':
    runraw(args)
elif args.operation == 'parsable':
    runbasicdata(args)