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
    def getraw(dbname, proteinsfile, output):
        subprocess.call(['tblastn', '-db', dbname + '/' + dbname, '-query', proteinsfile], stdout=output)

    @staticmethod
    def getbasicdata(dbname, proteinsfile):
        with tempfile.TemporaryFile(mode='r+') as tmpfile:
            subprocess.call(['tblastn', '-db', dbname + '/' + dbname, '-query', proteinsfile], stdout=tmpfile)
            for record in NCBIXML.parse(tmpfile):






class Exonerate:
    def __init__(self):
        self.toolname = 'exonerate'

    @staticmethod
    def getraw(dbname, proteinsfile, output):
        subprocess.call(['exonerate', '--model', 'protein2genome', '--query',
                         proteinsfile, '--target', dbname2path(dbname)], stdout=output)

    @staticmethod
    def getbasicdata(dbname, proteinsfile):
        with tempfile.TemporaryFile(mode='r+') as tmpfile:
            subprocess.call(['exonerate', '--model', 'protein2genome', '--showquerygff', '--query',
                             proteinsfile, '--target', dbname2path(dbname)], stdout=tmpfile)
            tmpfile.seek(0)
            return tmpfile.read()


class Genewise:
    def __init__(self):
        self.toolname = 'genewise'

    @staticmethod
    def getraw(dbname, proteinsfile, output):
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
    def getraw(dbname, proteinsfile, output):
        subprocess.call(['genblasta', '-q', proteinsfile, '-t', dbname2path(dbname)], stdout=output,
                        env=Genblasta._prepareenv())



def runraw(args):
    print('running ' + args.toolname + ' with raw output...')
    with open(args.toolname + '.out', 'w') as output:
        tools[args.toolname].getraw(args.dbname, args.proteins, output)
    print(args.toolname + '.out was written')


def runbasicdata(args):
    print('running ' + args.toolname + ' with basic data...')
    print(tools[args.toolname].getbasicdata(args.dbname, args.proteins))

tools = dict()
for x in [Blast(), Genewise(), Genblasta(), Exonerate()]:
    tools[x.toolname] = x

parser = argparse.ArgumentParser(description="Sweet bioinformatics tool.")
parser.add_argument('operation', choices=['raw', 'basicdata'], help='Run mode')
parser.add_argument('toolname', choices=tools.keys(), help='Name of a tool to be used')
parser.add_argument('proteins', help='Input proteins file name')
parser.add_argument('dbname', help='Sequence database name (created with initialize.py)')

args = parser.parse_args()

if args.operation == 'raw':
    runraw(args)
elif args.operation == 'basicdata':
    runbasicdata(args)