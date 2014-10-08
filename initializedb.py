#!/usr/bin/env python

import argparse
import shutil
import subprocess
import os

parser = argparse.ArgumentParser(description='Sweet database initialisation.')
parser.add_argument(metavar='name', dest='name', help='Database name')
parser.add_argument(metavar='sequence', dest='sequence', help='Fasta sequence file.')

args = parser.parse_args()
os.mkdir(args.name)
shutil.copyfile(args.sequence, args.name + '/' + args.name + '.fasta')
subprocess.call(['makeblastdb', '-in', args.sequence, '-dbtype', 'nucl', '-out', args.name + "/" + args.name])
