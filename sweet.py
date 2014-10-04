import argparse
import os
import subprocess
import configparser
import sys

if not os.path.exists('sweet.ini'):
    initial_cfg = configparser.ConfigParser()
    initial_cfg['general'] = { }
    initial_cfg['blast'] = {'dbname' : './blastdb/neurospora'}
    initial_cfg['genblasta'] = {}
    initial_cfg['exonerate'] = {}
    with open('sweet.ini', 'w') as configfile:
        initial_cfg.write(configfile)

cfg = configparser.ConfigParser()
cfg.read('sweet.ini')

general_cfg = cfg['general']
blast_cfg = cfg['blast']
genblasta_cfg = cfg['genblasta']
exonerate_cfg = cfg['exonerate']

def do_initblastdb(args):
    if not os.path.exists(os.path.dirname(args.dbname)):
        os.mkdir(os.path.dirname(args.dbname))
    print("Initializing BLAST database...")

    subprocess.call(['makeblastdb', '-in', args.seq, '-dbtype', 'nucl', '-out', args.dbname])

    print("Initialization BLAST database... finished.")

def do_search_genblasta(args):
    if genblasta_cfg.get('oldblastbin') is None:
        sys.exit("'oldblastbin' option for genblasta should be set in sweet.ini .")
    if args.seq is None:
         sys.exit('ERROR: Sequence argument is missing.')

    old_blast_env = os.environ.copy()
    old_blast_env["PATH"] = genblasta_cfg['oldblastbin'] + ":" + old_blast_env["PATH"]

    print("Running genblasta...")

    subprocess.call(['genblasta', '-q', args.prot, '-t', args.seq, '-o', 'genblasta.out'], env=old_blast_env)

    print("Writing to genblasta.out...")

def do_search_blast(args):
    if args.dbname is None:
        sys.exit('ERROR: Blast DB name argument is missing.')

    print("Running blast...")

    subprocess.call(['tblastn', '-query', args.prot, '-db', args.dbname, '-out', 'blast.out'])

    print("Writing to blast.out...")

def do_search_exonerate(args):
    if args.seq is None:
        sys.exit('ERROR: Sequence argument is missing.')

    print("Running exonerate...")

    with open('exonerate.out', 'w') as outfile:
        subprocess.call(['exonerate', '--model', 'protein2genome', args.prot, args.seq], stdout=outfile)

    print("Writing to exonerate.out..")

def do_search(args):
    print("Searching proteins...")
    if args.method == 'genblasta':
        do_search_genblasta(args)
    if args.method == 'blast':
        do_search_blast(args)
    if args.method == 'exonerate':
        do_search_exonerate(args)
    print("Searching proteins... finished")

parser = argparse.ArgumentParser(description="BI sweet annotation.")
subparsers = parser.add_subparsers(metavar='action', dest='action', help='What should I do? (initblastdb, search)')

def add_sequence_arg(parser, required = False):
    parser.add_argument('-s', '--sequence', metavar='sequence', dest='seq', default=general_cfg.get('sequence'),
                        required = general_cfg.get('sequence') is None if required else False,
                        help='Path to sequence in fasta format (required for genblasta and exonerate)')

def add_dbname_arg(parser, required = False):
    parser.add_argument('-d', '--dbname', metavar='dbname', dest='dbname', default=blast_cfg.get('dbname'),
                        required = blast_cfg.get('dbname') is None if required else False,
                        help='Database path (required for blast)')

parser_init = subparsers.add_parser('initblastdb')
add_sequence_arg(parser_init, True)
add_dbname_arg(parser_init, True)

parser_search = subparsers.add_parser('search')
parser_search.add_argument('method', choices=['genblasta', 'blast', 'exonerate'])
parser_search.add_argument('-p', '--proteins', metavar='proteins', dest='prot', required=general_cfg.get('proteins') is None,
                            default=general_cfg.get('proteins'), help='Path to proteins in fasta format')
add_sequence_arg(parser_search)
add_dbname_arg(parser_search)

args = parser.parse_args()
if args.action == 'initblastdb':
    do_initblastdb(args)
if args.action == 'search':
    do_search(args)