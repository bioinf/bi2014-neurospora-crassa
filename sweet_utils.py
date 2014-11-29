import os


def dbpath2seq(dbpath):
    return os.path.join(dbpath, 'sequence.fasta')


def get_run_id(run_path):
    return os.path.basename(run_path)[:-len('.parseable.out')]


def get_tool_id(run_path):
    return get_run_id(run_path).split('.')[-1]

def get_run_format(run_path):
    return {'blast': 'blast-xml', 'exonerate': 'exonerate-vulgar'}[get_tool_id(run_path)]

searchio_formats = {'blast': 'blast-xml', 'exonerate': 'exonerate-vulgar'}
