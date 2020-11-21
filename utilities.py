# Utilities file


def colored(seq):
    """Give each nucleotide a different color"""
    bcolors = {
        'A': '\033[92m',
        'C': '\033[94m',
        'G': '\033[93m',
        'T': '\033[91m',
        'U': '\033[91m',
        'reset': '\033[0;0m'
    }

    temp_str = ""

    for nuc in seq:
        if nuc in bcolors:
            temp_str += bcolors[nuc] + nuc
        else:
            temp_str += bcolors['reset'] + nuc

    return temp_str + '\033[0;0m'


def read_text_file(file_path):
    """Open text file and read line by line"""
    with open(file_path, 'r') as f:
        return ''.join([line.strip() for line in f.readlines()])


def write_text_file(file_path, seq, mode='a'):
    """Write given string to given file; mode='a' by default"""
    with open(file_path, mode) as f:
        f.write(seq + '\n')


def convert_FASTA(file_path):
    """Convert FASTA file to dictionary."""
    with open(file_path, 'r') as f:
        label = ''
        res_dict = {}
        for line in f:
            if '>' in line:
                label = line.strip('\n')
                res_dict[label] = ''
            else:
                res_dict[label] += line.strip('\n')

    return res_dict

