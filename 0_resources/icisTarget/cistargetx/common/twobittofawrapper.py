import os
import sys

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

from externalcallerror import ExternalCallError


def convert_to_fasta(feature_list, two_bit_filename, command):
    exec_command = [command, '-bed=stdin', two_bit_filename, 'stdout']
    try:
        procedure = subprocess.Popen(exec_command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        stdout = procedure.communicate(str(feature_list))[0]
    except OSError, msg:
        raise ExternalCallError("Error during execution of '" + ' '.join(exec_command) + "': " + str(msg))

    id2seq = dict()
    cur_seq = ''
    cur_id = None

    for line in stdout.split('\n'):
        if line.startswith('>'):
            if cur_id:
                id2seq[cur_id] = cur_seq
            cur_seq = ''
            cur_id = line.rstrip()[1:]
        else:
            cur_seq += line.rstrip()

    if cur_id:
        id2seq[cur_id] = cur_seq

    return id2seq

