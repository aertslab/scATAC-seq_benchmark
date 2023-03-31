import tempfile
import os
import zipfile
import re

from cistargetx.common.twobittofawrapper import convert_to_fasta
import clusterbusterwrapper


class MotifCollection(object):
    @staticmethod
    def load_from_folder(foldername, extension='cb'):
        def _derive_motif_id(motif_filename):
            idx = motif_filename.rfind("." + extension)

            if idx < 0:
                return motif_filename

            return motif_filename[:idx]

        return MotifCollection(dict((_derive_motif_id(motif_filename), os.path.join(foldername, motif_filename))
                                    for motif_filename in os.listdir(foldername)
                                    if motif_filename.endswith("." + extension)))

    @staticmethod
    def load_from_zip_archive(zip_filename):
        motif_id2filename = dict()

        with zipfile.ZipFile(zip_filename, mode='r', allowZip64=True) as zip_input_fh:
            with zip_input_fh.open('motifs.tsv', 'r') as input_fh:
                for line in input_fh:
                    columns = re.split('[ \t]+', line.strip())
                    motif_id = columns[3]
                    motif_filename = 'singletons/{0:s}.cb'.format(motif_id)
                    motif_id2filename[motif_id] = motif_filename

        return MotifCollection(motif_id2filename, zip_filename)

    def __init__(self, motif_id2filename, zip_filename=None):
        self.motif_id2filename = motif_id2filename
        self.zip_filename = zip_filename

    def __len__(self):
        return len(self.motif_id2filename)

    def __str__(self):
        return "MotifCollection(#{0:d})".format(len(self))

    def reduce(self, motif_ids):
        set_motif_ids = set(motif_ids)

        return MotifCollection(dict((motif_id, motif_filename)
                                    for motif_id, motif_filename in self.motif_id2filename.iteritems()
                                    if motif_id in set_motif_ids),
                               self.zip_filename)

    def write(self, motifs_filename):
        if self.zip_filename:
            with open(motifs_filename, 'w') as output_fh:
                with zipfile.ZipFile(self.zip_filename, mode='r', allowZip64=True) as zip_input_fh:
                    for motif_id, motif_filename in self.motif_id2filename.iteritems():
                        with zip_input_fh.open(motif_filename, 'r') as input_fh:
                            output_fh.write(input_fh.read())
        else:
            with open(motifs_filename, 'w') as output_fh:
                for motif_id, motif_filename in self.motif_id2filename.iteritems():
                    with open(motif_filename, 'r') as input_fh:
                        output_fh.write(input_fh.read())


def scan(feature_list, motifs, cbust_command, two_bit_to_fa_command, genome_2bit_filename):
    fasta_filename = tempfile.mktemp()

    with open(fasta_filename, 'w') as output_fh:
        for sequence_id, sequence_data in convert_to_fasta(feature_list,
                                                           genome_2bit_filename,
                                                           two_bit_to_fa_command).iteritems():
            output_fh.write(">{0:s}\n{1:s}\n".format(sequence_id, sequence_data))

    motifs_filename = tempfile.mktemp()
    motifs.write(motifs_filename)

    features = clusterbusterwrapper.scan(motifs_filename=motifs_filename,
                                         fasta_filename=fasta_filename,
                                         command=cbust_command,
                                         locations=feature_list,
                                         cluster_threshold=0.0,
                                         motif_threshold=3.0)

    os.remove(fasta_filename)
    os.remove(motifs_filename)

    for feature in features:
        yield feature
