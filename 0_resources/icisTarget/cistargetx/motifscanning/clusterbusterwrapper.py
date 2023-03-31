import os
import re
import sys

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

from cistargetx.common.externalcallerror import ExternalCallError
from cistargetx.conversion.featurelist import Feature


def scan(motifs_filename, fasta_filename, command, locations, cluster_threshold=5.0, motif_threshold=6.0):
    # Running Cluster Buster.
    statement = [command, '-f0', '-c{0:f}'.format(cluster_threshold), '-m{0:f}'.format(motif_threshold),
                 motifs_filename, fasta_filename]
    try:
        stdout, stderr = subprocess.Popen(statement, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    except OSError, msg:
        # CAVE: do not check stderr, because this stream is used by Cluster Buster for outputting progress information.
        raise ExternalCallError("Error during execution of '" + ' '.join(statement) + "': " + str(msg))

    class MotifFeature(object):
        def __init__(self, name, interval, strand, score):
            self.name = name
            self.interval = interval
            self.strand = strand
            self.score = score

    class ClusterFeature(object):
        UNKNOWN_LOCATION = Feature("?", 0, 0)

        def __init__(self, sequence_id, number):
            self.sequence_id = sequence_id
            self.number = number
            self.interval = tuple()
            self.score = float('nan')
            self.motifs = []

        def add_motif(self, motif):
            self.motifs.append(motif)

        @property
        def crm_feature(self):
            feature = locations.get(self.sequence_id, ClusterFeature.UNKNOWN_LOCATION)
            offset = feature.interval[0]

            # Interval returned by Cluster Buster is a closed interval. The locations are 1-based.
            interval = (offset + self.interval[0] - 1, offset + self.interval[1])

            return Feature(feature.chromosome, interval[0], interval[1], "CRM", self.score)

        @property
        def motif_features(self):
            feature = locations.get(self.sequence_id, ClusterFeature.UNKNOWN_LOCATION)
            offset = feature.interval[0]

            for motif in self.motifs:
                # Interval returned by Cluster Buster is a closed interval. The locations are 1-based.
                interval = (offset + motif.interval[0] - 1, offset + motif.interval[1])

                yield Feature(feature.chromosome, interval[0], interval[1], motif.name, motif.score, motif.strand)

        @property
        def features(self):
            result = [self.crm_feature]
            result.extend(self.motif_features)

            return result

    features = []
    cur_sequence_id = None
    cur_cluster = None

    for line in stdout.split('\n'):
        if line.startswith(">"):
            cur_sequence_id = line[1:line.rindex("(") - 1]
        elif re.match("^CLUSTER\s+(\d+)$", line):
            if cur_cluster:
                features.extend(cur_cluster.features)

            cur_cluster = ClusterFeature(cur_sequence_id, int(re.match("^CLUSTER\s+(\d+)$", line).group(1)))
        elif line.startswith("Location:"):
            cur_cluster.interval = map(int, re.match("^Location:\s+(\d+)\s+to\s+(\d+)$", line).groups())
        elif line.startswith("Score:"):
            cur_cluster.score = float(re.match("^Score:\s+(\d+(\.\d+)?(e(-)?(\d+))?)$", line).group(1))
        elif re.match("^(\S+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+(\.\d+)?)\s+\S+$", line):
            match = re.match("^(\S+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\d+(\.\d+)?)\s+\S+$", line)
            motif_name = match.group(1)
            motif_interval = (int(match.group(2)), int(match.group(3)))
            motif_strand = match.group(4)
            motif_score = float(match.group(5))

            cur_cluster.add_motif(MotifFeature(motif_name, motif_interval, motif_strand, motif_score))
        else:
            pass

    if cur_cluster:
        features.extend(cur_cluster.features)

    return features
