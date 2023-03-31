import logging
import operator
import os
import re
import sys
import tempfile

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess


COLORS = ['51CC8C', '51CCCC', '337F7F', '8ECC51', '597F33', '8E51CC', 'CCAD51', '7F6C33',
          '51CC70', '337F46', '5170CC', 'CC51AD', '7F336C', 'CC7F51', '7F4F33', 'CC5151',
          'BCCC51', '757F33', '60CC51', '3C7F33', '51CC9E', '337F62', '519ECC', '33627F',
          '6051CC', 'BC51CC', 'CC517F', 'CC6851', 'CC9651', '7F5E33', 'CCC451', '7F7A33',
          'A5CC51', '677F33', '77CC51', '4A7F33', '337F37', '51CC87', '337F54', '51CCB5',
          '337F71', '51B5CC', '33717F', '5187CC', '33547F', '5159CC', '7751CC', 'A551CC',
          'CC51C4', '7F337A', 'CC5196', 'CC5168', 'CC5D51', '7F3A33', 'CC7451', '7F4833',
          'CC8A51', 'CCB851', '9ACC51', '55CC51', '51CC92', '51C0CC', '517BCC', '6C51CC',
          'B151CC', 'CC51A1', 'CC515D', 'CC6E51', 'CC9051', 'CCB351', 'C2CC51', 'A0CC51',
          '7DCC51', '5BCC51', '51C6CC', '51A3CC', '7F335E', '7F3341']

NR_OF_COLORS = len(COLORS)

OUTPUT_CLUSTERS_FILENAME = "out_tree_clusters.txt"
OUTPUT_TREE_FILENAME = "out.tree"
OUTPUT_FBP_FILENAME = "outFBP.txt"


def get_color_string(idx):
    return COLORS[idx % NR_OF_COLORS]


def create_temp_filename(extension):
    return tempfile.mkstemp(suffix=extension, text=True)[1]


def create_temp_foldername(extension):
    return tempfile.mkdtemp(suffix=extension)


def reduce_database(motif_ids, full_motif_database_filename, reduced_motif_database_filename, has_clusters=None):
    """
    Create a motif database file in TRANSFAC format for usage with STAMP with only the motifs of interest.

    STAMP (v1.2 or lower) only parses the following specific TRANSFAC format correctly:
      - motif starts with tag "DE  motif_name".
      - "DE" tag is directly followed by matrix lines (no "P0" or "PO" tag allowed).
      - matrix lines are directly followed by a "XX" tag line.
      - if other tags are present between "DE" and matrix lines or between matrix lines and "XX",
        STAMP adds for each of those lines "0.0  0.0  0.0  0.0" to the matrix.

    DE  motif_name
    01     74     12     13      0      A
    02     73      3     19      3      A
    03      0     80      2     16      C
    04      2     83      0     13      C
    05     54     36      1      7      M
    06     27      3     64      5      G
    07     11     52      6     29      Y
    08     60     33      1      4      M
    XX

    """

    motif_ids = set(motif_ids)
    found_motifs = set()

    cur_motif = None
    matrix_start = False

    # Match "AC" tag lines.
    header_pattern = re.compile(r"^AC\s+(\S+)\s*$")
    # Match "P0  A  C  G  T" and "PO  A  C  G  T" lines.
    matrix_header_pattern = re.compile(r"^P[0O]\s+A\s+C\s+G\s+T\s*$")
    # Match matrix lines.
    matrix_pattern = re.compile(r"^([0-9]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)")

    with open(reduced_motif_database_filename, 'w') as output_fh:
        with open(full_motif_database_filename, 'r') as input_fh:
            for line in input_fh:
                line = line.rstrip()
                header_match = header_pattern.match(line)

                if header_match:
                    matrix_start = False
                    cur_motif = None

                    # Check if the motif name after the "AC" tag matches any of the wanted motifs.
                    if header_match.group(1) in motif_ids:
                        cur_motif = header_match.group(1)
                        found_motifs.add(cur_motif)

                        if len(found_motifs) == 1:
                            # For the first found motif, replace the "AC" tag with the "DE" tag.
                            print >> output_fh, 'DE' + line[2:]
                        else:
                            # For the rest of the motifs, replace the "AC" tag with the "DE" tag
                            # and add "XX" before this line.
                            print >> output_fh, 'XX\nDE' + line[2:]

                    continue

                if cur_motif:
                    # Check if this line contains the "P0" or "PO" tag.
                    matrix_header_match = matrix_header_pattern.match(line)

                    if matrix_header_match:
                        matrix_start = True

                        continue

                    # Check if this is a matrix line.
                    matrix_match = matrix_pattern.match(line)

                    if matrix_start and matrix_match:
                        # Print matrix line if the line contains 5 columns with numbers:
                        #   - column 1: matrix position
                        #   - column 2: A frequency (integer or floating point value)
                        #   - column 3: C frequency (integer or floating point value)
                        #   - column 4: G frequency (integer or floating point value)
                        #   - column 5: T frequency (integer or floating point value)
                        #   - column 6: optional consensus letter
                        print >> output_fh, line

            # Print "XX" after the last matrix.
            print >> output_fh, "XX"

    if len(found_motifs) != len(motif_ids):
        if has_clusters:
            logging.warning("STAMP: clustering {0:d} motifs found in motifCollection.".format(len(found_motifs)))
        else:
            logging.warning(
                "STAMP: only {0:d} out of {1:d} motifs were found in '{2:s}'.".format(len(found_motifs),
                                                                                      len(motif_ids),
                                                                                      full_motif_database_filename)
            )

    return motif_ids - found_motifs


def cluster_iterator(filename):
    # XX    Cluster_Members:    MA0335.1-MET4    PL0012.1-hlh-2__hlh-8    M00277-V-LMO2COM_01
    pattern = re.compile(r"^XX\s+Cluster_Members:\s+(.*)$")

    with open(filename, 'r') as input_fh:
        for line in input_fh:
            match = pattern.match(line)

            if not match:
                continue

            yield re.split('\s+', match.group(1).strip())


def cluster(full_motif_database_filename, score_distribution_filename, motif_names, cmd='STAMP', has_clusters=None):
    # Create subset of TRANSFAC database with only the requested motifs.
    reduced_motif_database_filename = create_temp_filename('.transfac')
    unknown_motif_names = reduce_database(motif_names, full_motif_database_filename, reduced_motif_database_filename,
                                          has_clusters)

    if (len(motif_names) - len(unknown_motif_names)) <= 1:
        logging.warning("STAMP: no motifs where found so skipping clustering.")
        os.remove(reduced_motif_database_filename)

        return [], list(unknown_motif_names)

    # Create temporary folder in which STAMP can write its result files.
    results_foldername = create_temp_foldername('.stamp')

    # Running STAMP.
    statement = [cmd, '-tf', reduced_motif_database_filename,
                 '-cc', 'SSD',  #'-align', 'SWA', '-go', '1.0', '-ge', '1.5', '-overlapalign'
                 '-sd', score_distribution_filename,
                 '-chp']
    try:
        stdout, stderr = subprocess.Popen(statement, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                          cwd=results_foldername).communicate()
    except OSError, msg:
        logging.error("STAMP: error during execution of '" + ' '.join(statement) + "': " + str(msg))
        sys.exit(1)

    os.remove(reduced_motif_database_filename)

    # Parse STAMP output.
    clusters = [motifs for motifs in cluster_iterator(os.path.join(results_foldername, OUTPUT_CLUSTERS_FILENAME))]

    os.remove(os.path.join(results_foldername, OUTPUT_CLUSTERS_FILENAME))
    os.remove(os.path.join(results_foldername, OUTPUT_TREE_FILENAME))
    os.remove(os.path.join(results_foldername, OUTPUT_FBP_FILENAME))
    os.rmdir(results_foldername)

    return clusters, list(unknown_motif_names)


def cluster_table_generator(ranked_motif_names, motif_clusters):
    for idx, cluster in enumerate(motif_clusters):
        ranks = (ranked_motif_names.index(motif) + 1 for motif in cluster)

        for rank, motif in sorted(zip(ranks, cluster), key=operator.itemgetter(0)):
            yield idx + 1, rank, motif


def get_color_mapping(ranked_motif_names, motif_clusters):
    if not motif_clusters:
        return dict()

    if NR_OF_COLORS < len(motif_clusters):
        logging.warning("STAMP: there are more clusters than color codes.")

    clusters = []

    for cluster in motif_clusters:
        clusters.append((min(ranked_motif_names.index(motif) + 1 for motif in cluster), cluster))

    motif_name2color = dict()

    for idx, cluster in enumerate(map(operator.itemgetter(1), sorted(clusters, key=operator.itemgetter(0)))):
        color = get_color_string(idx)

        for motif in cluster:
            motif_name2color[motif] = color

    return motif_name2color
