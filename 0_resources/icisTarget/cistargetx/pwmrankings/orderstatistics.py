import operator
import os
import sys

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

ORDERSTATISTICS_CMD = os.path.join(os.path.dirname(__file__), 'OrderStatistics.jar')


def combined_rankings_iterator(cfg, ranking_filenames, gene2regions=None):
    """ Iterator over (region, rank) tuples (or (gene, rank) tuples if gene2regions is specified). """

    # Running orderstatistics.
    exec_command = ['java', '-jar', ORDERSTATISTICS_CMD]
    exec_command += ranking_filenames

    try:
        stdout, stderr = subprocess.Popen(exec_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    except OSError, msg:
        print >> sys.stderr, "orderstatistics: error during execution of '" + ' '.join(exec_command) + "': " + str(msg)
        sys.exit(1)

    if stderr and stderr != "":
        print >> sys.stderr, "orderstatistics: error during execution of '" + ' '.join(exec_command) + "': " + str(
            stderr)
        sys.exit(1)

    # Parsing output of orderstatistics.
    region2rank = dict()

    for line in stdout.split('\n'):
        columns = line.split('\t')

        # Skip lines that don't contain a ranking.
        if len(columns) != 2:
            continue

        region = columns[0];
        score = float(columns[1])

        # CAVE: orderstatistics outputs only one ranking for each region.
        region2rank[region] = score

    # Save os file in debug mode.
    if cfg.has_option('log', 'debug_level') and cfg.getint('log', 'debug_level') > 0:
        output_folder, filename = os.path.split(ranking_filenames[0])
        pwm = '.'.join(filename.split('.')[0:-2])
        orderstat_filename = os.path.join(output_folder, pwm + ".os")

        with open(orderstat_filename, 'w') as output_fh:
            for region, rank in region2rank.iteritems():
                # CAVE: Rank must be saved using g format specifier.
                print >> output_fh, "{0:s}\t{1:.16g}".format(region, rank)

    # Mapping of regions to genes if necessary.
    if gene2regions:
        gene2rank = dict()

        for gene in gene2regions.keys():
            gene2rank[gene] = min(region2rank[region] for region in gene2regions[gene])

        ranked_gene_or_region_ids = [gene for gene, score in sorted(gene2rank.items(),
                                                                    reverse=False,
                                                                    key=operator.itemgetter(1))]
    else:
        ranked_gene_or_region_ids = [region for region, score in sorted(region2rank.items(),
                                                                        reverse=False,
                                                                        key=operator.itemgetter(1))]
    for idx, gene_or_region_id in enumerate(ranked_gene_or_region_ids):
        yield gene_or_region_id, idx + 1
