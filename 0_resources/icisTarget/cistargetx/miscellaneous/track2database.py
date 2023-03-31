import getopt
import gzip
import numpy
import operator
import os
import random
import sys
from functools import partial

from bx.bbi.bigwig_file import BigWigFile

from cistargetx.common.rankingsdatabase import RankingsDatabase
from cistargetx.common.utils import derive_dtype
from cistargetx.conversion.featurelist import Algorithm, Feature, FeatureList, \
    BedFileType, BigBedFileType, \
    BedGraphFileType, WigFileType, BigWigFileType, \
    BroadPeakFileType, BroadPeakBigBedFileType, NarrowPeakFileType, NarrowPeakBigBedFileType, \
    SgrFileType
from cistargetx.conversion.tools import REGIONS_BED_FILENAME


NEGATIVE_INFINITY = float('-inf')


class ScoringAlgorithm(object):
    AVERAGE = "avg"
    MAXIMUM = "max"


def write_database(database_filename, region_ids, feature_ids, feature2regions2r,
                   reduce_multiple_regions_to_one_gene=False):
    if reduce_multiple_regions_to_one_gene:
        # Extract the gene names from the region names by removing "#<number>" at the end of the region name.
        gene_or_region_ids = sorted(set([region_id.rsplit('#', 1)[0] for region_id in region_ids]))
    else:
        gene_or_region_ids = sorted(region_ids)

    gene_or_region_id_count = len(gene_or_region_ids)

    feature_ids = sorted(feature_ids)

    dtype = derive_dtype(gene_or_region_id_count)
    # Fill rankings array with zeros, but do not use numpy.zeros as it does not preallocate the whole array.
    rankings = numpy.full(shape=(gene_or_region_id_count, len(feature_ids)), fill_value=0, dtype=dtype)
    col_idx = 0

    for feature_id in feature_ids:
        rankings[:, col_idx] = [feature2regions2r[feature_id][gene_or_region_id] for gene_or_region_id in gene_or_region_ids]
        col_idx += 1

    RankingsDatabase.create(os.path.basename(database_filename),
                            feature_ids,
                            gene_or_region_ids,
                            rankings).write(database_filename)


def write_txt(database_filename, region_ids, feature_ids, feature2regions2r, reduce_multiple_regions_to_one_gene=False):
    if reduce_multiple_regions_to_one_gene:
        # Extract the gene names from the region names by removing "#<number>" at the end of the region name.
        gene_or_region_ids = sorted(set([region_id.rsplit('#', 1)[0] for region_id in region_ids]))
    else:
        gene_or_region_ids = sorted(region_ids)

    feature_ids = sorted(feature_ids)

    with open(database_filename, 'w') if database_filename != "-" else sys.stdout as output_fh:
        if len(feature_ids) > 1:
            print >> output_fh, ".\t{0:s}".format("\t".join(feature_ids))

        for gene_or_region_id in gene_or_region_ids:
            print >> output_fh, "{0:s}\t{1:s}".format(gene_or_region_id, "\t".join(
                map(str, (feature2regions2r[feature_id][gene_or_region_id] + 1 for feature_id in feature_ids))))


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, "Usage: {0:s} [-r <bed_filename>] [<optional parameters>] <filename1> ... <filenameN> <database_filename>\n".format(
        command_name)
    print >> output_fh, "Optional parameters:"
    print >> output_fh, " -d                            : Dump scores to stdout."
    print >> output_fh, " -f <format>                   : Format of the input files (default: 'auto'):"
    print >> output_fh, "                                   - auto (autodetect format on file extension)"
    print >> output_fh, "                                   - bed, bigbed"
    print >> output_fh, "                                   - bedgraph, wig, bigwig"
    print >> output_fh, "                                   - bwaob (output of bigWigAverageOverBed)"
    print >> output_fh, "                                   - broadpeak, broadpeak_bigbed, narrowpeak, narrowpeak_bigbed"
    print >> output_fh, "                                   - sgr"
    print >> output_fh, " -r <bed_filename>             : BED file that defines the regions to be scored (default: 'dm3-regions')."
    print >> output_fh, " -g                            : Reduce scores of multiple regions around a gene to one score per gene."
    print >> output_fh, " -a avg|max                    : Type of scoring algorithm to use: avg or max (default: 'max')."
    print >> output_fh, " -A [true|yes|1]|[false|no|0]  : Rescale score by the overlap fraction of the region with the feature (default: 'false')."
    print >> output_fh, " -o db|txt                     : Type of output: db or txt (default db)."
    print >> output_fh, " -s <parts_of_feature_name>    : Remove parts of feature name (remove_this_part|remove_this_other_part)."


def feature_iter(region, cont_features):
    for overlap in cont_features.find_overlap_with_feature(region):
        yield overlap


def bigwig_feature_iter(region, bigwig):
    elements = bigwig.get(region.chromosome, region.interval[0], region.interval[1])

    for start, stop, score in (elements if elements else []):
        yield Feature(region.chromosome, int(start), int(stop), score=score)


def score_avg(region, features, rescale_score_by_overlap_fraction):
    sum_overlap = 0.0
    total_n_overlap = 0

    for overlap in features:
        n_overlap = region.get_overlap_in_bp_with(overlap)
        total_n_overlap += n_overlap

        if rescale_score_by_overlap_fraction:
            sum_overlap += overlap.score * (float(n_overlap) / len(overlap))
        else:
            sum_overlap += overlap.score

    return sum_overlap / total_n_overlap if total_n_overlap > 0 else NEGATIVE_INFINITY


def score_max(region, features, rescale_score_by_overlap_fraction):
    scores = []

    for overlap in features:
        n_overlap = region.get_overlap_in_bp_with(overlap)

        if rescale_score_by_overlap_fraction:
            scores.append(overlap.score * (float(n_overlap) / len(overlap)))
        else:
            scores.append(overlap.score)

    return max(scores) if scores else NEGATIVE_INFINITY


def score_regions(regions, scoring_algorithm, features4regions_fnc, rescale_score_by_overlap_fraction,
                  reduce_multiple_regions_to_one_gene, dump2stdout):
    region_or_gene_to_score = dict()

    # Caveat: for this algorithm to work, the score for all inquired regions in
    # the high-throughput experiment must be provided. Even if that score is 0.
    if scoring_algorithm == ScoringAlgorithm.AVERAGE:
        for region in regions:
            if reduce_multiple_regions_to_one_gene:
                # Extract the gene name from the region name by removing "#<number>" at the end of the region name.
                gene_name = region.name.rsplit('#', 1)[0]

                score_avg_value = score_avg(region,
                                            features4regions_fnc(region),
                                            rescale_score_by_overlap_fraction)

                # Update score for gene name when the current value is lower than the old one.
                region_or_gene_to_score[gene_name] = max(region_or_gene_to_score.get(gene_name, score_avg_value),
                                                         score_avg_value)

                if dump2stdout:
                    print >> sys.stdout, gene_name + "\t" + str(region_or_gene_to_score[gene_name])
            else:
                region_or_gene_to_score[region.name] = score_avg(region,
                                                                 features4regions_fnc(region),
                                                                 rescale_score_by_overlap_fraction)

                if dump2stdout:
                    print >> sys.stdout, region.name + "\t" + str(region_or_gene_to_score[region.name])
    else:
        for region in regions:
            if reduce_multiple_regions_to_one_gene:
                # Extract the gene name from the region name by removing "#<number>" at the end of the region name.
                gene_name = region.name.rsplit('#', 1)[0]

                score_max_value = score_max(region,
                                            features4regions_fnc(region),
                                            rescale_score_by_overlap_fraction)

                # Update score for gene name when the current value is lower than the old one.
                region_or_gene_to_score[gene_name] = max(region_or_gene_to_score.get(gene_name, score_max_value),
                                                         score_max_value)

                if dump2stdout:
                    print >> sys.stdout, gene_name + "\t" + str(region_or_gene_to_score[gene_name])
            else:
                region_or_gene_to_score[region.name] = score_max(region,
                                                                 features4regions_fnc(region),
                                                                 rescale_score_by_overlap_fraction)

                if dump2stdout:
                    print >> sys.stdout, region.name + "\t" + str(region_or_gene_to_score[region.name])

    return region_or_gene_to_score


def get_score_regions_bigwig_average_over_bed(bigwig_average_over_bed_tsv_filename, regions, scoring_algorithm,
                                              dump2stdout=False, reduce_multiple_regions_to_one_gene=False):
    """
    Read TSV output file of bigWigAverageOverBed generated from bigWig file
    and same region BED file used as -r parameter for this script.

    The TSV output file of bigWigAverageOverBed has the following columns:
        - name:     name field from bed, which should be unique
        - size:     size of bed (sum of exon sizes
        - covered:  # bases within exons covered by bigWig
        - sum:      sum of values over all bases covered
        - mean0:    average over bases with non-covered bases counting as zeroes
        - mean:     average over just covered bases
        - min:      minimum coverage value (if '-minMax' option was specified)
        - max:      maximum coverage value (if '-minMax' option was specified)
    """

    region_or_gene_to_score = dict()

    if bigwig_average_over_bed_tsv_filename.lower().endswith('.gz'):
        open_fnc = gzip.open
        open_mode = 'rb'
    else:
        open_fnc = open
        open_mode = 'r'

    # Set starts scores for each region to minus infinity as not all chromosomes
    # for all regions are necessary available in the original bigWig file.
    for region in regions:
        if reduce_multiple_regions_to_one_gene:
            # Extract the gene name from the region name by removing "#<number>" at the end of the region name.
            gene_name = region.name.rsplit('#', 1)[0]

            region_or_gene_to_score[gene_name] = NEGATIVE_INFINITY
        else:
            region_or_gene_to_score[region.name] = NEGATIVE_INFINITY

    with open_fnc(bigwig_average_over_bed_tsv_filename, open_mode) as input_fh:
        for line in input_fh:
            if line.startswith('#'):
                continue

            columns = line.rstrip().split('\t')

            if len(columns) == 6:
                if scoring_algorithm == ScoringAlgorithm.MAXIMUM:
                    print >> sys.stderr, "bigWigAverageOverBed TSV output file '{0:s}' should have 8 columns (run with '-minMax') when using '{1:s}' as scoring algorithm.".format(bigwig_average_over_bed_tsv_filename, ScoringAlgorithm.MAXIMUM)
                    sys.exit(1)

                region_name = columns[0]
                region_mean0_coverage = float(columns[4])

                score_avg_value = region_mean0_coverage if region_mean0_coverage != 0 else NEGATIVE_INFINITY

                if reduce_multiple_regions_to_one_gene:
                    # Extract the gene name from the region name by removing "#<number>" at the end of the region name.
                    gene_name = region_name.rsplit('#', 1)[0]

                    # Update score for gene name when the current value is lower than the old one.
                    region_or_gene_to_score[gene_name] = max(region_or_gene_to_score.get(gene_name, score_avg_value),
                                                             score_avg_value)

                    if dump2stdout:
                        print >> sys.stdout, gene_name + "\t" + str(region_or_gene_to_score[gene_name])
                else:
                    region_or_gene_to_score[region_name] = score_avg_value

                    if dump2stdout:
                        print >> sys.stdout, region_name + "\t" + str(region_or_gene_to_score[region_name])
            elif len(columns) == 8:
                region_name = columns[0]

                if scoring_algorithm == ScoringAlgorithm.AVERAGE:
                    region_mean0_coverage = float(columns[4])

                    score_avg_or_max_value = region_mean0_coverage if region_mean0_coverage != 0 else NEGATIVE_INFINITY
                else:
                    region_max_coverage = float(columns[7])

                    score_avg_or_max_value = region_max_coverage if region_max_coverage != 0 else NEGATIVE_INFINITY

                if reduce_multiple_regions_to_one_gene:
                    # Extract the gene name from the region name by removing "#<number>" at the end of the region name.
                    gene_name = region_name.rsplit('#', 1)[0]

                    # Update score for gene name when the current value is lower than the old one.
                    region_or_gene_to_score[gene_name] = max(region_or_gene_to_score.get(gene_name,
                                                                                         score_avg_or_max_value),
                                                             score_avg_or_max_value)

                    if dump2stdout:
                        print >> sys.stdout, gene_name + "\t" + str(region_or_gene_to_score[gene_name])
                else:
                    region_or_gene_to_score[region_name] = score_avg_or_max_value

                    if dump2stdout:
                        print >> sys.stdout, region_name + "\t" + str(region_or_gene_to_score[region_name])
            else:
                continue

    return region_or_gene_to_score


def compare(i1, i2):
    if i1 == NEGATIVE_INFINITY and i2 == NEGATIVE_INFINITY:
        return random.randint(-1, 1)
    if i1 < i2:
        return -1
    if i1 > i2:
        return 1

    return 0


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "dgha:A:f:o:r:s:")
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        display_usage(sys.stderr)
        sys.exit(2)

    if len(args) < 2:
        print >> sys.stderr, "Wrong number of input arguments."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    autodetect_input_format_filetype = True
    regions_bed_filename = REGIONS_BED_FILENAME
    reduce_multiple_regions_to_one_gene = False
    dump2stdout = False
    scoring_algorithm = ScoringAlgorithm.MAXIMUM
    rescale_score_by_overlap_fraction = False
    output_fnc = write_database
    remove_parts_of_feature_id = None

    for o, a in opts:
        if o == "-f":
            a_lower = a.lower()

            if a_lower == 'auto':
                autodetect_input_format_filetype = True
            else:
                autodetect_input_format_filetype = False

                if a_lower == 'bed':
                    input_format_filetype = BedFileType()
                elif a_lower == 'bigbed':
                    input_format_filetype = BigBedFileType()

                elif a_lower == 'bedgraph':
                    input_format_filetype = BedGraphFileType()
                elif a_lower == 'wig':
                    input_format_filetype = WigFileType()
                elif a_lower == 'bigwig_old':
                    input_format_filetype = None
                    input_format_filetype_id = 'bigwig_old'
                elif a_lower == 'bigwig':
                    input_format_filetype = BigWigFileType()

                elif a_lower == 'bwaob':
                    input_format_filetype = None
                    input_format_filetype_id = 'bwaob'

                elif a_lower == 'broadpeak':
                    input_format_filetype = BroadPeakFileType()
                elif a_lower == 'broadpeak_bigbed':
                    input_format_filetype = BroadPeakBigBedFileType()
                elif a_lower == 'narrowpeak':
                    input_format_filetype = NarrowPeakFileType()
                elif a_lower == 'narrowpeak_bigbed':
                    input_format_filetype = NarrowPeakBigBedFileType()

                elif a_lower == 'sgr':
                    input_format_filetype = SgrFileType()

                else:
                    print >> sys.stderr, "'{0:s}' is not a valid input format.".format(a)
                    display_usage(sys.stderr)
                    sys.exit(1)
        elif o == "-o":
            a_lower = a.lower()

            if a_lower == 'db':
                output_fnc = write_database
            elif a_lower == 'txt':
                output_fnc = write_txt
            else:
                print >> sys.stderr, "'{0:s}' is not a valid input format.".format(a)
                display_usage(sys.stderr)
                sys.exit(1)
        elif o == "-r":
            regions_bed_filename = a
        elif o == "-g":
            reduce_multiple_regions_to_one_gene = True
        elif o == "-d":
            dump2stdout = True
        elif o == "-h":
            display_usage()
            sys.exit()
        elif o == "-a":
            a_lower = a.lower()

            if a_lower == 'avg':
                scoring_algorithm = ScoringAlgorithm.AVERAGE
            elif a_lower == 'max':
                scoring_algorithm = ScoringAlgorithm.MAXIMUM
            else:
                print >> sys.stderr, "'{0:s}' is not a valid type of algorithm.".format(a)
                display_usage(sys.stderr)
                sys.exit(1)
        elif o == "-A":
            a_lower = a.lower()

            if a_lower == '1' or a_lower == "true" or a_lower == "yes":
                rescale_score_by_overlap_fraction = True
            elif a_lower == '0' or a_lower == "false" or a_lower == "no":
                rescale_score_by_overlap_fraction = False
            else:
                print >> sys.stderr, "'{0:s}' is not a valid value for rescale score by overlap fraction.".format(a)
                display_usage(sys.stderr)
                sys.exit(1)
        elif o == "-s":
            remove_parts_of_feature_id = a
        else:
            assert False, "Unhandled option"

    input_filenames = args[:-1]
    database_filename = args[-1]

    if os.path.exists(database_filename):
        print >> sys.stderr, "'{0:s}' already exists.".format(database_filename)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    regions = FeatureList.from_file(regions_bed_filename, fast_algorithms=Algorithm(False, False, True))
    region_ids = [feature.name for feature in regions]
    feature_ids = []
    feature2regions2r = dict()

    for filename in input_filenames:
        filename_basename = os.path.basename(filename)
        filename_basename_lowercase = filename_basename.lower()

        feature_id_lowercase, extension = (os.path.splitext(filename_basename_lowercase[:-3])
                                           if filename_basename_lowercase.endswith('.gz')
                                           else os.path.splitext(filename_basename_lowercase))

        feature_id_length = len(feature_id_lowercase)

        if extension == '.bed':
            autodetected_input_format_filetype = BedFileType()

            if feature_id_lowercase.endswith('.broadpeak'):
                autodetected_input_format_filetype = BroadPeakFileType()
                feature_id_length -= len('.broadpeak')
            elif feature_id_lowercase.endswith('.narrowpeak'):
                autodetected_input_format_filetype = NarrowPeakFileType()
                feature_id_length -= len('.narrowpeak')

        elif extension == '.bigbed' or extension == '.bb':
            autodetected_input_format_filetype = BigBedFileType()

            if feature_id_lowercase.endswith('.broadpeak'):
                autodetected_input_format_filetype = BroadPeakBigBedFileType()
                feature_id_length -= len('.broadpeak')
            elif feature_id_lowercase.endswith('.narrowpeak'):
                autodetected_input_format_filetype = NarrowPeakBigBedFileType()
                feature_id_length -= len('.narrowpeak')

        elif extension == '.bedgraph':
            autodetected_input_format_filetype = BedGraphFileType()
        elif extension == '.wig':
            autodetected_input_format_filetype = WigFileType()
        elif extension == '.bigwig' or extension == '.bw':
            autodetected_input_format_filetype = BigWigFileType()

        elif extension == '.bwaob':
            autodetected_input_format_filetype = None
            autodetected_input_format_filetype_id = 'bwaob'

        elif extension == '.broadpeak':
            autodetected_input_format_filetype = BroadPeakFileType()
        elif extension == '.narrowpeak':
            autodetected_input_format_filetype = NarrowPeakFileType()

        elif extension == '.sgr':
            autodetected_input_format_filetype = SgrFileType()

        elif autodetect_input_format_filetype:
            print >> sys.stderr, "File format could not automatically be detected for '{0:s}'.".format(filename)
            display_usage(sys.stderr)
            sys.exit(1)

        if autodetect_input_format_filetype:
            input_format_filetype = autodetected_input_format_filetype

            if input_format_filetype is None:
                input_format_filetype_id = autodetected_input_format_filetype_id

        # Cleaned feature ID (without filename extensions).
        feature_id = filename_basename[:feature_id_length]

        if remove_parts_of_feature_id:
            for remove_part_of_feature_id in remove_parts_of_feature_id.split('|'):
                # Remove additional strings from feature ID.
                feature_id = feature_id.replace(remove_part_of_feature_id, '')

        if input_format_filetype:
            input_format_filetype_id = input_format_filetype.filetype_id

        print >> sys.stderr, "Analyzing feature '{0:s}' (format: {1:s}) ... ".format(feature_id,
                                                                                     input_format_filetype_id)

        if dump2stdout:
            print >> sys.stdout, "#" + feature_id

        try:
            if not input_format_filetype:
                if input_format_filetype_id == 'bwaob':
                    region2score = get_score_regions_bigwig_average_over_bed(
                        bigwig_average_over_bed_tsv_filename=filename,
                        regions=regions,
                        scoring_algorithm=scoring_algorithm,
                        dump2stdout=dump2stdout,
                        reduce_multiple_regions_to_one_gene=reduce_multiple_regions_to_one_gene
                    )
                else:
                    with open(filename, 'rb') as input_fh:
                        region2score = score_regions(
                            regions=regions,
                            scoring_algorithm=scoring_algorithm,
                            features4regions_fnc=partial(bigwig_feature_iter, bigwig=BigWigFile(input_fh)),
                            rescale_score_by_overlap_fraction=rescale_score_by_overlap_fraction,
                            reduce_multiple_regions_to_one_gene=reduce_multiple_regions_to_one_gene,
                            dump2stdout=dump2stdout
                        )
            else:
                cont_features = FeatureList.from_file(filename,
                                                      filetype=input_format_filetype,
                                                      fast_algorithms=Algorithm(True, False, False))

                region2score = score_regions(regions=regions,
                                             scoring_algorithm=scoring_algorithm,
                                             features4regions_fnc=partial(feature_iter, cont_features=cont_features),
                                             rescale_score_by_overlap_fraction=rescale_score_by_overlap_fraction,
                                             reduce_multiple_regions_to_one_gene=reduce_multiple_regions_to_one_gene,
                                             dump2stdout=dump2stdout)

                del cont_features

            feature_ids.append(feature_id)
            feature2regions2r[feature_id] = dict((region_id, idx)
                                                 for idx, region_id in enumerate(map(operator.itemgetter(0),
                                                                                     sorted(region2score.iteritems(),
                                                                                            cmp=compare,
                                                                                            key=operator.itemgetter(1),
                                                                                            reverse=True))))
        except Exception as e:
            print >> sys.stderr, "Error while scoring {0:s}: {1:s}".format(feature_id, str(e))

    output_fnc(database_filename, region_ids, feature_ids, feature2regions2r, reduce_multiple_regions_to_one_gene)


if __name__ == "__main__":
    main()
