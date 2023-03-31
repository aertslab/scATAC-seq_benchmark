import gzip
import re

import numpy
import pyBigWig

# http://biostar.stackexchange.com/questions/99/fast-interval-intersection-methodologies
from bx.intervals.intersection import Intersecter, Interval


class FileType(object):
    def __init__(self, filetype_id):
        self.filetype_id = filetype_id

    def feature_iterator(self, filename, id_transform=lambda x: x):
        is_bigwig_or_bigbed = False

        if (self.filetype_id == 'bed'
                or self.filetype_id == 'bedgraph'
                or self.filetype_id == 'wig'
                or self.filetype_id == 'broadpeak'
                or self.filetype_id == 'narrowpeak'
                or self.filetype_id == 'sgr'):
            is_bigwig_or_bigbed = False
        elif (self.filetype_id == 'bigbed'
                or self.filetype_id == 'bigwig'
                or self.filetype_id == 'broadpeak_bigbed'
                or self.filetype_id == 'narrowpeak_big_bed'):
            is_bigwig_or_bigbed = True
        else:
            raise('Unsupported filetype ID "{0:s}"'.format(self.filetype_id))

        if is_bigwig_or_bigbed:
            for feature in self.feature_iterator_for_file(filename, id_transform):
                yield feature
        else:
            if filename.lower().endswith('.gz'):
                open_fnc = gzip.open
                open_mode = 'rb'
            else:
                open_fnc = open
                open_mode = 'r'

            with open_fnc(filename, open_mode) as input_fh:
                for feature in self.feature_iterator_for_stream(input_fh, id_transform):
                    yield feature

    def feature_iterator_for_stream(self, input_fh, id_transform):
        pass

    def feature_iterator_for_file(self, filename, id_transform):
        pass


class BedFileType(FileType):
    """
    Read a BED file and yield Feature objects.

    BED format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format1
    """

    def __init__(self):
        FileType.__init__(self, "bed")

    def feature_iterator_for_stream(self, input_fh, id_transform=lambda x: x):
        for line in input_fh:
            yield Feature.from_string(line, id_transform)


class BigBedFileType(FileType):
    """
    Read a bigBed file and yield Feature objects.

    bigBed format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format1.5
    """

    def __init__(self):
        FileType.__init__(self, "bigbed")

    def feature_iterator_for_file(self, filename, id_transform=lambda x: x):
        def convert_bigbed_entry_to_feature(chromosome, bigbed_entry):
            # Each bigBed entry is a tuple of 3 elements:
            #   - start
            #   - end
            #   - string with the rest of the fields separated by TABs
            start, stop = bigbed_entry[0:2]
            name, score, strand = bigbed_entry[2].split('\t')[0:3]

            # Use score value as score.
            return Feature(chromosome, start, stop, name, score=float(score), strand=strand)

        with pyBigWig.open(filename, 'rb') as bigbed:
            # Only proceed if we have a bigBed file.
            if bigbed.isBigBed():
                chromosomes = bigbed.chroms()

                # Get data for all chromosomes.
                for chromosome in chromosomes:
                    # Get all bigBed entries for the current chromosome.
                    bigbed_entries = bigbed.entries(chromosome, 0, chromosomes[chromosome])

                    if bigbed_entries:
                        for bigbed_entry in bigbed_entries:
                            yield convert_bigbed_entry_to_feature(chromosome=chromosome,
                                                                  bigbed_entry=bigbed_entry)


class BedGraphFileType(FileType):
    """
    Read a bedGraph file and yield Feature objects.

    bedGraph format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format1.8
    """

    def __init__(self):
        FileType.__init__(self, "bedgraph")

    def feature_iterator_for_stream(self, input_fh, id_transform=lambda x: x):
        for line in input_fh:
            if line.startswith('track') or line.startswith('browser') or line.startswith('#'):
                continue

            columns = re.split('[\t ]+', line.rstrip())
            name = id_transform(Feature.DUMMY_NAME)
            score = float(re.sub(',', '.', columns[3]))

            yield Feature(columns[0], int(columns[1]), int(columns[2]), name, score)


class WigFileType(FileType):
    """
    Read a WIG file and yield Feature objects.

    WIG format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format6
    """

    def __init__(self):
        FileType.__init__(self, "wig")

    def feature_iterator_for_stream(self, input_fh, id_transform=lambda x: x):
        fixed_step_pattern = re.compile(r"^fixedStep\s+chrom=(\w+)\s+start=(\d+)\s+step=(\d+)(\s+span=(\d+))?$")
        variable_step_pattern = re.compile(r"^variableStep\s+chrom=(\w+)(\s+span=(\d+))?$")
        in_variable_region = False

        for line in input_fh:
            if line.startswith("track"):
                continue

            match = fixed_step_pattern.match(line)

            if match:
                in_variable_region = False
                chromosome = match.group(1)
                start = int(match.group(2))
                step = int(match.group(3))
                span = int(match.group(5)) if match.group(5) else 1
                count = 0
            elif variable_step_pattern.match(line):
                match = variable_step_pattern.match(line)
                in_variable_region = True
                chromosome = match.group(1)
                span = int(match.group(3)) if match.group(3) else 1
            elif in_variable_region:
                start_position, score = line.rstrip().split()[0:2]
                start_position = int(start_position)
                score = float(score)
                yield Feature(chromosome, start_position, start_position + span, name=id_transform(Feature.DUMMY_NAME),
                              score=score)
            else:
                # Fixed step region.
                start_position = start + (count * step)
                score = float(line.rstrip())
                yield Feature(chromosome, start_position, start_position + span, name=id_transform(Feature.DUMMY_NAME),
                              score=score)
                count += 1


class BigWigFileType(FileType):
    """
    Read a bigWig file and yield Feature objects.

    bigWig format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format6.1
    """

    def __init__(self):
        FileType.__init__(self, "bigwig")

    def feature_iterator_for_file(self, filename, id_transform=lambda x: x):
        with pyBigWig.open(filename, 'rb') as bigwig:
            # Only proceed if we have a bigWig file.
            if bigwig.isBigWig():
                chromosomes = bigwig.chroms()

                for chromosome in chromosomes:
                    chrom_size = chromosomes[chromosome]
                    step_size = 10000

                    # Retrieve the values from the bigWig file in chucks, so we
                    # don't need to retrieve all values per chromosome at once.
                    for start_chunk_pos in range(0, chrom_size, step_size):
                        # Calculate end chunk position and check if if is not larger than the chromosome size.
                        end_chunk_pos = (start_chunk_pos + step_size
                                         if start_chunk_pos + step_size < chrom_size
                                         else chrom_size)

                        # Get all bigWig values for the current chunk.
                        bigwig_values_chunk = numpy.array(bigwig.values(chromosome, start_chunk_pos, end_chunk_pos))

                        # Get all positions which have a value (not "nan").
                        bigwig_values_chunk_idx_not_nan = numpy.where(numpy.logical_not(numpy.isnan(bigwig_values_chunk)))

                        # Loop over all positions in the chunk which have a value (not "nan").
                        for score_idx, score in zip(bigwig_values_chunk_idx_not_nan,
                                                    bigwig_values_chunk[bigwig_values_chunk_idx_not_nan]):
                            # Calculate the start and end positions.
                            start = start_chunk_pos + int(score_idx)
                            stop = start + 1

                            yield Feature(chromosome, start, stop, score=float(score))


class BroadPeakFileType(FileType):
    """
    Read a broadPeak file and yield Feature objects.

    ENCODE broadPeak: Broad Peaks (or Regions) format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format13
    """

    def __init__(self):
        FileType.__init__(self, "broadpeak")

    def feature_iterator_for_stream(self, input_fh, id_transform=lambda x: x):
        for line in input_fh:
            if line.startswith('#'):
                continue

            chromosome, start, stop, name, score, strand, signal_value = line.rstrip().split('\t')[0:7]

            # Use signal value as score.
            yield Feature(chromosome, int(start), int(stop), name, score=float(signal_value), strand=strand)


class BroadPeakBigBedFileType(FileType):
    """
    Read a broadPeak bigBed file and yield Feature objects.

    ENCODE broadPeak: Broad Peaks (or Regions) format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format13
    bigBed format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format1.5
    """

    def __init__(self):
        FileType.__init__(self, "broadpeak_bigbed")

    def feature_iterator_for_file(self, filename, id_transform=lambda x: x):
        def convert_bigbed_entry_to_feature(chromosome, bigbed_entry):
            # Each bigBed entry is a tuple of 3 elements:
            #   - start
            #   - end
            #   - string with the rest of the fields separated by TABs
            start, stop = bigbed_entry[0:2]
            name, score, strand, signal_value = bigbed_entry[2].split('\t')[0:4]

            # Use signal value as score.
            return Feature(chromosome, start, stop, name, score=float(signal_value), strand=strand)

        with pyBigWig.open(filename, 'rb') as bigbed:
            # Only proceed if we have a bigBed file.
            if bigbed.isBigBed():
                chromosomes = bigbed.chroms()

                # Get data for all chromosomes.
                for chromosome in chromosomes:
                    # Get all bigBed entries for the current chromosome.
                    bigbed_entries = bigbed.entries(chromosome, 0, chromosomes[chromosome])

                    if bigbed_entries:
                        for bigbed_entry in bigbed_entries:
                            yield convert_bigbed_entry_to_feature(chromosome=chromosome,
                                                                  bigbed_entry=bigbed_entry)


class NarrowPeakFileType(FileType):
    """
    Read a narrowPeak file and yield Feature objects.

    ENCODE narrowPeak: Narrow (or Point-Source) Peaks format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format12
    """

    def __init__(self):
        FileType.__init__(self, "narrowpeak")

    def feature_iterator_for_stream(self, input_fh, id_transform=lambda x: x):
        for line in input_fh:
            if line.startswith('#'):
                continue

            chromosome, start, stop, name, score, strand, signal_value = line.rstrip().split('\t')[0:7]

            # Use signal value as score.
            yield Feature(chromosome, int(start), int(stop), name, score=float(signal_value), strand=strand)


class NarrowPeakBigBedFileType(FileType):
    """
    Read a narrowPeak bigBed file and yield Feature objects.

    ENCODE narrowPeak: Narrow (or Point-Source) Peaks format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format12
    bigBed format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#format1.5
    """

    def __init__(self):
        FileType.__init__(self, "narrowpeak_bigbed")

    def feature_iterator_for_file(self, filename, id_transform=lambda x: x):
        def convert_bigbed_entry_to_feature(chromosome, bigbed_entry):
            # Each bigBed entry is a tuple of 3 elements:
            #   - start
            #   - end
            #   - string with the rest of the fields separated by TABs
            start, stop = bigbed_entry[0:2]
            name, score, strand, signal_value = bigbed_entry[2].split('\t')[0:4]

            # Use signal value as score.
            return Feature(chromosome, start, stop, name, score=float(signal_value), strand=strand)

        with pyBigWig.open(filename, 'rb') as bigbed:
            # Only proceed if we have a bigBed file.
            if bigbed.isBigBed():
                chromosomes = bigbed.chroms()

                # Get data for all chromosomes.
                for chromosome in chromosomes:
                    # Get all bigBed entries for the current chromosome.
                    bigbed_entries = bigbed.entries(chromosome, 0, chromosomes[chromosome])

                    if bigbed_entries:
                        for bigbed_entry in bigbed_entries:
                            yield convert_bigbed_entry_to_feature(chromosome=chromosome,
                                                                  bigbed_entry=bigbed_entry)


class SgrFileType(FileType):
    """
    Read a SGR file and yield Feature objects.

    SGR format:

        A SGR file has three columns of data:
          - chromosome
          - position
          - score

        https://code.google.com/p/biotoolbox/wiki/DataFileFormats
    """

    def __init__(self):
        FileType.__init__(self, "sgr")

    def feature_iterator_for_stream(self, input_fh, id_transform=lambda x: x):
        for line in input_fh:
            chromosome, position, score = line.rstrip().split()[0:3]
            position = int(position)
            score = float(score)
            yield Feature(chromosome, position, position + 1, name=id_transform(Feature.DUMMY_NAME), score=score)


# Better implemented via Strategy Pattern on FeatureList, but in the interest of time
# faster implementation is opted for.
# TODO: Provide implementations for all operations in the absence of fast data structures.
class Algorithm(object):
    def __init__(self, fast_overlap_select=True, fast_name_lookup=True, fast_random_access=True):
        self.fast_overlap_select = fast_overlap_select
        self.fast_name_lookup = fast_name_lookup
        self.fast_random_access = fast_random_access


class FeatureList(object):
    """ Current implementation has a double memory load, this allows for
        fast overlap select operations.
        IDs of features in a list do not have to be unique.
        CAVE: Doesn't take strandedness of features into account.
    """

    @staticmethod
    def is_bed_file(data):
        for line in re.split('[\n\r]+', data):
            if not line or re.match('^[ \t]+$', line) or line.startswith("track"):
                continue

            return re.match('^\w+[ \t]+[0-9]+[ \t]+[0-9]+([ \t].*)?$', line) is not None

        return False

    # Deprecated: better use FeatureList.from_file(filename, transform)
    @staticmethod
    def from_bed_file(filename, transform=lambda x: x):
        def _feature_iterator(filename):
            with open(filename, 'r') as input_fh:
                for line in input_fh:
                    yield Feature.from_string(line, transform)

        return FeatureList(_feature_iterator(filename))

    @staticmethod
    def from_file(filename, filetype=BedFileType(), transform=lambda x: x, fast_algorithms=Algorithm()):
        return FeatureList(filetype.feature_iterator(filename, transform), fast_algorithms)

    @staticmethod
    def from_stream(input_fh, filetype=BedFileType(), transform=lambda x: x, fast_algorithms=Algorithm()):
        return FeatureList(filetype.feature_iterator_for_stream(input_fh, transform), fast_algorithms)

    @staticmethod
    def from_string(data, transform=lambda x: x):
        def _feature_iterator(data):
            for line in re.split('[\n\r]+', data):
                if not line or re.match('^[ \t]+$', line):
                    continue
                if not line.startswith("track"):
                    yield Feature.from_string(line, transform)

        return FeatureList(_feature_iterator(data))

    def __init__(self, features_iterator, fast_algorithms=Algorithm()):
        self.features = []
        self.chromosome2tree = dict()
        self.name2features = dict()

        for feature in features_iterator:
            name = feature.name
            chromosome = feature.chromosome
            start, end = feature.interval

            if fast_algorithms.fast_random_access:
                # 1. Data structure required for fast random access to features.
                self.features.append(feature)
            if fast_algorithms.fast_name_lookup:
                # 2. Data structure required for fast access on name of feature.
                self.name2features.setdefault(name, []).append(feature)
            if fast_algorithms.fast_overlap_select:
                # 3. Data structure required for fast overlap select.
                if chromosome not in self.chromosome2tree.keys():
                    self.chromosome2tree[chromosome] = Intersecter()

                self.chromosome2tree[chromosome].add_interval(
                    Interval(start, end, {'name': name, 'score': feature.score}))

    def __len__(self):
        return len(self.features)

    def __getitem__(self, param):
        if isinstance(param, str):
            if param not in self.name2features:
                raise KeyError

            features = self.name2features[param]

            if len(features) == 1:
                return features[0]
            else:
                raise KeyError
        elif isinstance(param, int):
            if param < 0 or param >= len(self):
                raise IndexError

            return self.features[param]
        else:
            raise ValueError

    def get(self, param, default=None):
        try:
            return self.__getitem__(param)
        except KeyError:
            return default

    def __iter__(self):
        return self.features.__iter__()

    def __str__(self):
        return "\n".join(map(str, self))

    def find_overlap_with_feature(self, feature, fraction=None):
        if feature.chromosome in self.chromosome2tree.keys():
            start, end = feature.interval
            overlap = self.chromosome2tree[feature.chromosome].find(start, end)

            # In the following code snippets 'overlap_feature' is the feature that
            # corresponds to a conserved region and 'feature' represents a peak.

            # This version of the fraction overlap filter is equivalent to:
            #   bedtools intersect -wa -a <BED file with scored regions> -b <BED file with peaks> -f <overlap_fraction> | cut -f 4 | sort -u
            # (The -f parameter is always relative to the bed file denoted as 'a')
            def fraction_filter1(overlap_feature):
                return (not fraction
                        or (float(overlap_feature.get_overlap_in_bp_with(feature)) / len(overlap_feature)) >= fraction)

            # This version of the fraction overlap filter is equivalent to:
            #   bedtools intersect -wo -e -f 0.4 -F 0.4 -a <BED file with scored regions> -b <BED file with peaks | cut -f 4 | sort -u
            def fraction_filter2(overlap_feature):
                if not fraction:
                    return True

                overlap_in_bp = float(overlap_feature.get_overlap_in_bp_with(feature))

                if len(feature) == 0:
                    overlap_fraction_relative_to_feature = 0.0
                else:
                    overlap_fraction_relative_to_feature = overlap_in_bp / len(feature)

                if len(overlap_feature) == 0:
                    overlap_fraction_relative_to_overlap_feature = 0.0
                else:
                    overlap_fraction_relative_to_overlap_feature = overlap_in_bp / len(overlap_feature)

                return max(overlap_fraction_relative_to_feature,
                           overlap_fraction_relative_to_overlap_feature) >= fraction

            # This version of the fraction overlap filter is equivalent to:
            #   bedtools intersect -wo -a <BED file with peaks> -b <BED file with scored regions> -f <overlap_fraction> | cut -f 8 | sort -u
            #   bedtools intersect -wa -a <BED file with scored regions> -b <BED file with peaks> -f 1.0 | cut -f 4 | sort -u
            # (The -f parameter is always relative to the bed file denoted as 'a')
            def fraction_filter3(overlap_feature):
                return (not fraction
                        or overlap_feature in feature
                        or (float(overlap_feature.get_overlap_in_bp_with(feature)) / len(feature)) >= fraction)

            #            for i in overlap:
            #                print >>sys.stderr, "{0:s}\t{1:d}\t{2:d}\t{3:s}\t{4:s}\t{5:d}\t{6:d}\t{7:s}\t{8:d}\t{9:d}\t{10:f}".format(
            #                    feature.chromosome, feature.interval[0], feature.interval[1], feature.name,
            #                    feature.chromosome, i.start, i.end, i.value, len(feature),
            #                    Feature(feature.chromosome, i.start, i.end, i.value).get_overlap_in_bp_with(feature),
            #                    float(Feature(feature.chromosome, i.start, i.end, i.value).get_overlap_in_bp_with(feature))/len(feature))

            return filter(fraction_filter2,
                          [Feature(feature.chromosome, i.start, i.end, i.value['name'], i.value['score'])
                           for i in overlap])
        else:
            return []

    def find_overlap_with(self, other, fraction=None):
        def _feature_iterator(self, other):
            for feature1 in other:
                for feature2 in self.find_overlap_with_feature(feature1, fraction):
                    yield feature2

        return FeatureList(_feature_iterator(self, other))

    def find_overlap_with_and_return_other(self, other, fraction=None):
        def _feature_iterator(self, other):
            for feature1 in other:
                for feature2 in self.find_overlap_with_feature(feature1, fraction):
                    yield feature1

        return FeatureList(_feature_iterator(self, other))

    def filter_by_name(self, names):
        def _feature_iterator(featureList, names):
            for name in names:
                for feature in featureList.name2features.get(name, []):
                    yield feature

        # Avoid duplicates.
        return FeatureList(_feature_iterator(self, set(names)))

    def write(self, filename, track_name='', description=''):
        with open(filename, 'w') as output_fh:
            if track_name:
                output_fh.write('track name={0:s} description="{1:s}" useScore=1\n'.format(track_name, description))

            output_fh.write(str(self))

    @property
    def feature_ids(self):
        return self.name2features.keys()


class Feature:
    """ A class of genomic features defined by a half-open interval. Locations are 0-based."""

    DUMMY_NAME = '.'

    @staticmethod
    def from_string(line, transform=lambda x: x):
        columns = re.split('[\t ]+', line.rstrip())

        assert len(columns) >= 3, "Invalid BED file supplied: at least three columns are expected. Please, check carefully that the correct input type was selected."
        assert re.match("[0-9]+", columns[1]), "Invalid BED file supplied: second column must contain integers."
        assert re.match("[0-9]+", columns[2]), "Invalid BED file supplied: third column must contain integers."

        name = Feature.DUMMY_NAME if len(columns) == 3 else transform(columns[3])

        try:
            score = float(re.sub(',', '.', columns[4])) if len(columns) >= 5 else None
        except ValueError:
            raise AssertionError("Invalid BED file supplied: fifth column must contain floating point numbers (score).")

        strand = columns[5] if len(columns) >= 6 else None

        assert not strand or strand in (
            '+', '-', '.', '?'), "Invalid BED file supplied: sixth column must contain strand (+/-/?)."

        return Feature(columns[0], int(columns[1]), int(columns[2]), name, score, strand)

    def __init__(self, chromosome, start, end, name=DUMMY_NAME, score=None, strand=None):
        assert chromosome.strip() != ""
        assert end >= start
        assert name.strip() != ""

        self.chromosome = chromosome
        self.interval = (start, end)
        self.name = name
        self.score = score
        self.strand = strand

    def __str__(self):
        r = "{0:s}\t{1:d}\t{2:d}\t{3:s}".format(self.chromosome, self.interval[0], self.interval[1], self.name)

        if self.score and self.strand:
            r += "\t{0:f}\t{1:s}".format(self.score, self.strand)
        elif self.score:
            r += "\t{0:f}".format(self.score)
        elif self.strand:
            r += "\t0.0\t{0:s}".format(self.strand)
        return r

    def __len__(self):
        """ Length of feature in base pairs. """
        return self.interval[1] - self.interval[0]

    def has_overlap_with(self, other):
        return (self.chromosome == other.chromosome
                and self.interval[0] < other.interval[1]
                and self.interval[1] > other.interval[0])

    def __contains__(self, other):
        return (self.chromosome == other.chromosome
                and other.interval[0] >= self.interval[0]
                and other.interval[1] <= self.interval[1])

    def get_overlap_in_bp_with(self, other):
        if not self.has_overlap_with(other):
            return 0

        return min(self.interval[1], other.interval[1]) - max(self.interval[0], other.interval[0])

    @property
    def start(self):
        return self.interval[0]

    @property
    def end(self):
        return self.interval[1]
