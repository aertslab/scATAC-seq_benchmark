import collections
import cStringIO
import gzip
import hashlib
import itertools
import logging
import math
import numpy
import operator
import os
import re
import shutil
import zipfile

import matplotlib


# Anti-Grain Graphics library.
matplotlib.use('agg')

import pylab
import wand.image
import wand.exceptions

from motifclustering import get_color_mapping, cluster, cluster_table_generator
from recoverycurves import RecoveryCurves
from trackclustering import track_cluster
from utils import ContextWrapper, write_full_statistics_as_tsv, write_overview_as_html, write_overview_as_tsv, \
    write_template


INI_FILENAME = 'parameters.ini'
# Must end in ".txt".
QUERY_IDS_FILENAME = 'query.id.txt'
STATISTICS_FILENAME = 'statistics.tbl'
BATCH_RESULTS_FILENAME = 'report.tsv.gz'
MOTIFS_CB_FILENAME = 'motifs.cb'
CLUSTERS_FILENAME = 'clusters.tbl'
HTML_REPORT_FILENAME = 'report.html'
HIST_FILENAME = 'aucdistributions.png'
OVERVIEW_FILENAME = 'rccoverview.png'
ARCHIVE_FILENAME = 'archive.zip'
AVERAGE_RECOVERY_CURVE_FILENAME = 'average.rcc'
AVG2STD_RECOVERY_CURVE_FILENAME = 'avg2std.rcc'

HTML_EXTENSION = 'html'
TSV_EXTENSION = 'tsv'
BED_EXTENSION = 'bed'
RECOVERY_CURVE_EXTENSION = 'rcc'
RANK_TABLE_EXTENSION = 'tbl'
FIGURE_EXTENSION = 'png'
DATA_EXTENSION = 'tsv'

GROUP_COLORS = ['r', 'g', 'b', 'c', 'y', 'm']
GROUP_COLORS_LENGTH = len(GROUP_COLORS)


def write_gene_or_region_ids(gene_or_region_ids, filename):
    with open(filename, 'w') as output_fh:
        for gene_or_region_id in gene_or_region_ids:
            print >> output_fh, gene_or_region_id


def install_output_folder(output_folder, ini_filename):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    logos_folder = os.path.join(output_folder, 'logos')

    if not os.path.exists(logos_folder):
        os.mkdir(logos_folder)

    if ini_filename:
        shutil.copy(ini_filename, os.path.join(output_folder, INI_FILENAME))


def create_stepfunction(x, y, dtype):
    n = 2 * len(x)
    xstep = numpy.empty(shape=n, dtype=dtype)
    # Must be zeros.
    ystep = numpy.zeros(shape=n, dtype=dtype)

    xstep[0::2] = x
    xstep[1::2] = x
    ystep[1::2] = y
    ystep[2::2] = y[:-1]

    xstep = numpy.append(xstep, xstep[-1] + 1)
    ystep = numpy.append(ystep, ystep[-1])

    return xstep, ystep


def write_rcc_data(filename, rankthreshold, rcc):
    ranknrs = numpy.arange(1, rankthreshold + 1)
    xs, ys = create_stepfunction(ranknrs, rcc, numpy.float_)

    with open(filename, 'w') as output_fh:
        prev_x, prev_y = [None] * 2
        output_fh.write("{0:d}\t{1:f}\n".format(1, 0))

        for x, y in zip(xs, ys):
            if prev_x == 1 and prev_y == 0 and x == 1 and y == 1:
                output_fh.write("{0:d}\t{1:f}\n".format(int(x), y))
            elif x == prev_x and y != prev_y:
                output_fh.write("{0:d}\t{1:f}\n".format(int(prev_x), prev_y))
                output_fh.write("{0:d}\t{1:f}\n".format(int(x), y))

            prev_x, prev_y = x, y

        output_fh.write("{0:d}\t{1:f}\n".format(int(rankthreshold + 1), y))


def write_figure(rankthreshold, dtype, rcc, avg, avg2std, p, filename, ylabel):
    fig = pylab.figure(num=None, figsize=(8, 4), dpi=200)
    ranknrs = numpy.arange(1, rankthreshold + 1)

    xstep, rccstep = create_stepfunction(ranknrs, rcc, dtype)
    plot1 = pylab.plot(xstep, rccstep, 'b-')

    xstep, rccstep = create_stepfunction(ranknrs, avg, numpy.float_)
    plot2 = pylab.plot(xstep, rccstep, 'r-')

    xstep, rccstep = create_stepfunction(ranknrs, avg2std, numpy.float_)
    plot3 = pylab.plot(xstep, rccstep, 'g-')

    pylab.plot([p[0], p[0], 1], [0, p[1], p[1]], 'k-.')
    pylab.xlabel("rank")
    pylab.ylabel("# {0:s}".format(ylabel))
    pylab.xlim([1, rankthreshold + 1])
    pylab.grid(True)
    # pylab.legend([plot1, plot2, plot3], ['Recovery curve', 'AVG', '2* STD'], loc=2, prop = {'size':12})

    figure_output = cStringIO.StringIO()

    pylab.savefig(figure_output, format='png')
    pylab.close()

    try:
        # Save figure in PNG format with a color palette (PseudoClass) instead of using truecolor (DirectClass).
        # Saving the PNG file with a color palette results in a filesize that is 3 times smaller.
        with wand.image.Image(blob=figure_output.getvalue()) as figure_image:
            figure_image.type = 'palette'
            figure_image.save(filename=filename)
    except wand.exceptions.CoderError:
        # Sometimes Wand throws the following error:
        #     CoderError: WriteBlob Failed `filename.png' @ error/png.c/MagickPNGErrorHandler/1728
        # If this happens, just write the unoptimized PNG file.
        with open(filename, 'w') as fh:
            fh.write(figure_output.getvalue())


def write_overview_figure(rccs, names, colors, rankthreshold, aucthreshold, filename, ylabel):
    fig = pylab.figure(num=None, figsize=(8, 4), dpi=200)
    ranknrs = numpy.arange(1, rankthreshold + 1)
    n = rccs.shape[1]
    plots = []

    for col_idx in numpy.arange(0, n, 2):
        xstep, rccstep = create_stepfunction(ranknrs, rccs[:, col_idx], numpy.int_)
        p = pylab.plot(xstep, rccstep, linestyle='-', color=colors[col_idx])
        plots.append(p[0])

        xstep, rccstep = create_stepfunction(ranknrs, rccs[:, col_idx + 1], numpy.float_)
        p = pylab.plot(xstep, rccstep, linestyle='-', color=colors[col_idx + 1], linewidth=1.5)
        plots.append(p[0])

    maxn = rccs.max()

    pylab.plot([aucthreshold, aucthreshold], [0, maxn], 'k-.')
    pylab.xlabel("rank")
    pylab.ylabel("# {0:s}".format(ylabel))
    pylab.xlim([1, rankthreshold + 1])
    pylab.ylim([0, maxn])
    pylab.grid(True)
    pylab.legend(plots, names, loc=1, prop={'size': 8})

    figure_output = cStringIO.StringIO()

    pylab.savefig(figure_output, format='png')
    pylab.clf()

    try:
        # Save figure in PNG format with a color palette (PseudoClass) instead of using truecolor (DirectClass).
        # Saving the PNG file with a color palette results in a filesize that is 3 times smaller.
        with wand.image.Image(blob=figure_output.getvalue()) as figure_image:
            figure_image.type = 'palette'
            figure_image.save(filename=filename)
    except wand.exceptions.CoderError:
        # Sometimes Wand throws the following error:
        #     CoderError: WriteBlob Failed `filename.png' @ error/png.c/MagickPNGErrorHandler/1728
        # If this happens, just write the unoptimized PNG file.
        with open(filename, 'w') as fh:
            fh.write(figure_output.getvalue())


def write_histogram(curves, curve_attributes, nes_threshold, filename, n_bins=100, normalize=True):
    datasets = []
    aucs = numpy.empty(shape=0, dtype=numpy.float_)

    for curve in curves:
        datasets.append(curve.aucs)
        aucs = numpy.append(aucs, curve.aucs)

    n = len(curves)
    bin_edges = numpy.linspace(start=aucs.min(), stop=aucs.max(), num=n_bins + 1, endpoint=True)
    bin_width = bin_edges[1] - bin_edges[0]

    histogram = numpy.full(shape=(n_bins, n), fill_value=0, dtype=numpy.float_ if normalize else numpy.int_)
    means = numpy.empty(shape=n)
    stdevs = numpy.empty(shape=n)

    for idx in range(n):
        ns, dummy = numpy.histogram(datasets[idx], bins=bin_edges)
        means[idx] = datasets[idx].mean()
        stdevs[idx] = datasets[idx].std()

        if normalize:
            histogram[:, idx] = ns / numpy.float_(len(datasets[idx]))
        else:
            histogram[:, idx] = ns

    max_n = histogram.max()

    plots = []

    for idx in numpy.arange(n):
        p = pylab.bar(bin_edges[:-1], histogram[:, idx], color=curve_attributes[idx].color, width=bin_width, alpha=0.5)
        auc_threshold = (nes_threshold * stdevs[idx]) + means[idx]
        pylab.plot([auc_threshold, auc_threshold], [0, max_n], color=curve_attributes[idx].color, linestyle=':')
        plots.append(p[0])

    pylab.xlabel('AUC')
    pylab.ylabel('p' if normalize else 'N')
    pylab.xlim(bin_edges[[0, -1]])
    pylab.ylim([0, max_n])
    pylab.legend(plots,
                 ["{0:s}\n$\mu={1:.4g},\sigma={2:.4g}$".format(curve_attributes[idx].name,
                                                               means[idx],
                                                               stdevs[idx])
                  for idx in numpy.arange(n)],
                 loc=1,
                 prop={'size': 8})

    figure_output = cStringIO.StringIO()

    pylab.savefig(figure_output, format='png')
    pylab.clf()

    try:
        # Save figure in PNG format with a color palette (PseudoClass) instead of using truecolor (DirectClass).
        # Saving the PNG file with a color palette results in a filesize that is 3 times smaller.
        with wand.image.Image(blob=figure_output.getvalue()) as figure_image:
            figure_image.type = 'palette'
            figure_image.save(filename=filename)
    except wand.exceptions.CoderError:
        # Sometimes Wand throws the following error:
        #     CoderError: WriteBlob Failed `filename.png' @ error/png.c/MagickPNGErrorHandler/1728
        # If this happens, just write the unoptimized PNG file.
        with open(filename, 'w') as fh:
            fh.write(figure_output.getvalue())


def write_clusters(motifs, clusters, filename):
    if not clusters:
        return

    with open(filename, 'w') as output_fh:
        for row in cluster_table_generator(motifs, clusters):
            print >> output_fh, "{0:d}\t{1:d}\t{2:s}".format(*row)


def get_all_subdirs_and_filenames_for_directory(root_directory_name):
    subdirs_list = list()
    filenames_list = list()

    # Strip off trailing slashes.
    root_directory_name = os.path.normpath(root_directory_name)

    # Get the length of the root directory name and add 1 for the trailing slash.
    length_root_directory_name = len(root_directory_name) + 1

    for folder, subdirs, filenames in os.walk(root_directory_name):
        # Make the folder relative to the results_folder.
        folder = folder[length_root_directory_name:]

        for subdir in subdirs:
            subdirs_list.append(os.path.join(folder, subdir))

        for filename in filenames:
            filenames_list.append(os.path.join(folder, filename))

    return root_directory_name, subdirs_list, filenames_list


def create_archive(folder_name):
    filenames = get_all_subdirs_and_filenames_for_directory(folder_name)[2]

    zip_filename = os.path.join(folder_name, ARCHIVE_FILENAME)

    with zipfile.ZipFile(zip_filename, mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True) as zip_output_fh:
        for filename in filenames:
            # Do not include archive.zip in archive.zip in case it already would exist, which could happen if a job
            # is restarted after a crash of i-cisTarget.
            if filename != ARCHIVE_FILENAME:
                zip_output_fh.write(os.path.join(folder_name, filename), arcname=os.path.join('icistarget', filename))


def write_tsv(filename, ranked_gene_or_region_ids_with_description):
    with open(filename, 'w') as output_fh:
        for ranked_gene_or_region_id_with_description in ranked_gene_or_region_ids_with_description:
            rank, gene_or_region_id, gene_or_region_id_description = ranked_gene_or_region_id_with_description
            print >> output_fh, "{0:d}\t{1:s}\t{2:s}".format(int(ranked_gene_or_region_id_with_description['rank']),
                                                             ranked_gene_or_region_id_with_description['id'],
                                                             ranked_gene_or_region_id_with_description['description'])


def generate_filename(folder, filename, extension):
    full_filename = filename + "." + extension
    full_path = os.path.join(folder, full_filename)

    if len(full_path) <= 255:
        return full_path, full_filename
    else:
        full_filename = hashlib.sha1(filename).hexdigest() + "." + extension

        return os.path.join(folder, full_filename), full_filename


def generate_report(curves,
                    context,
                    output_folder,
                    title=None,
                    true_target_ids=[],
                    ini_filename=None,
                    elements_to_highlight=None,
                    template_variables=dict()):
    assert len(curves) > 0

    # Setup output folder.
    gene_or_region_ids = list(curves[0].database.gene_ids)
    install_output_folder(output_folder, ini_filename)
    logos_folder = os.path.join(output_folder, 'logos')

    # Write gene or region IDs to disk.
    write_gene_or_region_ids(gene_or_region_ids, os.path.join(output_folder, QUERY_IDS_FILENAME))

    # Write statistics to disk.
    with open(os.path.join(output_folder, STATISTICS_FILENAME), 'w') as output_fh:
        name = re.sub('\s+', '_', title.strip())
        logging.info("Writing statistics table")
        write_full_statistics_as_tsv(name,
                                     curves,
                                     output_fh,
                                     feature_id2description=ContextWrapper(context),
                                     feature_id2gene_annotations=ContextWrapper(context, 'gene_ids', False))

    # Assemble parameters.
    class Parameters:
        pass

    params = Parameters()
    params.total_feature_count = sum(map(lambda curve: curve.database.feature_count, curves))
    params.nes_threshold = curves[0].nes_threshold
    params.enriched_feature_count = sum(map(lambda curve: curve.enriched_feature_count_without_randoms, curves))
    params.total_gene_id_count = curves[0].database.total_gene_count
    params.gene_id_count = curves[0].database.gene_count
    params.gene_ids_filename = QUERY_IDS_FILENAME
    params.nes_threshold = curves[0].nes_threshold
    params.auc_threshold = curves[0].auc_threshold
    params.rank_cutoff = curves[0].rank_cutoff
    params.rank_threshold = curves[0].rank_threshold
    params.has_true_target_ids = len(true_target_ids) > 0

    # Generate histogram of AUC distributions.
    curve_attributes = []
    CurveAttributes = collections.namedtuple('CurveAttributes', 'color name')

    for idx, curve in enumerate(curves):
        curve_attributes.append(
            CurveAttributes(color=GROUP_COLORS[idx % GROUP_COLORS_LENGTH], name=curves[idx].database.name))

    histogram_filename = os.path.join(output_folder, HIST_FILENAME)
    write_histogram(curves, curve_attributes, params.nes_threshold, histogram_filename)

    # Clustering analysis.
    enriched_features = []

    for curve in curves:
        enriched_features += list(curve.enriched_features)

    enriched_features = map(operator.itemgetter(0), sorted(enriched_features, key=operator.itemgetter(2), reverse=True))

    context.write_clusterbuster_motifs(feature_ids=enriched_features,
                                       cb_motifs_output_filename=os.path.join(output_folder, MOTIFS_CB_FILENAME),
                                       include_transfac_pro=False)

    if context.has_stamp():
        logging.info("clustering motifs")
        clusters, missing_features = cluster(context.get_stamp_motif_database_filename(),
                                             context.get_stamp_score_distribution_filename(),
                                             enriched_features,
                                             cmd=context.stamp_command,
                                             has_clusters=context.has_cluster())

        if context.has_cluster():
            logging.info("clustering tracks")
            track_clusters, track_garbage = track_cluster(context.get_track_cluster_filename(), missing_features)
            clusters.extend(track_clusters)

        write_clusters(enriched_features, clusters, os.path.join(output_folder, CLUSTERS_FILENAME))
        feature_id2color = get_color_mapping(enriched_features, clusters)
        context.remove_stamp_intermediate_files()
    elif context.has_cluster():
        logging.info("clustering tracks")
        clusters, garbage = track_cluster(context.get_track_cluster_filename(), enriched_features)
        write_clusters(enriched_features, clusters, os.path.join(output_folder, CLUSTERS_FILENAME))
        feature_id2color = get_color_mapping(enriched_features, clusters)
    else:
        feature_id2color = dict()

    # Generate overview RCC curves.
    rcc_overview_filename = os.path.join(output_folder, OVERVIEW_FILENAME)
    best_rccs = numpy.empty(shape=(params.rank_threshold, 0), dtype=numpy.int_)
    best_feature_ids = []
    best_colors = []

    for idx, curve in enumerate(curves):
        if not curve.enriched_feature_count:
            continue

        best_feature_id = max(curve.enriched_features, key=operator.itemgetter(2))[0]
        curve.get_recovery_curve(best_feature_id)

        best_rccs = numpy.append(best_rccs, curve.get_recovery_curve(best_feature_id)[:, numpy.newaxis], axis=1)
        best_feature_ids.append(best_feature_id)
        best_colors.append(curve_attributes[idx].color)

        best_rccs = numpy.append(best_rccs, curve.avgrcc[:, numpy.newaxis], axis=1)
        best_feature_ids.append("avg(" + best_feature_id + ")")
        best_colors.append(curve_attributes[idx].color)
    if numpy.prod(best_rccs.shape):
        write_overview_figure(best_rccs, best_feature_ids, best_colors, params.rank_threshold, params.rank_cutoff,
                              rcc_overview_filename, context.ylabel)

    # Write recovery curves for average and avg2std.
    for idx, curve in enumerate(curves):
        name = curve_attributes[idx].name

        if len(name) > 150:
            name = 'database{0:d}'.format(idx + 1)

        rcc_data_filename = os.path.join(output_folder, name + ".avg." + DATA_EXTENSION)
        write_rcc_data(rcc_data_filename, params.rank_threshold, curves[idx].avgrcc)

        rcc_data_filename = os.path.join(output_folder, name + ".avg2std." + DATA_EXTENSION)
        write_rcc_data(rcc_data_filename, params.rank_threshold, curves[idx].avg2stdrcc)

    class Statistics:
        pass

    stats = Statistics()
    stats.histogram_filename = HIST_FILENAME
    stats.rcc_overview_filename = OVERVIEW_FILENAME

    class Row:
        def __init__(self, rank, feature_id, nes, database):
            self.rank = rank
            self.feature_id = feature_id
            self.feature_description = context.get_description_for_feature(feature_id)
            self.feature_gene_annotations = ", ".join(context.get_gene_annotations_for_feature(feature_id))
            self.nes = format(nes, '.5f')
            self.database = re.sub('>', '&gt;', re.sub('<', '&lt;', database))

    data = []
    iterators = map(RecoveryCurves.iterate_enriched_features, curves)
    database_name2nes_threshold = dict((curve.database.name, curve.adjusted_enrichement_threshold)
                                       for curve in curves
                                       if curve.has_no_enriched_features)
    name2database_in_mem = dict()

    for idx, (feature_id, auc, nes, curve) in enumerate(sorted(itertools.chain(*iterators),
                                                               key=operator.itemgetter(2),
                                                               reverse=True)):
        logging.info("analysis of significant ranking {0:s} (NES = {1:.10g}, {2:d}/{3:d})".format(
            feature_id,
            nes,
            idx + 1,
            params.enriched_feature_count)
        )

        database_name = curve.database.get_group_name(feature_id)
        database_output_folder = os.path.join(output_folder, database_name)

        if not os.path.isdir(database_output_folder):
            os.mkdir(database_output_folder)

        if curve.has_enriched_random_features and database_name not in database_name2nes_threshold:
            database_name2nes_threshold[database_name] = curve.adjusted_enrichement_threshold

        row = Row(idx + 1, feature_id, nes, database_name)

        candidate_target_gene_or_region_ids = curve.get_candidate_target_gene_ids(feature_id)
        top_ranked_region_or_gene_ids = curve.get_top_ranked_gene_ids(feature_id)

        rcc = curve.get_recovery_curve(feature_id)
        cp = curve.get_critical_point(feature_id)

        row.color = feature_id2color.get(feature_id, 'FFFFFF')

        row.logo_filename = context.copy_logo_for_feature_to(feature_id, logos_folder)

        full_path, full_filename = generate_filename(database_output_folder, feature_id, FIGURE_EXTENSION)

        write_figure(params.rank_threshold, curve.database.dtype, rcc, curve.avgrcc, curve.avg2stdrcc, cp, full_path,
                     context.ylabel)

        row.rcc_figure_filename = full_filename

        full_path = generate_filename(database_output_folder, feature_id, DATA_EXTENSION)[0]
        write_rcc_data(full_path, params.rank_threshold, rcc)

        if context.gene_or_region_id_locations:
            full_path, candidate_targets_bed_filename = generate_filename(database_output_folder,
                                                                          feature_id,
                                                                          "targets." + BED_EXTENSION)
            track_name = "{0:s}_candidate_targets".format(feature_id)

            context.gene_or_region_id_locations.filter_by_name(
                map(operator.itemgetter(1), candidate_target_gene_or_region_ids)).write(full_path, track_name)
        else:
            candidate_targets_bed_filename = None

        namespace = {
            'params': {
                'recovered_ids_count': cp[1],
                'leading_edge_rank_cutoff': cp[0],
                'total_id_count': params.total_gene_id_count
            },
            'data': context.get_ranked_gene_or_region_ids_for_template(candidate_target_gene_or_region_ids,
                                                                       elements_to_highlight),
            'bed': context.get_locations_for_template(candidate_target_gene_or_region_ids),
            'database_name': database_name,
            'bedfilename': candidate_targets_bed_filename,
            'feature_id': feature_id,
            'feature_description': context.get_description_for_feature(feature_id)
        }

        namespace = dict(namespace.items() + template_variables.items())
        full_path, candidate_targets_filename = generate_filename(database_output_folder,
                                                                  feature_id,
                                                                  "targets." + HTML_EXTENSION)
        write_template(context.candidate_targets_template_filename, full_path, namespace)
        row.candidate_targets_filename = candidate_targets_filename
        full_path = generate_filename(database_output_folder, feature_id, "targets." + TSV_EXTENSION)[0]

        write_tsv(full_path, context.get_ranked_gene_or_region_ids_for_tsv(candidate_target_gene_or_region_ids))

        if context.gene_or_region_id_locations:
            full_path, top_ranked_genes_bed_filename = generate_filename(
                database_output_folder,
                feature_id,
                "top" + str(params.rank_threshold) + "." + BED_EXTENSION)

            track_name = "{0:s}_top{1:d}_targets".format(feature_id, params.rank_threshold)

            context.gene_or_region_id_locations.filter_by_name(map(operator.itemgetter(1),
                                                                   top_ranked_region_or_gene_ids)
                                                               ).write(full_path, track_name)
        else:
            top_ranked_genes_bed_filename = None

        namespace = {
            'data': context.get_ranked_gene_or_region_ids_for_template(top_ranked_region_or_gene_ids,
                                                                       elements_to_highlight),
            'bed': context.get_locations_for_template(top_ranked_region_or_gene_ids),
            'database_name': database_name,
            'bedfilename': top_ranked_genes_bed_filename,
            'feature_id': feature_id,
            'feature_description': context.get_description_for_feature(feature_id)
        }

        namespace = dict(namespace.items() + template_variables.items())
        full_path, top_ranked_genes_filename = generate_filename(
            database_output_folder,
            feature_id,
            "top" + str(params.rank_threshold) + "." + HTML_EXTENSION)

        write_template(context.top_ranked_targets_template_filename, full_path, namespace)

        row.top_ranked_genes_filename = top_ranked_genes_filename
        full_path = generate_filename(database_output_folder,
                                      feature_id,
                                      "top" + str(params.rank_threshold) + "." + TSV_EXTENSION)[0]

        write_tsv(full_path, context.get_ranked_gene_or_region_ids_for_tsv(top_ranked_region_or_gene_ids))

        if true_target_ids:
            ppv, tpr, tp, fn, fp = curve.get_ppv_tpr(feature_id, true_target_ids)

            row.ppv = format(ppv, '.5f') if not math.isnan(ppv) else '-'
            row.ppv_numerator = tp
            row.ppv_denominator = tp + fp

            row.tpr = format(tpr, '.5f') if not math.isnan(tpr) else '-'
            row.tpr_numerator = tp
            row.tpr_denominator = tp + fn

        data.append(row)

    for name, database in name2database_in_mem.iteritems():
        logging.info("Clearing top ranked genes of database '{0:s}' from memory".format(name))
        database.clear_top_regions()

    # Generate warning.
    if database_name2nes_threshold:
        def to_string(database_name_nes_threshold):
            database_name, nes_threshold = database_name_nes_threshold

            return '{0:s} ({1:.5f})'.format(database_name, nes_threshold)

        warning_str = 'For some databases, randomly generated features were enriched. ' \
                      + 'Also, non-random features in these same databases that have a lower enrichment ' \
                      + 'than the first randomly generated enriched feature were removed from the report. ' \
                      + 'The involved databases and their adjusted NES thresholds are: {0:s}.'.format(
            ', '.join(itertools.imap(to_string, database_name2nes_threshold.iteritems())))
    else:
        warning_str = ''

    # Fill in velocity template.
    namespace = {'title': title,
                 'archive': ARCHIVE_FILENAME,
                 'params': params,
                 'stats': stats,
                 'data': data,
                 'general': template_variables,
                 'warning': warning_str}

    write_template(context.main_template_filename,
                   os.path.join(output_folder, HTML_REPORT_FILENAME),
                   namespace)

    # Create archive.
    create_archive(output_folder)


def generate_batch_report(name2curves, output_folder, context, template_variables=dict()):
    variables = dict(template_variables.iteritems())
    variables["statsfilename"] = BATCH_RESULTS_FILENAME

    with open(os.path.join(output_folder, HTML_REPORT_FILENAME), 'w') as output_fh:
        write_overview_as_html(name2curves.iteritems(),
                               output_fh,
                               context.batch_report_template_filename,
                               template_variables=variables,
                               feature_id2description=ContextWrapper(context),
                               feature_id2gene_annotations=ContextWrapper(context, 'gene_ids', False))

    with gzip.open(os.path.join(output_folder, BATCH_RESULTS_FILENAME), 'wb') as gzip_output_fh:
        write_overview_as_tsv(name2curves.iteritems(),
                              gzip_output_fh,
                              feature_id2description=ContextWrapper(context),
                              feature_id2gene_annotations=ContextWrapper(context, 'gene_ids', False))

    # Create archive.
    create_archive(output_folder)
