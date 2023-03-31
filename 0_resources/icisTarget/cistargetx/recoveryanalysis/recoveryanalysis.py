import ConfigParser
import getopt
import logging
import operator
import os
import sys
import timeit
from multiprocessing import Pool, Queue
from signal import signal, SIGPIPE, SIG_DFL

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

from cistargetx.common.rankingsdatabase import load_databases
from cistargetx.common.signatureformats import LEGACY, GMT_COMMA, GMT_SEMICOLON, GMT_TAB, GeneSignature
from recoverycurves import RecoveryCurves
from reportcontext import ReportContext
from reportgeneration import generate_report
from utils import ContextWrapper, write_overview_as_html, write_overview_as_tsv


INI_TEMPLATE = os.path.join(os.path.dirname(__file__), 'recoveryanalysis-template.ini')


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, "Usage: python {0:s} [-o html_report] <inifile> [<genesetfile>] <outputdir>".format(
        command_name)
    print >> output_fh, "       python {0:s} -o [tsv_overview|html_overview] <inifile> <genesetfile1> <genesetfile2> ...".format(
        command_name)
    print >> output_fh, "       python {0:s} -o filter <inifile> <genesetfile1> <genesetfile2> ...".format(command_name)
    print >> output_fh
    print >> output_fh, "Optional parameters: -n <NES-threshold> (overrides value in ini-file)"
    print >> output_fh, "                     -r <rank-threshold> (overrides value in ini-file)"
    print >> output_fh, "                     -a <AUC-percentage-threshold> (overrides value in ini-file)"
    print >> output_fh, "                     -c <feature1>;<feature2>;... (overrides value in ini-file)"
    print >> output_fh, "                     -t <report-title> : the title to be used in the HTML report."
    print >> output_fh, "                     -f gmt(=gmt_comma)|gmt_tab|gmt_comma|gmt_semicolon|legacy|auto : the gene signature input format; default is auto detection of format."
    print >> output_fh, "When supplying multiple filenames with genesets the results table will have an additional column containing the ID of the geneset. The ID is taken as the filename without extension."
    print >> output_fh, "One of the geneset filenames can be 'stdin'."
    print >> output_fh
    print >> output_fh, "Example ini file:"
    print >> output_fh

    with open(INI_TEMPLATE, 'r') as input_fh:
        print >> output_fh, input_fh.read()


def get_option(cfg, section, options):
    for option in options:
        if cfg.has_option(section, option):
            return cfg.get(section, option)

    return None


def has_option(cfg, section, options):
    return any(cfg.has_option(section, option) for option in options)


class Parameters:
    pass


def parse_arguments(args):
    try:
        opts, args = getopt.getopt(args, "ho:n:r:a:t:c:f:l")
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err)
        display_usage(sys.stderr)
        sys.exit(2)

    if len(args) < 2:
        print >> sys.stderr, "Wrong number of input arguments."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(2)

    p = Parameters()

    p.inifile = args[0]

    if not os.path.isfile(p.inifile):
        print >> sys.stderr, "'{0:s}' file doesn't exist.".format(p.inifile)
        sys.exit(1)

    cfg = ConfigParser.RawConfigParser()
    cfg.read(p.inifile)

    p.outputformat = 'html_report'
    p.rank_threshold = cfg.getint('params', 'rank_threshold')
    p.auc_threshold = cfg.getfloat('params', 'auc_threshold')
    p.nes_threshold = (cfg.getfloat('params', 'nes_threshold')
                       if cfg.has_option('params', 'nes_threshold')
                       else cfg.getfloat('params', 'enrichment_score_threshold'))
    p.title = cfg.get('report', 'title') if cfg.has_option('report', 'title') else ''
    p.motifcombinations = (get_option(cfg, 'params', ('feature_combinations', 'motif_recombinations')).split(';')
                           if has_option(cfg, 'params', ('feature_combinations', 'motif_recombinations'))
                           else [])

    code2format = {"auto": None,
                   "gmt": GMT_COMMA,
                   "legacy": LEGACY,
                   "gmt_tab": GMT_TAB,
                   "gmt_comma": GMT_COMMA,
                   "gmt_semicolon": GMT_SEMICOLON}
    # Detect gene signature input format automatically.
    p.format = code2format["auto"]

    p.rccdata = False

    for o, a in opts:
        if o == "-o":
            if a.lower() not in ('filter', 'html_report', 'html_overview', 'tsv_overview'):
                print >> sys.stderr, "'{0:s}' is not a valid output format.".format(a)
                display_usage(sys.stderr)
                sys.exit(1)
            p.outputformat = a
        elif o == "-h":
            display_usage()
            sys.exit()
        elif o == "-r":
            p.rank_threshold = int(a)
        elif o == "-a":
            p.auc_threshold = float(a)
        elif o == "-n":
            p.nes_threshold = float(a)
        elif o == "-t":
            p.title = a
        elif o == "-c":
            p.motifcombinations = a.split(';')
        elif o == "-f":
            if a not in code2format:
                print >> sys.stderr, "'{0:s}' is not a valid gene signature input format.".format(a)
                display_usage(sys.stderr)
                sys.exit(1)

            p.format = code2format[a]
        else:
            assert False, "Unhandled option"

    if p.outputformat == 'html_report':
        if len(args) == 2:
            p.genesetfilenames = [get_option(cfg, 'data', ('ids_filename', 'gene_set_filename'))]
            p.outputdir = args[1]
        elif len(args) == 3:
            p.genesetfilenames = [args[1]]
            p.outputdir = args[2]
        else:
            print >> sys.stderr, "Wrong number of input arguments."
            display_usage(sys.stderr)
            sys.exit(2)
    else:
        p.outputdir = os.path.curdir
        p.genesetfilenames = args[1:]

    return p, cfg


def init_logging(log_filename=None):
    logging.basicConfig(stream=sys.stderr, format='%(message)s', level=logging.DEBUG)

    if log_filename:
        log_file_handler = logging.FileHandler(log_filename)
        log_file_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        log_file_handler.setFormatter(formatter)
        logging.root.addHandler(log_file_handler)


def load_gene_ids(filename):
    with open(filename, 'r') if filename != 'stdin' else sys.stdin as input_fh:
        return set(line.rstrip() for line in input_fh.readlines())


def derive_id(filename):
    if filename == 'stdin':
        return filename

    return os.path.basename(os.path.splitext(filename)[0])


def get_databases(cfg, p, gene_ids):
    filenames = [os.path.expandvars(db.strip().rstrip(';'))
                 for db in cfg.get('data', 'ranking_dbs').split('\n')
                 if db and not (db.strip().startswith('#') or db.strip().startswith(';'))]

    if (has_option(cfg, 'params', ('feature_masks', 'motif_masks'))
            and get_option(cfg, 'params', ('feature_masks', 'motif_masks')) != '*'):

        feature_masks = get_option(cfg, 'params', ('feature_masks', 'motif_masks')).rstrip().split(";")
    else:
        feature_masks = []

    features_to_combine = p.motifcombinations
    enrichment_within_db = (cfg.has_option('params', 'enrichment_within_db')
                            and cfg.getboolean('params', 'enrichment_within_db'))
    control_db = (os.path.expandvars(get_option(cfg, 'data', ('control_db', 'recombinations_control_db')))
                  if has_option(cfg, 'data', ('control_db', 'recombinations_control_db'))
                  else None)
    combination_method = (cfg.get('params', 'combination_method')
                          if cfg.has_option('params', 'combination_method')
                          else 'pairwise')

    return load_databases(filenames,
                          gene_ids,
                          feature_masks,
                          features_to_combine,
                          enrichment_within_db,
                          control_db,
                          name_translation=derive_id,
                          combination_method=combination_method)


def generate_context(cfg):
    if cfg.has_option('data', 'rankings_metadata'):
        context = ReportContext.create_from_zip_archive(
            (cfg.get('report', 'id_type')
             if cfg.has_option('report', 'id_type')
             else 'regions'),

            (get_option(cfg, 'report', ('id_description_filename', 'gene_description_file'))
             if has_option(cfg, 'report', ('id_description_filename', 'gene_description_file'))
             else None),

            (cfg.get('report', 'id_description_highlighting_filename')
             if cfg.has_option('report', 'id_description_highlighting_filename')
             else None),

            (cfg.get('report', 'id_location_bed_filename')
             if cfg.has_option('report', 'id_location_bed_filename')
             else None),

            os.path.expandvars(cfg.get('data', 'rankings_metadata')),

            track_cluster_filename=(os.path.expandvars(cfg.get('data', 'track_cluster_filename'))
                                    if cfg.has_option('data', 'track_cluster_filename')
                                    else None),

            stamp_command=(cfg.get('STAMP', 'command')
                           if cfg.has_option('STAMP', 'command')
                           else None),

            motif_annotations_filename=(os.path.expandvars(cfg.get('data', 'motif_annotations_filename'))
                                        if cfg.has_option('data', 'motif_annotations_filename')
                                        else None),

            track_annotations_filename=(os.path.expandvars(cfg.get('data', 'track_annotations_filename'))
                                        if cfg.has_option('data', 'track_annotations_filename')
                                        else None)
        )
    else:
        context = ReportContext.create_from_files(
            (cfg.get('report', 'id_type') if cfg.has_option('report', 'id_type') else 'regions'),

            (get_option(cfg, 'report', ('id_description_filename', 'gene_description_file'))
             if has_option(cfg, 'report', ('id_description_filename', 'gene_description_file'))
             else None),

            (cfg.get('report', 'id_description_highlighting_filename')
             if cfg.has_option('report', 'id_description_highlighting_filename')
             else None),

            (cfg.get('report', 'id_location_bed_filename')
             if cfg.has_option('report', 'id_location_bed_filename')
             else None),

            (os.path.expandvars(cfg.get('report', 'logo_folder'))
             if cfg.has_option('report', 'logo_folder')
             else None),

            (cfg.get('report', 'logo_extension')
             if cfg.has_option('report', 'logo_extension')
             else None),

            track_cluster_filename=(os.path.expandvars(cfg.get('data', 'track_cluster_filename'))
                                    if cfg.has_option('data', 'track_cluster_filename')
                                    else None),

            stamp_command=(cfg.get('STAMP', 'command')
                           if cfg.has_option('STAMP', 'command')
                           else None),

            stamp_motif_database_filename=(cfg.get('STAMP', 'transfac_motif_database')
                                           if cfg.has_option('STAMP', 'transfac_motif_database')
                                           else None),

            stamp_score_distribution_filename=(cfg.get('STAMP', 'score_database')
                                               if cfg.has_option('STAMP', 'score_database')
                                               else None),

            motif_annotations_filename=(os.path.expandvars(cfg.get('data', 'motif_annotations_filename'))
                                        if cfg.has_option('data', 'motif_annotations_filename')
                                        else None),

            track_annotations_filename=(os.path.expandvars(cfg.get('data', 'track_annotations_filename'))
                                        if cfg.has_option('data', 'track_annotations_filename')
                                        else None)
        )

    return context


def get_true_target_ids(cfg):
    if has_option(cfg, 'data', ('true_target_ids_filename', 'crm_regions_filename')):
        true_target_ids = load_gene_ids(os.path.expandvars(get_option(cfg, 'data', ('true_target_ids_filename',
                                                                                    'crm_regions_filename'))))
    else:
        true_target_ids = set()

    return true_target_ids


def write_filtered_gene_ids(named_curves_iterable, output_fh=sys.stdout):
    for name, curves in named_curves_iterable:
        gene_ids = set()

        for curve in curves:
            for feature_id, auc, nes, curve in curve.iterate_enriched_features():
                gene_ids |= set(map(operator.itemgetter(1), curve.get_candidate_target_gene_ids(feature_id)))

        for gene_id in gene_ids:
            print >> output_fh, "{0:s}\t{1:s}".format(name, gene_id)


def generate_named_curves(id2gene_signatures, p, cfg):
    for name, gene_ids in id2gene_signatures.iteritems():
        logging.info("{0:s}: Analyzing ...".format(name))

        try:
            curves = [RecoveryCurves(db, p.auc_threshold, p.rank_threshold, p.nes_threshold)
                      for db in get_databases(cfg, p, gene_ids)]
            logging.info("{0:s}: Fraction of genes/regions mappable to the database(s) = {1:.2f}%".format(
                name,
                min(curve.fraction_of_mapped_gene_ids(gene_ids) for curve in curves) * 100.0))

            yield name, curves
        except AssertionError, msg:
            logging.warning("{0:s}: The following error occurred '{1:s}'".format(name, msg))


def produce_queue_item(name, gene_ids, niceness=19):
    os.nice(niceness)

    try:
        for name, curves in generate_named_curves({name: gene_ids}, produce_queue_item.p, produce_queue_item.cfg):
            produce_queue_item.named_curves_queue.put((name, curves))
    except:
        # Make sure an item is put on the queue.
        produce_queue_item.named_curves_queue.put((name, []))


def produce_queue_item_initializer(p, cfg, named_curves_queue):
    produce_queue_item.p = p
    produce_queue_item.cfg = cfg
    produce_queue_item.named_curves_queue = named_curves_queue


def main():
    start_time = timeit.default_timer()

    # Ignore SIG_PIPE and don't throw exceptions on it:
    #   http://docs.python.org/library/signal.html
    signal(SIGPIPE, SIG_DFL)

    #
    # Parsing of input parameters and assembling all metadata into a context.
    #
    p, cfg = parse_arguments(sys.argv[1:])
    init_logging(cfg.get('params', 'log_filename') if cfg.has_option('params', 'log_filename') else None)
    context = generate_context(cfg) if p.outputformat in ('html_report', 'html_overview', 'tsv_overview') else None

    #
    # Loading the gene signatures into memory.
    #
    id2gene_signatures = dict()
    for filename in p.genesetfilenames:
        with open(filename, 'r') if filename != 'stdin' else sys.stdin as input_fh:
            for signature in GeneSignature.load_from_stream(input_fh, derive_id(filename), p.format):
                if signature.name in id2gene_signatures:
                    assert False, "ID '{0:s}' for gene signature is not unique.".format(signature.name)
                else:
                    id2gene_signatures[signature.name] = signature.identifiers

    #
    # Actual analysis.
    #
    if p.outputformat == 'html_report':
        if len(id2gene_signatures) != 1:
            logging.error("When requesting a HTML report only one signature can be provided.")
            sys.exit(1)

        name = id2gene_signatures.keys()[0]
        gene_ids = id2gene_signatures[name]
        logging.info("{0:s}: Analyzing ...".format(name))
        curves = [RecoveryCurves(db, p.auc_threshold, p.rank_threshold, p.nes_threshold)
                  for db in get_databases(cfg, p, gene_ids)]
        logging.info("{0:s}: Fraction of genes/regions mappable to the database(s) = {1:.2f}%".format(
            name,
            min(curve.fraction_of_mapped_gene_ids(gene_ids) for curve in curves) * 100.0))

        generate_report(curves,
                        context,
                        p.outputdir,
                        title=p.title,
                        true_target_ids=get_true_target_ids(cfg),
                        ini_filename=p.inifile)
    elif p.outputformat in ("html_overview", "tsv_overview", "filter"):
        n_cores = cfg.getint('params', 'n_cores') if cfg.has_option('params', 'n_cores') else 1
        niceness = cfg.getint('params', 'process_niceness') if cfg.has_option('params', 'process_niceness') else 19

        if n_cores > 1:
            # One of the imported python modules sets the CPU affinity to one core.
            # If we do not recent the CPU affinity all processes (python and cbust) wil run
            # on the same core, defeating the whole purpose of the multiprocessing module.
            #
            # Most likely one of the imported modules in linked with OpenBLAS:
            #     http://stackoverflow.com/questions/23537716/importing-scipy-breaks-multiprocessing-support-in-python
            #
            # If we ever upgrade to Python 3, the call to taskset can be replaced with:
            #     os.sched_setaffinity(pid, mask)
            DEVNULL = open(os.devnull, 'w')
            subprocess.call(['taskset','-p', '0xFFFFFFFF', str(os.getpid())], stdout=DEVNULL)
            DEVNULL.close()

            # To avoid memory overload, the size of the queue is limited.
            named_curves_queue = Queue(maxsize=100)
            pool = Pool(processes=n_cores - 1,
                        initializer=produce_queue_item_initializer,
                        initargs=(p, cfg, named_curves_queue))

            for name, gene_signature in id2gene_signatures.iteritems():
                pool.apply_async(produce_queue_item, (name, gene_signature, niceness))

            pool.close()

            # Because the number of elements that need to be removed from queue is known in advance,
            # the parent process can be used to consume results from the queue.
            total_signature_count = len(id2gene_signatures)

            for idx in range(total_signature_count):
                name, curves = named_curves_queue.get()

                try:
                    logging.warning("{0:s}: Current items in the queue = {1:d}".format(name,
                                                                                       named_curves_queue.qsize()))
                except NotImplementedError:
                    pass

                if not curves:
                    continue

                if p.outputformat == 'tsv_overview':
                    write_overview_as_tsv([(name, curves)],
                                          sys.stdout,
                                          feature_id2description=ContextWrapper(context),
                                          feature_id2gene_annotations=ContextWrapper(context, 'gene_ids', False))
                elif p.outputformat == 'html_overview':
                    write_overview_as_html([(name, curves)],
                                           sys.stdout,
                                           context.batch_report_template_filename,
                                           feature_id2description=ContextWrapper(context),
                                           feature_id2gene_annotations=ContextWrapper(context, 'gene_ids', False))
                else:
                    write_filtered_gene_ids([(name, curves)])

            pool.join()
        else:
            # To avoid memory overload when analyzing large datasets, a producer-consumer pattern is used.
            # This is implemented via a python generator function.
            # Decide which consumer function to use.
            if p.outputformat == 'tsv_overview':
                def _wrapper(named_curves_iterable):
                    write_overview_as_tsv(named_curves_iterable,
                                          sys.stdout,
                                          feature_id2description=ContextWrapper(context),
                                          feature_id2gene_annotations=ContextWrapper(context, 'gene_ids', False))

                named_curve_consumer = _wrapper
            elif p.outputformat == 'html_overview':
                def _wrapper(named_curves_iterable):
                    write_overview_as_html(named_curves_iterable,
                                           sys.stdout,
                                           context.batch_report_template_filename,
                                           feature_id2description=ContextWrapper(context),
                                           feature_id2gene_annotations=ContextWrapper(context, 'gene_ids', False))

                named_curve_consumer = _wrapper
            else:
                named_curve_consumer = write_filtered_gene_ids

            named_curve_consumer(generate_named_curves(id2gene_signatures, p, cfg))
    else:
        logging.error("Invalid output format: '{0:s}'.".format(p.outputformat))

    logging.info("Elapsed time = {0:f}s".format(timeit.default_timer() - start_time))


if __name__ == "__main__":
    main()
