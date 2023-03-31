import StringIO
import getopt
import logging
import os
import sys
import timeit
import traceback
from ConfigParser import RawConfigParser
from multiprocessing import Pool, Queue, Process, Value, current_process

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

import airspeed

from bunch import Bunch
from cistargetx.iregulon.daemon import load_nomenclature2code
from cistargetx.recoveryanalysis.recoveryanalysis import get_databases, init_logging, get_option, has_option
from cistargetx.recoveryanalysis.recoverycurves import RecoveryCurves
from gmtparser import MetadataEnrichedGeneSignature
from mysql import create_db_connection, insert_collection_metadata, insert_gene_signature, insert_analysis, \
    remove_collection, CollectionAlreadyPresentException


DATABASE_CONFIGURATION_FILE = os.path.join(os.path.dirname(__file__), 'seq-srv-01.ini')
DATABASE_DATAMODEL_TEMPLATE_FILE = os.path.join(os.path.dirname(__file__), 'datamodel.vm')

# TODO: When uploading using multiple cores: some analysis results are not inserted into the DB??
# TODO: Known loss is: Only a single gene could be mapped to the database.


def display_usage(output_fh=sys.stderr, command_name=os.path.basename(sys.argv[0])):
    print >> output_fh, "Usage: python {0:s} [<options>] <genesetfile>".format(command_name)
    print >> output_fh, "       python {0:s} [<options>] init".format(command_name)
    print >> output_fh, "Option: --cfg=<database_configuration_file> : The database configuration file to use."
    print >> output_fh, "        --operation=<new|add|rm|complete>   : Operation to perform for supplied collection (default is new)."
    print >> output_fh, "        --duplicates=<rm_old|rm_new|none>   : Operation to perform when duplicate signature based on PubMed ID is encountered (default is none)."
    print >> output_fh, "        --motif2tf=<tf_motif_mapping|tf_motif_mapping_v6_1>   : motif2tf database (default is tf_motif_mapping corresponding to motif collection v3)."
    print >> output_fh
    print >> output_fh, "WARNING: DO NOT USE MULTIPLE CORES. Some analysis results are not inserted into the DB"
    print >> output_fh
    print >> output_fh, "Example configuration file:"
    print >> output_fh
    with open(DATABASE_CONFIGURATION_FILE, 'r') as input_fh:
        print >> output_fh, input_fh.read()


def parse_arguments(input_args):
    opts, args = getopt.getopt(input_args, '', ['cfg=', 'operation=', 'duplicates=', 'motif2tf='])
    if len(args) != 1:
        print >> sys.stderr, "Wrong number of input arguments."
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    gmt_filename = args[0]
    if not os.path.isfile(gmt_filename) and not gmt_filename == 'init':
        print >> sys.stderr, "'{0:s}' file doesn't exist.".format(gmt_filename)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    cfg_filename = DATABASE_CONFIGURATION_FILE
    operation = 'new'
    duplicate = 'none'
    motif2tf = 'tf_motif_mapping'
    for o, a in opts:
        if o == '--cfg':
            cfg_filename = a
        elif o == '--operation':
            operation = a
        elif o == '--duplicates':
            duplicate = a
        elif o == '--motif2tf':
            motif2tf = a

    if not os.path.isfile(cfg_filename):
        print >> sys.stderr, "'{0:s}' file doesn't exist.".format(cfg_filename)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    db_cfg = RawConfigParser()
    db_cfg.read(cfg_filename)

    if operation not in ('new', 'add', 'rm', 'complete'):
        print >> sys.stderr, "'{0:s}' is an invalid operation.".format(operation)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    if duplicate not in ('rm_old', 'rm_new', 'none'):
        print >> sys.stderr, "'{0:s}' is an invalid duplicate strategy.".format(duplicate)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    if motif2tf not in ('tf_motif_mapping', 'tf_motif_mapping_v6_1'):
        print >> sys.stderr, "'{0:s}' is an invalid motif2tf database.".format(motif2tf)
        print >> sys.stderr
        display_usage(sys.stderr)
        sys.exit(1)

    nomenclature2cfg = dict()
    nomenclature2p = dict()
    for nomenclature, filename in db_cfg.items('configurations'):
        if not os.path.isfile(filename):
            print >> sys.stderr, "'{0:s}' file doesn't exist.".format(filename)
            sys.exit(1)

        cfg = RawConfigParser()
        cfg.read(filename)
        nomenclature2cfg[nomenclature.upper()] = cfg

        with open(filename, 'r') as input_fh:
            rawdata = input_fh.read()

        p = Bunch(rank_threshold=cfg.getint('params', 'rank_threshold'),
                  auc_threshold=cfg.getfloat('params', 'auc_threshold'),
                  nes_threshold=cfg.getfloat('params', 'nes_threshold'),
                  motifcombinations=(
                      get_option(cfg, 'params', ('feature_combinations', 'motif_recombinations')).split(';')
                      if has_option(cfg, 'params', ('feature_combinations', 'motif_recombinations'))
                      else []),
                  inifile=rawdata)

        nomenclature2p[nomenclature.upper()] = p

    return db_cfg, gmt_filename, nomenclature2cfg, nomenclature2p, operation, duplicate, motif2tf


def generate_named_curves(id2gene_signatures, nomenclature2p, nomenclature2cfg):
    for name, signature in id2gene_signatures.iteritems():
        logging.info("{0:s}: Analyzing ...".format(name))
        p = nomenclature2p[signature.nomenclature]
        cfg = nomenclature2cfg[signature.nomenclature]
        try:
            curves = [RecoveryCurves(db, p.auc_threshold, p.rank_threshold, p.nes_threshold)
                      for db in get_databases(cfg, p, signature.identifiers)]

            logging.info("{0:s}: Fraction of genes/regions mappable to the database(s) = {1:.2f}%".format(
                name,
                min(curve.fraction_of_mapped_gene_ids(signature.identifiers) for curve in curves) * 100.0)
            )

            yield name, curves
        except AssertionError, msg:
            logging.warning("{0:s}: The following error occured '{1:s}'".format(name, msg))


def consume_queue_item(named_curves_queue, n_remaining_processes, name2database_id, name2signature, nomenclature2p):
    while True:
        # Blocking get ...
        name, curves = named_curves_queue.get()
        try:
            logging.error("{0:s}: Current items in the queue = {1:d}".format(name, named_curves_queue.qsize()))
        except NotImplementedError:
            pass

        if not curves:
            return

        p = nomenclature2p[name2signature[name].nomenclature]
        insert_analysis(current_process().connection, name2database_id[name], name2signature[name], p, curves)
        n_remaining_processes.value -= 1

        logging.info('Number of remainig signatures: {0:d}/{1:d}'.format(n_remaining_processes.value,
                                                                         current_process().n_signatures))


def produce_queue_item(name, signature):
    os.nice(produce_queue_item.niceness)
    try:
        for name, curves in generate_named_curves({name: signature},
                                                  produce_queue_item.nomenclature2p,
                                                  produce_queue_item.nomenclature2cfg):
            produce_queue_item.named_curves_queue.put((name, curves))
    except BaseException as e:
        # Make sure an item is put on the queue ...
        logging.error(e.message)
        dump_stacktrace2log()
        produce_queue_item.named_curves_queue.put((name, []))


def produce_queue_item_initializer(nomenclature2p, nomenclature2cfg, named_curves_queue, niceness):
    produce_queue_item.nomenclature2p = nomenclature2p
    produce_queue_item.nomenclature2cfg = nomenclature2cfg
    produce_queue_item.named_curves_queue = named_curves_queue
    produce_queue_item.niceness = niceness


def reset_database(connection, db_cfg, motif2tf):
    template_base_folder = os.path.dirname(DATABASE_DATAMODEL_TEMPLATE_FILE)
    namespace = {'operation': 'init', 'database_name': db_cfg.get('dbserver', 'database'), 'motif2tf': motif2tf}

    with open(DATABASE_DATAMODEL_TEMPLATE_FILE, 'r') as input_fh:
        template = airspeed.Template(input_fh.read())

    ddl = str(template.merge(namespace, loader=airspeed.CachingFileLoader(template_base_folder)))

    cursor = connection.cursor()
    cursor.execute(ddl)
    connection.commit()
    cursor.close()


def dump_stacktrace2log():
    exc_type, exc_value, exc_traceback = sys.exc_info()
    output = StringIO.StringIO()
    traceback.print_exception(exc_type, exc_value, exc_traceback, file=output)
    logging.error(output.getvalue())


def main():
    start_time = timeit.default_timer()

    #
    # Parsing of input parameters and assembling all metadata into a context.
    #

    db_cfg, gmt_filename, nomenclature2cfg, nomenclature2p, operation, duplicate_strategy, motif2tf = parse_arguments(
        sys.argv[1:])

    init_logging(db_cfg.get('params', 'log_filename') if db_cfg.has_option('params', 'log_filename') else None)

    #
    # Open connection to database.
    #

    if gmt_filename == 'init':
        connection = create_db_connection(db_cfg, True)
        reset_database(connection, db_cfg, motif2tf)
        connection.close()
        sys.exit()
    else:
        connection = create_db_connection(db_cfg)
        nomenclature2code = load_nomenclature2code(connection)

    #
    # Loading the gene signatures and metadata into memory.
    #

    with open(gmt_filename, 'r') if gmt_filename != 'stdin' else sys.stdin as input_fh:
        gene_signature_rawdata = input_fh.readlines()

    collection_metadata = MetadataEnrichedGeneSignature.load_collection_metadata(gene_signature_rawdata)
    if not collection_metadata.name:
        logging.error("No name for gene signature collection provided.")
        sys.exit(2)
    if operation == 'rm':
        remove_collection(connection, collection_metadata)
        connection.commit()
        connection.close()
        sys.exit(0)
    try:
        collection_id = insert_collection_metadata(connection, collection_metadata)
    except CollectionAlreadyPresentException as e:
        if operation in ('add', 'complete'):
            collection_id = e.collection_id
        else:
            logging.error("Collection with name '{0:s}' is already in the database.".format(collection_metadata.name))
            sys.exit(1)
    name2signature = dict()
    name2database_id = dict()
    for signature in MetadataEnrichedGeneSignature.load_signatures(gene_signature_rawdata):
        if signature.nomenclature not in nomenclature2code:
            logging.error("Gene signature with ID '{0:s}' has an unknown gene nomenclature ID: '{1:s}'.".format(
                signature.name,
                signature.nomenclature)
            )
        elif signature.name in name2signature:
            logging.error("ID '{0:s}' for gene signature is not unique.".format(signature.name))
        else:
            signature_id, new, analysis_present = insert_gene_signature(connection,
                                                                        collection_id,
                                                                        collection_metadata,
                                                                        signature,
                                                                        nomenclature2code,
                                                                        duplicate_strategy)

            if signature_id:  # Signatures can be skipped based on duplicate PubMedID ...
                if operation == 'complete':
                    if not analysis_present:
                        name2signature[signature.name] = signature
                        name2database_id[signature.name] = signature_id
                else:
                    name2signature[signature.name] = signature
                    name2database_id[signature.name] = signature_id
                if not new and operation not in ('add', 'complete'):
                    logging.warning(
                        "Gene signature with ID '{0:s}' is already present in database.".format(signature.name))
                elif new and operation in ('add', 'complete'):
                    logging.warning(
                        "Gene signature with ID '{0:s}' was not yet present in database.".format(signature.name))

    #
    # Doing recovery analysis.
    #

    n_cores = db_cfg.getint('params', 'n_cores') if db_cfg.has_option('params', 'n_cores') else 1
    niceness = db_cfg.getint('params', 'process_niceness') if db_cfg.has_option('params', 'process_niceness') else 19

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

        # To avoid memory overload the size of the queue is limited.
        named_curves_queue = Queue(maxsize=10)

        # More consumers than producers, because dumping in the database is the bottleneck step.
        n_producers = max(1, int(float(n_cores - 1) / 3.0))
        logging.info("Number of producers = {0:d}.".format(n_producers))
        producer_pool = Pool(processes=n_producers,
                             initializer=produce_queue_item_initializer,
                             initargs=(nomenclature2p, nomenclature2cfg, named_curves_queue, niceness))

        for name, gene_signature in name2signature.iteritems():
            producer_pool.apply_async(produce_queue_item, (name, gene_signature))

        producer_pool.close()

        consumer_pool = []
        n_consumers = max(1, n_cores - n_producers - 1)
        logging.info("Number of consumers = {0:d}.".format(n_consumers))
        n_remaining_processes = Value('i', len(name2signature))
        for idx in range(n_consumers):
            process = Process(target=consume_queue_item,
                              args=(named_curves_queue,
                                    n_remaining_processes,
                                    name2database_id,
                                    name2signature,
                                    nomenclature2p))
            process.n_signatures = len(name2signature)
            process.connection = create_db_connection(db_cfg)
            consumer_pool.append(process)
            process.start()

        producer_pool.join()
        for process in consumer_pool:
            process.terminate()
            process.connection.close()
    else:
        # To avoid memory overload when analyzing large datasets, a producer-consumer pattern is used. This
        # is implemented via a python generator function.
        n_signatures = len(name2signature)
        for idx, (name, curves) in enumerate(generate_named_curves(name2signature, nomenclature2p, nomenclature2cfg)):
            p = nomenclature2p[name2signature[name].nomenclature]
            logging.info("Progress: {0:d}/{1:d}".format(idx + 1, n_signatures))
            insert_analysis(connection, name2database_id[name], name2signature[name], p, curves)

    #
    # Size of database.
    #

    # Try to reconnect to the server if necessary (to avoid (2006, 'MySQL server has gone away')).
    connection.ping(True)

    cursor = connection.cursor()
    cursor.execute("""SELECT * FROM sizeOfDatabase;""")
    sizeInMb = cursor.fetchone()[0]
    cursor.close()
    print >> sys.stderr, "Total size of database = {0:f}Mb".format(sizeInMb)

    #
    # Close connection to database.
    #
    connection.close()

    logging.info("Elapsed time = {0:f}s".format(timeit.default_timer() - start_time))


if __name__ == "__main__":
    main()
