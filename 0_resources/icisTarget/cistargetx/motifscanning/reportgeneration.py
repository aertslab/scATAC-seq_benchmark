import os
import itertools

import airspeed

from cistargetx.recoveryanalysis.reportgeneration import create_archive


BED_FILENAME = 'crms.bed'
HTML_REPORT_FILENAME = 'report.html'


def write_template(template_filename, output_filename, namespace):
    with open(output_filename, 'w') as output_fh:
        write_template2file_handle(template_filename, output_fh, namespace)


def write_template2file_handle(template_filename, output_fh, namespace):
    template_base_folder = os.path.dirname(template_filename)

    with open(template_filename, 'r') as input_fh:
        template = airspeed.Template(input_fh.read())

    output_fh.write(str(template.merge(namespace, loader=airspeed.CachingFileLoader(template_base_folder))))


def generate_report(output_folder, feature_iterator, template_filename, template_variables):
    i1, i2 = itertools.tee(feature_iterator)

    with open(os.path.join(output_folder, BED_FILENAME), 'w') as output_fh:
        for feature in i1:
            print >> output_fh, str(feature)

    namespace = {
        'bedfilename': BED_FILENAME,
        'features': list(i2)
    }

    namespace = dict(namespace.items() + template_variables.items())

    write_template(template_filename, os.path.join(output_folder, HTML_REPORT_FILENAME), namespace)

    # Create archive.
    create_archive(output_folder)
