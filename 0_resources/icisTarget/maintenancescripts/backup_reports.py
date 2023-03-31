#!/usr/bin/env python
import base64
import os
import StringIO
import zipfile
import zlib

try: from ftplib import FTP_TSL as FTP
except ImportError: from ftplib import FTP


def upload_files(cfg, job_unique_id, filenames):
    ftp = FTP()
    ftp.connect(cfg.get('ftpserver', 'servername'), port=cfg.getint('ftpserver', 'port') if cfg.has_option('ftpserver', 'port') else 21)
    ftp.login(cfg.get('ftpserver', 'username'), base64.b64decode(cfg.get('ftpserver', 'password')))
    if cfg.get('ftpserver','folder'):
        ftp.cwd(os.path.expandvars(cfg.get('ftpserver','folder')))
    ftp.mkd(str(job_unique_id)); ftp.cwd(str(job_unique_id))
    for filename in filenames:
        with open(filename, 'rb') as input_fh:
            ftp.storbinary('STOR {0:s}'.format(os.path.basename(filename)), input_fh)
    ftp.quit()


def create_zip_output_handle(filename):
    return zipfile.ZipFile(filename, 'w', zipfile.ZIP_DEFLATED)


def add2zip(zip_output_fh, filename, data_as_string):
    data = StringIO.StringIO(data_as_string)
    info = zip_output_fh.zipinfo()
    info.name = filename
    info.size = data.len
    zip_output_fh.writestr(info, data)


def handle_download(block):
    file.write(block)

def archive_folder(session, foldername, output_foldername):
    """

    """
    session.cwd(foldername)

    with open(filename, 'w')
    session.retrbinary('RETR {0:s}'.format('archive.zip'), )
    output_fh = open(filename, 'w')
    def _handle_download(block): output_fh.write(block)
    for filename in session.nlst():
        if


    pass

def main():
    pass


if __name__ == "__main__":
    main()