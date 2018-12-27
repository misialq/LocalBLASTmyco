import glob
import pprint
import shutil
import threading
from datetime import datetime, timedelta
import os
import time
import types
from tqdm import tqdm

import configobj
import pkg_resources
import validate
import logging
from multiprocessing.pool import Pool

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML

# BLAST installation instructions can be found here:
# https://www.blaststation.com/intl/members/en/howtoblastmac.html
# or even better: http://www.wernerlab.org/software/macqiime/macqiime-installation/installing-blast-in-os-x
from LocalBLAST import exceptions
from LocalBLAST.exceptions import BLASTDataMissing
from LocalBLAST.helpers import find_gene_name, find_uniID, generate_uni_link, get_uniprot_mass, get_uniprot_locus

tqdm.pandas()
logger = logging.getLogger()


class LocalSingleBLAST:
    def __init__(self, filename, config, xml_result_loc):
        self.filename = filename
        self.config = config
        #self.thread_id = str(threading.get_ident())
        self.thread_id = os.getpid()
        self.xml_result_path = os.path.join(xml_result_loc, f'BLAST_temp_result_{self.thread_id}.xml')
        self.blast_dbase = self.config['blastdb']
        self.sequence_len_threshold = self.config['sequence_len_threshold']

        self.seq_len = None
        self.seq_qc = False
        self.current_result = None
        self.final_result = None
        self.gene_id = None
        self.blast_title = None
        self.uniprot_id = None
        self.e_value = None

    def run_blast(self):
        handle = open(self.filename, "rU")
        if 'fas' in self.config['input_format'].lower():
            format_type = 'fasta'
        else:
            logger.error('Formats other than FASTA are currently not supported.')
            raise NotImplementedError
        # read the sequence from fas file and save it to 'sequence' variable
        for record in SeqIO.parse(handle, format_type):
            self.seq_len = len(record.seq)
        handle.close()

        # remove existing XML result file, if any
        if os.path.exists(self.xml_result_path):
            os.remove(self.xml_result_path)
            logger.debug(f'Removing temporary BLAST result ({self.xml_result_path})...')

        # perform a local BLAST, provided that sequence length is higher than given threshold
        if self.seq_len > self.sequence_len_threshold:
            if self.config['blast_engine'].lower() == 'blastx':
                from Bio.Blast.Applications import NcbiblastxCommandline
                comline = NcbiblastxCommandline(query=self.filename,
                                                db=self.blast_dbase,
                                                outfmt=5,
                                                out=self.xml_result_path)
            elif self.config['blast_engine'].lower() == 'blastp':
                from Bio.Blast.Applications import NcbiblastpCommandline
                comline = NcbiblastpCommandline(query=self.filename,
                                                db=self.blast_dbase,
                                                outfmt=5,
                                                out=self.xml_result_path)
            elif self.config['blast_engine'].lower() == 'blastn':
                from Bio.Blast.Applications import NcbiblastnCommandline
                comline = NcbiblastnCommandline(query=self.filename,
                                                db=self.blast_dbase,
                                                outfmt=5,
                                                out=self.xml_result_path)
            else:
                logger.error(f'BLASTing using {self.config["blast_engine"]} is not yet supported.')
                raise NotImplementedError
            os.system(str(comline))
            self.seq_qc = True

    def read_result(self):
        try:
            result = open(self.xml_result_path)
            records = NCBIXML.parse(result)
            blast_record = next(records)

            self.current_result = blast_record.descriptions[0]
            result.close()

        except ValueError:
            # ValueError arises when BLAST did not yield any result, i.e. XML file was empty
            if not self.seq_qc:
                logger.error(f'[{self.filename}] ValueError: Low sequencing quality')
            else:
                logger.error(f'[{self.filename}] ValueError')
        except FileNotFoundError:
            # FileNotFoundError arises when BLAST didn't yield a result because of too short query sequence
            logger.error(f'[{self.filename}] FileNotFoundError: BLAST result is missing')
        except IndexError:
            logger.error(f'[{self.filename}] IndexError')
        except UnboundLocalError:
            pass

        finally:
            if self.current_result:
                self.blast_title = self.current_result.title
                self.e_value = self.current_result.e
                self.gene_id = find_gene_name(self.blast_title)
                self.uniprot_id = find_uniID(self.blast_title)
            else:
                self.blast_title, self.e_value, self.gene_id, self.uniprot_id = np.nan, np.nan, np.nan, np.nan

            self.final_result = {
                'filename': self.filename,
                'seq_length': self.seq_len,
                'e_value': self.e_value,
                'gene_id': self.gene_id,
                'uniprot_id': self.uniprot_id
            }
            logger.debug(f'[{self.filename}] Result: {self.final_result}')


class LocalMultiBLAST:
    def __init__(self, config_file, use_uniprot, clean_up=False):
        self.timestamp_start = datetime.now()
        self.timestamp_end = None
        self.completion_time = None
        self.use_uniprot = use_uniprot
        self.clean_up = clean_up
        logger.info('Initializing BLAST object...')
        self.cwd = os.path.dirname(os.path.realpath(__file__))
        validator = validate.Validator()
        if config_file is None:
            self.config = {}
        self.config = configobj.ConfigObj(
            config_file,
            configspec=pkg_resources.resource_filename(
                __name__,
                '../config/LocalBLAST.spec'))
        self.validated = self.config.validate(validator)
        if isinstance(self.validated, dict) or not self.validated:
            logger.error(pprint.pformat(self.validated))
            raise exceptions.ConfigNotValid('Configuration file contains error')
        if self.cwd != "":
            os.chdir(self.cwd)

        self.blast_dbase = self.config['blastdb']
        logger.info(f'{self.blast_dbase.upper()} BLAST local database will be used.')

        self.sequence_len_threshold = self.config['sequence_len_threshold']
        logger.info(f'Only sequences longer than {self.sequence_len_threshold} will be considered.')

        # create a temporary location for blast results
        self.blast_tmp_dir = os.path.join(self.cwd, 'blast_temp')
        if not os.path.isdir(self.blast_tmp_dir):
            os.mkdir(self.blast_tmp_dir)
        self.val_err = False
        self.ind_err = False
        self.input_sequences_loc = self.config['sequences_loc']
        self.input_format = self.config['input_format'].lower()
        self.max_proc_count = self.config['max_proc_count']

        # other inits:
        self.filenames = list()
        self.ident_seqs = list()
        self.file_count = None
        self.all_files = None
        self.all_blast_results = None
        self.results_df = None
        self.results_loc = self.config['results_loc']
        if not os.path.isdir(self.results_loc):
            os.mkdir(self.results_loc)
        self.result_file_path = os.path.join(self.results_loc, 'blast_results_summary.csv')

    def gather_files(self):
        self.all_files = glob.glob(os.path.join(self.input_sequences_loc, f'*.{self.input_format}'))
        self.file_count = len(self.all_files)
        logger.info(f'Found {self.file_count} {self.input_format.upper()} files in {self.input_sequences_loc}')

    def process_file(self, filename):
        single_blast_obj = LocalSingleBLAST(filename, self.config, self.blast_tmp_dir)
        single_blast_obj.run_blast()
        single_blast_obj.read_result()
        return single_blast_obj.final_result

    def generate_df_report(self):
        if self.all_blast_results:
            self.results_df = pd.DataFrame(self.all_blast_results)
            logger.info('Generating UniProt links...')
            self.results_df['uniprot_link'] = self.results_df['uniprot_id'].progress_apply(generate_uni_link)
            self.results_df.fillna('N/A', inplace=True)
            self.results_df.sort_values('filename', ascending=True, inplace=True)
            if self.use_uniprot:
                logger.info('[UniProt] Fetching protein masses...')
                self.results_df['mass'] = self.results_df['uniprot_link'].progress_apply(get_uniprot_mass)
                logger.info('[UniProt] Fetching gene loci...')
                self.results_df['locus_tag'] = self.results_df['uniprot_link'].progress_apply(get_uniprot_locus)
            self.results_df.reset_index(inplace=True)

            # export result but check if file exists
            if os.path.isfile(self.result_file_path):
                current_ts = datetime.now().strftime('%Y%m%d%H%M%S')
                new_filename = f'blast_results_summary_{current_ts}.csv'
                logger.warning(f'Result file already exists in "{self.results_loc}". '
                               f'Existing file will be renamed to "{new_filename}"".')
                os.rename(os.path.join(self.results_loc, 'blast_results_summary.csv'),
                          os.path.join(self.results_loc, new_filename))
            self.results_df.to_csv(self.result_file_path)
        else:
            logger.error('No data to generate DataFrame from. Aborting.')
            raise BLASTDataMissing

    def run(self):
        self.gather_files()
        if self.file_count:
            pool = Pool(processes=self.max_proc_count)
            jobs = [pool.apply_async(self.process_file, (filename,)) for filename in self.all_files]
            while not all([job.ready() for job in jobs]):
                logger.info('BLASTing - waiting for job completion...')
                time.sleep(2)
            self.all_blast_results = [job.get() for job in jobs]
        self.generate_df_report()

        # clean up (remove temp BLAST results)
        if self.clean_up:
            logger.info('Cleaning up...')
            if os.path.isdir(self.blast_tmp_dir):
                shutil.rmtree(self.blast_tmp_dir)
                logger.info(f'Folder {self.blast_tmp_dir} and all its contents were removed.')
            else:
                logger.info(f'Folder {self.blast_tmp_dir} does not exist therefore nothing will be removed.')

        self.timestamp_end = datetime.now()
        self.completion_time = (self.timestamp_end - self.timestamp_start) / timedelta(minutes=1)
        logger.info('MultiBLAST completed in {:0.1f} minutes'.format(self.completion_time))

