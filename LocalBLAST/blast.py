import glob
import pprint
import threading
from datetime import datetime, timedelta
import os
import time
import types

import configobj
import pkg_resources
import validate
import logging
from queue import Empty, Queue

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML

# BLAST installation instructions can be found here:
# https://www.blaststation.com/intl/members/en/howtoblastmac.html
# or even better: http://www.wernerlab.org/software/macqiime/macqiime-installation/installing-blast-in-os-x
from LocalBLAST import exceptions
from LocalBLAST.helpers import find_gene_name, find_uniID

logger = logging.getLogger()


class Worker(threading.Thread):
    """ Thread executing tasks from a given tasks queue """

    def __init__(self, tasks):
        threading.Thread.__init__(self)
        self.tasks = tasks
        self.stopping_signal = threading.Event()
        self.start()

    def run(self):
        while not self.stopping_signal.is_set():
            try:
                func, args, kwargs = self.tasks.get(block=False, timeout=1)
            except Empty:
                continue
            try:
                func(*args, **kwargs)
            except Exception as e:
                # An exception happened in this thread
                logger.error(e)
            finally:
                # Mark this task as done, whether an exception happened or not
                self.tasks.task_done()


class ThreadPool:
    """ Pool of threads consuming tasks from a queue """

    def __init__(self, num_threads):
        self.tasks = Queue(num_threads)
        self.workers = []
        for _ in range(num_threads):
            self.workers.append(Worker(self.tasks))

    def add_task(self, func, *args, **kargs):
        """ Add a task to the queue """
        self.tasks.put((func, args, kargs))

    def map(self, func, args_list):
        """ Add a list of tasks to the queue """
        for args in args_list:
            self.add_task(func, args)

    def close(self):
        for worker in self.workers:
            worker.stopping_signal.set()

    def wait_completion(self):
        """ Wait for completion of all the tasks in the queue """
        self.tasks.join()
        self.close()


class LocalSingleBLAST:
    def __init__(self, filename, config, xml_result_loc):
        self.filename = filename
        self.config = config
        self.thread_id = str(threading.get_ident())
        self.xml_result_path = os.path.join(xml_result_loc, 'BLAST_temp_result_', self.thread_id, '.xml')

        self.blast_dbase = self.config['blastdb']
        logger.info(f'{self.blast_dbase.upper()} BLAST local database will be used.')

        self.sequence_len_threshold = self.config['sequence_len_threshold']
        logger.info(f'Only sequences longer than {self.sequence_len_threshold} will be considered.')

        self.seq_len = None
        self.seq_qc = False

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

            current_result = blast_record.descriptions[0]
            title = current_result.title
            e_val = current_result.e
            gene = find_gene_name(title)
            uniID = find_uniID(title)

        except ValueError:
            # ValueError arises when BLAST did not yield any result, i.e. XML file was empty
            if seq_qc == False:
                val_err = True
                err_count += 1
                log.write(str(filx[0:-4] + " Low sequencing quality, err_count = " + str(err_count)) + "\n")
                print("ValueError at {0} file. ".format(filx[0:-4]))
                my_utils.append_value([ident_seqs, expects, genes, uniIDs], np.nan)
            else:
                val_err = True
                err_count += 1
                log.write(str(filx[0:-4] + " ValueError, err_count = " + str(err_count)) + "\n")
                print("ValueError at {0} file. ".format(filx[0:-4]))
                my_utils.append_value([ident_seqs, expects, genes, uniIDs], np.nan)
            pass

        except FileNotFoundError:
            # FileNotFoundError arises when BLAST didn't yield a result because of too short query sequence
            err_count += 1
            log.write(str(filx[0:-4] + " BLAST result missing, err_count = " + str(err_count)) + "\n")
            print("FileNotFoundError at {0} file. ".format(filx[0:-4]))
            my_utils.append_value([ident_seqs, expects, genes, uniIDs], np.nan)
            pass

        except IndexError:
            err_count += 1
            log.write(str(filx[0:-4] + " IndexError, err_count = " + str(err_count)) + "\n")
            print("IndexError at {0} file. ".format(filx[0:-4]))
            my_utils.append_value([ident_seqs, expects, genes, uniIDs], np.nan)
            pass

        finally:
            result.close()




class LocalMultiBLAST:
    def __init__(self, config_file):
        self.timestamp_start = datetime.now()
        self.timestamp_end = None
        self.completion_time = None
        logger.info('Initializing BLAST object...')
        self.cwd = os.path.dirname(os.path.realpath(__file__))
        validator = validate.Validator()
        if config_file is None:
            self.config = {}
        self.config = configobj.ConfigObj(
            config_file,
            configspec=pkg_resources.resource_filename(
                __name__,
                '../config/LocalMultiBLAST.spec'))
        self.validated = self.config.validate(validator)
        if isinstance(self.validated, dict) or not self.validated:
            logger.error(pprint.pformat(self.validated))
            raise exceptions.ConfigNotValid('Configuration file contains error')
        if self.cwd != "":
            os.chdir(self.cwd)

        # create a temporary location for blast results
        self.blast_tmp_dir = os.path.join(self.cwd, 'blast_temp')
        if not os.path.isdir(self.blast_tmp_dir):
            os.mkdir(self.blast_tmp_dir)
        self.val_err = False
        self.ind_err = False
        self.input_sequences_loc = self.config['sequences_loc']
        self.input_format = self.config['input_format'].lower()

        # other inits:
        self.all_info_dict = dict()
        self.filenames = list()
        self.seq_lens = list()
        self.seq_qcs = list()
        self.ident_seqs = list()
        self.uni_ids = list()
        self.expects = list()
        self.genes = list()
        self.file_count = None
        self.all_files = None

    def gather_files(self):
        self.all_files = glob.glob(os.path.join(self.input_sequences_loc, f'*.{self.input_format}'))
        self.file_count = len(self.all_files)
        logger.info(f'Found {self.file_count} {self.input_format.upper()} files in {self.input_sequences_loc}')


    def run(self):
        self.gather_files()

def main():
    # go through all FASTA files in a selected location, save a name and perform a local blast, then
    # gets a GI number and adds it to a list

    counter = 0
    tot_err_cnt = 0
    seq_qc = True

    for root, dirs, files in os.walk(seq_loc):
        for i, filx in enumerate(files):
            if filx.endswith(".fas"):
                f_name = filx[0:-4]
                filenames.append(f_name)

                val_err = False

                # opens the FAS file to get the sequence length

                handle = open(os.path.join(root, filx), "rU")
                for record in SeqIO.parse(handle,
                                          "fasta"):  # read the sequence from fas file and save it to 'sequence' variable
                    seq_len = len(record.seq)
                handle.close()
                seq_lens.append(seq_len)


                # removes existing XML result file, if any
                os.remove(xml_result_path) if os.path.exists(xml_result_path) else None

                # performs a local BLAST, provided that sequence length is higher than given threshold
                if seq_len > seq_thr:
                    comline = NcbiblastxCommandline(query=os.path.join(root, filx), db=dbase, outfmt=5,
                                                    out=xml_result_path)
                    os.system(str(comline))
                    seq_qc = True
                else:
                    seq_qc = False
                 #   expects.append(np.nan)

                seq_qcs.append(seq_qc)
                err_count = 0

                # opens the XML file with BLAST result

                try:
                    result = open(xml_result_path)
                    records = NCBIXML.parse(result)
                    blast_record = next(records)

                    current_result = blast_record.descriptions[0]
                    title = current_result.title
                    e_val = current_result.e
                    gene = my_utils.find_gene_name(title)
                    uniID = my_utils.find_uniID(title)

                    ident_seqs.append(title)
                    uniIDs.append(uniID)
                    expects.append(e_val)
                    genes.append(gene)

                except ValueError:  # ValueError arises when BLAST did not yield any result, i.e. XML file was empty
                    if seq_qc == False:
                        val_err = True
                        err_count += 1
                        log.write(str(filx[0:-4] + " Low sequencing quality, err_count = " + str(err_count)) + "\n")
                        print("ValueError at {0} file. ".format(filx[0:-4]))
                        my_utils.append_value([ident_seqs, expects, genes, uniIDs], np.nan)
                    else:
                        val_err = True
                        err_count += 1
                        log.write(str(filx[0:-4] + " ValueError, err_count = " + str(err_count)) + "\n")
                        print("ValueError at {0} file. ".format(filx[0:-4]))
                        my_utils.append_value([ident_seqs, expects, genes, uniIDs], np.nan)
                    pass

                except FileNotFoundError:  # FileNotFoundError arises when BLAST didn't yield a result because of too short query sequence
                    err_count += 1
                    log.write(str(filx[0:-4] + " BLAST result missing, err_count = " + str(err_count)) + "\n")
                    print("FileNotFoundError at {0} file. ".format(filx[0:-4]))
                    my_utils.append_value([ident_seqs, expects, genes, uniIDs], np.nan)
                    pass

                except IndexError:
                    err_count += 1
                    log.write(str(filx[0:-4] + " IndexError, err_count = " + str(err_count)) + "\n")
                    print("IndexError at {0} file. ".format(filx[0:-4]))
                    my_utils.append_value([ident_seqs, expects, genes, uniIDs], np.nan)
                    pass

                assert len(uniIDs) == len(ident_seqs) == len(expects) == len(seq_lens), "Arrays are of unequal length"

                if i % 50 == 0:
                    current_progress = i * 100 / len(files)
                    print("Current progress: {0:1.1f}%".format(current_progress))

    results_all = pd.DataFrame({'Filename': filenames,
                                'Seq_len': seq_lens,
                                'UniProtID': uniIDs,
                                'E-value': expects,
                                'Quality': my_utils.map_quality(seq_lens, expects),
                                'Gene': genes,
                                'UniProt_link': [my_utils.generate_uni_link(x) for x in uniIDs]})

    results_all.fillna('N/A', inplace=True)
    results_all.sort_values('Filename', ascending=True, inplace=True)
    results_all.to_csv("../result/results_summary.csv")

    results_all['Mass'] = results_all['UniProt_link'].apply(my_utils.get_uniprot_mass)
    results_all['LocusTag'] = results_all['UniProt_link'].apply(my_utils.get_uniprot_locus)
    results_all.to_csv("../result/results_summary.csv")
