from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
import os, time, types
import pandas as pd
import my_utils
import numpy as np

# BLAST installation instructions can be found here:
# https://www.blaststation.com/intl/members/en/howtoblastmac.html
# or even better: http://www.wernerlab.org/software/macqiime/macqiime-installation/installing-blast-in-os-x

def main():
    file_path = os.path.dirname(__file__)
    if file_path != "":
        os.chdir(file_path)

    # create a log file
    with open("../result/error_log.txt", 'w') as log:

        # create all necessary lists and variables

        filenames = []
        locus_tag = ""
        locus_tags = []
        GIs = []
        seq_lens = []
        seq_qcs = []
        ident_seqs = []
        uniIDs = []
        expects = []
        genes = []
        plates = []
        backbones = []
        primers = []
        wells = []
        val_err = False
        ind_err = False
        xml_result_path = "../blast/BLAST_result_temp.xml"

        # enter BLAST details below:

        org = "Mycobacterium tuberculosis H37Rv"  # organism
        seq_loc = str("../sequences")  # path to where the sequences are located
        que = str(input("\nUse the default BLAST database (myco_proteins_uniprot)? Y/N: "))

        if que.upper() == "Y" or que.upper() == "":
            dbase = "../blast/db/myco_proteins_uniprot"
        else:
            dbase = str(input("Database to BLAST against: "))  # database to blast against

        seq_thr = input("Set the sequence length threshold to [bp]: (default: 100) ")

        if seq_thr == "":
            seq_thr = 100
        else:
            seq_thr = int(seq_thr)

        # go through all FASTA files in a selected location, save a name and perform a local blast, then
        # gets a GI number and adds it to a list

        timer1 = time.time()

        file_no = len([name for name in os.listdir(seq_loc) if name.endswith(".fas")])
        counter = 0
        tot_err_cnt = 0
        seq_qc = True

        print("\nDatabase used for BLAST: " + str(dbase.upper()))

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


if __name__ == '__main__':
    main()