import numpy as np
from urllib.request import urlopen

def find_uniID(title):
    title_parts = title.split('|')
    uniID = title_parts[3]
    return uniID

def find_gene_name(title):
    title_parts = title.split('|')
    main_name = title_parts[-1]
    main_name_split = main_name.split()
    for text in main_name_split:
        if text[:3] == "GN=":
            return text[3:]

def generate_uni_link(id):
    if id is not np.nan:
        return 'http://www.uniprot.org/uniprot/' +id + '.txt'
    else:
        return 'N/A'

def append_value(list_of_lists, value):
    return [l.append(value) for l in list_of_lists]

def map_quality(seq_len_list, eval_list):
    assert len(seq_len_list) == len(eval_list), "Lists of uneqal length"
    qualities = []
    for l, e in zip(seq_len_list, eval_list):
        if e < 0.001 and l > 50:
            qualities.append('ok')
        elif l == 0:
            qualities.append(np.nan)
        else:
            qualities.append('curate')
    return qualities

def split_plate_name(plate_name):
    split_name = plate_name.split('-')
    return split_name


def get_uniprot_locus(uni_link):
    if "http" in uni_link:
        try:
            lines = urlopen(uni_link).readlines()
            for line in lines:
                line = str(line).lstrip("b'")
                if "OrderedLocusNames" in line:
                    split1 = line.split(';')
                    for element in split1:
                        if "OrderedLocusNames" in element:
                            split2 = element.split()
                            if len(split2) == 1:
                                locusName = split2[0][18:25].strip('; ')
                            else:
                                locusName = split2[1][18:25].strip('; ')

                    break
                else:
                    locusName = "N/A"
            return locusName
        except Exception as e:
            print(e)
            return "N/A"
    else:
        return "N/A"


def get_uniprot_mass(uni_link):
    if "http" in uni_link:
        try:
            lines = urlopen(uni_link).readlines()
            for line in lines:
                line = str(line).lstrip("b'")
                if 'SQ' in line[:2]:
                    mass = int(line.split(';')[1].strip()[:-3]) / 1000
                    break
                else:
                    mass = "N/A"
            return mass
        except Exception as e:
            print(e)
            return "N/A"
    else:
        return "N/A"