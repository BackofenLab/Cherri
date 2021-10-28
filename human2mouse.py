import pandas as pd
import os
import gzip
import pickle

OUTFILE = os.path.join("mapping_data", "human_mouse_mapping.pckl")


HUMAN2MOUSE_FILE = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                               "mapping_data", "HMD_HumanPhenotype.rpt")
DATA = "/home/rabsch/Documents/denbi_mount/rabsch/DATA/"

HUMAN_GTF = os.path.join(DATA, "gencode.v38.annotation.gtf.gz")
assert os.path.exists(HUMAN_GTF)
MOUSE_GTF = os.path.join(DATA, "gencode.vM7.annotation.gtf.gz")
assert os.path.exists(MOUSE_GTF)


def mapping():
    df = pd.read_csv(
        HUMAN2MOUSE_FILE,
        sep="\t",
        names=["human", "bla", "mouse", "bla2"], index_col=False)
    human_g_name2g_id = gtf_mapping(HUMAN_GTF)
    mouse_g_name2g_id = gtf_mapping(MOUSE_GTF)
    mapping = {}
    for x, df_slice in df.iterrows():
        try:
            human = human_g_name2g_id[df_slice["human"]]
            mouse = mouse_g_name2g_id[df_slice["mouse"]]
            mapping[mouse] = human
            mapping[human] = mouse
        except KeyError:
            pass
    with open(OUTFILE, "wb") as handle:
        pickle.dump(mapping, handle)


def gtf_mapping(file: str):
    data = {}
    with gzip.open(file, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            line = line.rstrip().split("\t")
            if line[2] == "gene":
                gene_id = line[8].split('gene_id "')[-1].split('"')[0].split('.')[0]
                gene_name = line[8].split('gene_name "')[-1].split('"')[0]
                data[gene_name] = gene_id
    return data


def main():
    mapping()


if __name__ == '__main__':
    main()