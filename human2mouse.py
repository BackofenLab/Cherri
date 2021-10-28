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
    human_g_name2t_id = gtf_mapping(HUMAN_GTF)
    mouse_g_name2t_id = gtf_mapping(MOUSE_GTF)
    mapping = {}
    for x, df_slice in df.iterrows():
        try:
            human = df_slice["human"]
            mouse = df_slice["mouse"]
            human_t_ids = human_g_name2t_id[human]
            mouse_t_ids = mouse_g_name2t_id[mouse]
            for entry in human_t_ids:
                mapping[entry] = mouse_t_ids
            for entry in mouse_t_ids:
                mapping[entry] = human_t_ids
        except KeyError:
            pass

    with open(OUTFILE, "wb") as handle:
        pickle.dump(mapping, handle)


def gtf_mapping(file: str):
    gene_to_tanscripts = {}
    transcript_to_gene = {}
    with gzip.open(file, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            line = line.rstrip().split("\t")
            if line[2] == "transcript":
                gene_id = line[8].split('gene_name "')[-1].split('"')[0].split('.')[0]
                trans_id = line[8].split('transcript_id "')[-1].split('"')[0].split('.')[0]
                transcript_to_gene[trans_id] = gene_id
                if gene_id not in gene_to_tanscripts:
                    gene_to_tanscripts[gene_id] = set()
                gene_to_tanscripts[gene_id].add(trans_id)

    return gene_to_tanscripts




def main():
    mapping()


if __name__ == '__main__':
    main()