import os
import pandas as pd
import pickle

HUMAN_POS = os.path.join("mapping_data", "paris_HEK293T_context_150_pos_occ_pos.csv")
MOUSE_POS = os.path.join("mapping_data", "paris_mouse_context_150_pos_occ_pos.csv")
MAPPING = os.path.join("mapping_data", "human_mouse_mapping.pckl")


def main():
    mapping = load_mapping()
    human_df = pd.read_csv(HUMAN_POS, sep=",")
    mouse_df = pd.read_csv(MOUSE_POS, sep=",")
    human_df["FID"] = human_df["ID_1st"].str.split(";").str[0]
    mouse_df["FID"] = mouse_df["ID_1st"].str.split(";").str[0]
    counter = 0
    forbidden = set()
    for x, df_slice in human_df.iterrows():
        try:
            h_t_id = df_slice["FID"]
            if h_t_id in forbidden:
                continue
            mouse_t_ids = mapping[h_t_id]
            for m_id in mouse_t_ids:
                sub_df = mouse_df[mouse_df["FID"] == m_id]
                if len(sub_df) > 0:
                    counter += 1
            forbidden.add(h_t_id)
        except KeyError:
            pass
    print(counter)


def load_mapping():
    with open(MAPPING, "rb") as handle:
        mapping = pickle.load(handle)
    return mapping

if __name__ == '__main__':
    main()