import os
import numpy as np
import pandas as pd

ANI_THRESHOLD = 99.95
MIN_COV_THRESHOLD = 90.0
FASTANI_FILE = "fastani_output.txt"
META_FILE = "All_type.csv"

meta_df = pd.read_csv(META_FILE)
meta_df["GCA_ID"] = meta_df["GCA_ID"].astype(str).str.strip()
type_dict = dict(zip(meta_df["GCA_ID"], meta_df["Type"]))

ani_df = pd.read_csv(
    FASTANI_FILE,
    sep="\t",
    header=None,
    names=["Query", "Ref", "ANI", "Mapped", "Total"],
)


def extract_id(path):
    filename = os.path.basename(path)
    return filename[:-4] if filename.endswith(".fna") else filename


ani_df["Query_ID"] = ani_df["Query"].apply(extract_id)
ani_df["Ref_ID"] = ani_df["Ref"].apply(extract_id)

ani_df = ani_df[ani_df["Query_ID"] != ani_df["Ref_ID"]].copy()
ani_df["Coverage"] = (ani_df["Mapped"] / ani_df["Total"]) * 100

ani_df["Strain_A"] = np.where(
    ani_df["Query_ID"] < ani_df["Ref_ID"],
    ani_df["Query_ID"],
    ani_df["Ref_ID"],
)
ani_df["Strain_B"] = np.where(
    ani_df["Query_ID"] < ani_df["Ref_ID"],
    ani_df["Ref_ID"],
    ani_df["Query_ID"],
)

pair_stats = (
    ani_df.groupby(["Strain_A", "Strain_B"])
    .agg(
        Mean_ANI=("ANI", "mean"),
        Min_Coverage=("Coverage", "min"),
        Mean_Coverage=("Coverage", "mean"),
        Direction_Count=("ANI", "count"),
    )
    .reset_index()
)

qualified_unique = pair_stats[
    (pair_stats["Mean_ANI"] >= ANI_THRESHOLD)
    & (pair_stats["Min_Coverage"] >= MIN_COV_THRESHOLD)
    & (pair_stats["Direction_Count"] == 2)
].copy()

qualified_unique["Type_A"] = qualified_unique["Strain_A"].map(type_dict)
qualified_unique["Type_B"] = qualified_unique["Strain_B"].map(type_dict)
unique_pairs = qualified_unique.dropna(subset=["Type_A", "Type_B"]).copy()
unique_pairs["Is_Cross_Interface"] = (
    unique_pairs["Type_A"] != unique_pairs["Type_B"]
)

unique_pairs.to_csv(
    "unique_strain_sharing_pairs.csv", index=False, encoding="utf-8-sig"
)

cross_only = unique_pairs[unique_pairs["Is_Cross_Interface"]].copy()
cross_only.to_csv(
    "cross_interface_only_pairs.csv", index=False, encoding="utf-8-sig"
)

valid_pair_tuples = set(zip(unique_pairs["Strain_A"], unique_pairs["Strain_B"]))
ani_df["Pair_Tuple"] = list(zip(ani_df["Strain_A"], ani_df["Strain_B"]))
directional_raw = ani_df[ani_df["Pair_Tuple"].isin(valid_pair_tuples)].copy()
directional_raw["Query_Type"] = directional_raw["Query_ID"].map(type_dict)
directional_raw["Ref_Type"] = directional_raw["Ref_ID"].map(type_dict)
directional_raw["Is_Cross_Interface"] = (
    directional_raw["Query_Type"] != directional_raw["Ref_Type"]
)
directional_raw = directional_raw.drop(
    columns=["Strain_A", "Strain_B", "Pair_Tuple"]
)
directional_raw.to_csv(
    "cross_interface_sharing_directional_raw.csv",
    index=False,
    encoding="utf-8-sig",
)

all_types = sorted(list(set(meta_df["Type"].dropna().unique())))
df_sym = pd.concat(
    [
        unique_pairs[["Type_A", "Type_B"]].rename(
            columns={"Type_A": "Row", "Type_B": "Col"}
        ),
        unique_pairs[["Type_B", "Type_A"]].rename(
            columns={"Type_B": "Row", "Type_A": "Col"}
        ),
    ],
    ignore_index=True,
)

if not df_sym.empty:
    matrix_symmetric = pd.crosstab(df_sym["Row"], df_sym["Col"]).reindex(
        index=all_types, columns=all_types, fill_value=0
    )
    for t in all_types:
        matrix_symmetric.loc[t, t] = int(matrix_symmetric.loc[t, t] // 2)
else:
    matrix_symmetric = pd.DataFrame(0, index=all_types, columns=all_types)

matrix_symmetric.to_csv(
    "interface_sharing_matrix_symmetric.csv", encoding="utf-8-sig"
)