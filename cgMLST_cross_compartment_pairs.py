import re
import numpy as np
import pandas as pd

ALLELE_FILE = "results_alleles.tsv"
METADATA_FILE = "All_type.csv"
ALL_PAIRS_OUT = "all_pairs_under_100_shared_only.csv"
CROSS_NICHE_OUT = "cross_niche_under_10_results.csv"
INITIAL_THRESHOLD = 100
TRANSMISSION_THRESHOLD = 10


def extract_numeric(val):
    if pd.isna(val):
        return -1
    match = re.search(r"(\d+)", str(val))
    return int(match.group(1)) if match else -1


def main():
    df_meta = pd.read_csv(METADATA_FILE)
    meta_map = {
        str(row.iloc[0]).strip(): str(row.iloc[1]).strip()
        for _, row in df_meta.iterrows()
    }

    try:
        df_alleles = pd.read_csv(ALLELE_FILE, sep="\t", index_col=0)
    except Exception:
        df_alleles = pd.read_csv(ALLELE_FILE, sep="\t")
        df_alleles.set_index(df_alleles.columns[0], inplace=True)

    raw_ids = df_alleles.index.astype(str).str.strip().tolist()
    matrix = df_alleles.map(extract_numeric).values.astype(np.int32)

    n = len(raw_ids)
    all_pairs = []

    for i in range(n):
        target_matrix = matrix[i + 1 :]
        if len(target_matrix) == 0:
            continue

        both_have_data = (matrix[i] != -1) & (target_matrix != -1)
        diff_mask = both_have_data & (matrix[i] != target_matrix)
        dists = diff_mask.sum(axis=1)

        hits = np.where(dists <= INITIAL_THRESHOLD)[0]
        for idx in hits:
            j = i + 1 + idx
            type_i = meta_map.get(raw_ids[i], "Unknown")
            type_j = meta_map.get(raw_ids[j], "Unknown")

            all_pairs.append(
                {
                    "Source": raw_ids[i],
                    "Target": raw_ids[j],
                    "Distance": int(dists[idx]),
                    "Type_S": type_i,
                    "Type_T": type_j,
                    "Category": "Same_Niche" if type_i == type_j else "Cross_Niche",
                }
            )

    if not all_pairs:
        return

    df_all = pd.DataFrame(all_pairs)
    df_all.to_csv(ALL_PAIRS_OUT, index=False, encoding="utf-8-sig")

    mask = (
        (df_all["Distance"] <= TRANSMISSION_THRESHOLD)
        & (df_all["Category"] == "Cross_Niche")
        & (df_all["Type_S"] != "Unknown")
        & (df_all["Type_T"] != "Unknown")
    )
    df_filtered = df_all[mask].copy()

    if not df_filtered.empty:
        df_filtered["Path"] = df_filtered.apply(
            lambda row: " <-> ".join(sorted([str(row["Type_S"]), str(row["Type_T"])])),
            axis=1,
        )
        df_filtered = df_filtered.sort_values(["Path", "Distance"])
        df_filtered.to_csv(CROSS_NICHE_OUT, index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()