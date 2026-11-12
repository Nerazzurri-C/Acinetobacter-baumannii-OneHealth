import csv
import os
import sys
import pandas as pd

META_FILE = "Final_Metadata_Merged.csv"


def filter_multi_thresholds_snp(input_file):
    sep = "\t" if input_file.endswith((".tsv", ".txt")) else ","

    out_files = {
        10: "pairwise_snp_le10.csv",
        20: "pairwise_snp_le20.csv",
        40: "pairwise_snp_le40.csv",
        100: "pairwise_snp_le100.csv",
    }

    handles = {k: open(v, "w", newline="", encoding="utf-8") for k, v in out_files.items()}
    writers = {k: csv.writer(handles[k]) for k in out_files}
    header = ["Strain_1", "Strain_2", "SNP_Distance"]
    for w in writers.values():
        w.writerow(header)

    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split(sep)
            if len(parts) < 3:
                continue

            s_a, s_b, dist_str = parts[0].strip(), parts[1].strip(), parts[2].strip()
            if s_a == s_b or s_a > s_b:
                continue

            try:
                snp = int(float(dist_str))
            except ValueError:
                continue

            if snp <= 100:
                writers[100].writerow([s_a, s_b, snp])
                if snp <= 40:
                    writers[40].writerow([s_a, s_b, snp])
                    if snp <= 20:
                        writers[20].writerow([s_a, s_b, snp])
                        if snp <= 10:
                            writers[10].writerow([s_a, s_b, snp])

    for h in handles.values():
        h.close()


def process_cross_type(snp_file, out_file, meta_dict):
    if not os.path.exists(snp_file):
        return

    df_snp = pd.read_csv(snp_file)
    df_snp["Strain_1"] = df_snp["Strain_1"].astype(str).str.strip()
    df_snp["Strain_2"] = df_snp["Strain_2"].astype(str).str.strip()

    records = []
    for row in df_snp.itertuples(index=False):
        s1 = row.Strain_1
        s2 = row.Strain_2
        dist = row.SNP_Distance

        info1 = meta_dict.get(s1, {})
        info2 = meta_dict.get(s2, {})

        type1 = str(info1.get("Type", "Unknown")).strip()
        type2 = str(info2.get("Type", "Unknown")).strip()

        if type1 != "Unknown" and type2 != "Unknown" and type1 != type2:
            records.append(
                {
                    "Strain_1": s1,
                    "Type_1": type1,
                    "Source_1": info1.get("Source", "Unknown"),
                    "Country_1": info1.get("Country", "Unknown"),
                    "MLST_1": info1.get("MLST", "Unknown"),
                    "Year_1": info1.get("Year", "Unknown"),
                    "Strain_2": s2,
                    "Type_2": type2,
                    "Source_2": info2.get("Source", "Unknown"),
                    "Country_2": info2.get("Country", "Unknown"),
                    "MLST_2": info2.get("MLST", "Unknown"),
                    "Year_2": info2.get("Year", "Unknown"),
                    "SNP_Distance": dist,
                }
            )

    res_df = pd.DataFrame(records)
    if not res_df.empty:
        res_df = res_df.sort_values(by=["SNP_Distance", "Strain_1", "Strain_2"])
        res_df.to_csv(out_file, index=False, encoding="utf-8-sig")


def main():
    target_file = sys.argv[1] if len(sys.argv) >= 2 else "snp_matrix_square.tsv"
    if not os.path.exists(target_file) or not os.path.exists(META_FILE):
        sys.exit(1)

    filter_multi_thresholds_snp(target_file)

    meta_df = pd.read_csv(META_FILE)
    meta_df["GCA_ID"] = meta_df["GCA_ID"].astype(str).str.strip()
    meta_dict = meta_df.set_index("GCA_ID").to_dict(orient="index")

    for threshold in [10, 20, 40, 100]:
        process_cross_type(
            f"pairwise_snp_le{threshold}.csv",
            f"cross_type_transmission_le{threshold}.csv",
            meta_dict,
        )


if __name__ == "__main__":
    main()