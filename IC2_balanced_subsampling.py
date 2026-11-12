import os
import numpy as np
import pandas as pd

INPUT_CSV = "IC2_all.csv"
OUTPUT_CSV = "IC2_uniform_year_balanced.csv"
TARGET_CLINICAL_TOTAL = 70
RANDOM_SEED = 42


def main():
    if not os.path.exists(INPUT_CSV):
        return

    df = pd.read_csv(INPUT_CSV)
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
    df["Country"] = df["Country"].fillna("Unknown").astype(str).str.strip()
    df["Type"] = (
        df["Type"].fillna("Clinical").astype(str).str.strip().str.capitalize()
    )

    df = df[df["Year"].notna()].copy()
    df["Year"] = df["Year"].astype(int)

    kept_indices = []

    non_clinical = df[df["Type"].isin(["Environmental", "Animal"])]
    kept_indices.extend(non_clinical.index)

    clinical_df = df[df["Type"] == "Clinical"].copy()
    all_years = sorted(clinical_df["Year"].unique())

    year_counts = clinical_df["Year"].value_counts().to_dict()
    remaining_quota = TARGET_CLINICAL_TOTAL
    active_years = list(all_years)
    year_target_quotas = {}

    while active_years and remaining_quota > 0:
        base_quota = remaining_quota // len(active_years)
        remainder = remaining_quota % len(active_years)

        exhausted = False
        for i, y in enumerate(active_years):
            allocated = base_quota + (1 if i < remainder else 0)
            if year_counts[y] <= allocated:
                year_target_quotas[y] = year_counts[y]
                remaining_quota -= year_counts[y]
                active_years.remove(y)
                exhausted = True
                break

        if not exhausted:
            for i, y in enumerate(active_years):
                allocated = base_quota + (1 if i < remainder else 0)
                year_target_quotas[y] = allocated
                remaining_quota -= allocated
            break

    np.random.seed(RANDOM_SEED)
    sampled_clin_list = []

    for y in all_years:
        quota = year_target_quotas.get(y, 0)
        if quota == 0:
            continue

        y_group = clinical_df[clinical_df["Year"] == y]
        if len(y_group) <= quota:
            sampled_clin_list.append(y_group)
        else:
            countries = [
                c for c in y_group["Country"].unique() if c != "Unknown"
            ]
            num_c = len(countries) if len(countries) > 0 else 1
            per_c_quota = int(np.ceil(quota / num_c))

            c_sampled = []
            for _, c_group in y_group.groupby("Country"):
                if len(c_group) <= per_c_quota:
                    c_sampled.append(c_group)
                else:
                    c_sampled.append(
                        c_group.sample(n=per_c_quota, random_state=RANDOM_SEED)
                    )

            y_sampled = pd.concat(c_sampled)
            if len(y_sampled) > quota:
                y_sampled = y_sampled.sample(n=quota, random_state=RANDOM_SEED)
            sampled_clin_list.append(y_sampled)

    sampled_clinical_df = pd.concat(sampled_clin_list)
    kept_indices.extend(sampled_clinical_df.index)

    final_indices = list(set(kept_indices))
    subsampled_df = (
        df.loc[final_indices]
        .sort_values(by=["Year", "ID"])
        .reset_index(drop=True)
    )
    subsampled_df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()
