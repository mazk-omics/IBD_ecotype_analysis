#!/usr/bin/env python3

import gzip
import json
import math
import shlex
from pathlib import Path

import numpy as np


ROOT = Path("/home/mazekai/IBD_EcoTyper")

DATA = (
    ROOT
    / "01_raw_data"
    / "bulk"
    / "MSCCR_GSE193677"
)

RAW = DATA / "GSE193677_MSCCR_Biopsy_counts.txt.gz"
ADJ = DATA / "GSE193677_MSCCR_Biopsy_adjcounts.txt.gz"

OUT = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "08_msccr_bulk_audit"
)

OUT.mkdir(parents=True, exist_ok=True)

JSON_OUT = OUT / "MSCCR_bulk_schema_audit.json"
TXT_OUT = OUT / "MSCCR_bulk_schema_audit.txt"


def strip_quotes(x):
    if len(x) >= 2 and x[0] == '"' and x[-1] == '"':
        return x[1:-1]
    return x


def audit_matrix(path):

    print()
    print("=" * 72)
    print("AUDITING:", path.name)
    print("=" * 72)

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="strict"
    ) as fh:

        # ----------------------------------------------------
        # Header:
        # "sample1" "sample2" ... "sample2490"
        #
        # There is NO gene-column header token.
        # ----------------------------------------------------

        header_line = fh.readline().strip()

        if not header_line:
            raise RuntimeError("Empty header")

        sample_ids = shlex.split(header_line)

        n_samples = len(sample_ids)

        print("Header samples:", n_samples)

        sample_duplicate_count = (
            n_samples - len(set(sample_ids))
        )

        sample_blank_count = sum(
            x == "" for x in sample_ids
        )

        gene_ids = []
        gene_seen = set()

        duplicate_gene_count = 0
        blank_gene_count = 0

        valid_rows = 0
        bad_width_rows = 0
        nonnumeric_rows = 0

        total_values = 0
        negative_values = 0
        noninteger_values = 0
        zero_values = 0

        min_value = math.inf
        max_value = -math.inf

        sample_sums = np.zeros(
            n_samples,
            dtype=np.float64
        )

        first_genes = []

        for line_number, line in enumerate(
            fh,
            start=2
        ):

            line = line.strip()

            if not line:
                continue

            # Split ONCE:
            #
            # "ENSG..." <numeric payload>
            #
            pieces = line.split(
                maxsplit=1
            )

            if len(pieces) != 2:
                bad_width_rows += 1
                continue

            gene_id = strip_quotes(
                pieces[0]
            )

            numeric_text = pieces[1]

            values = np.fromstring(
                numeric_text,
                sep=" ",
                dtype=np.float64
            )

            if values.size != n_samples:
                bad_width_rows += 1
                continue

            if not np.all(
                np.isfinite(values)
            ):
                nonnumeric_rows += 1
                continue

            if gene_id == "":
                blank_gene_count += 1

            if gene_id in gene_seen:
                duplicate_gene_count += 1
            else:
                gene_seen.add(gene_id)

            gene_ids.append(gene_id)

            if len(first_genes) < 10:
                first_genes.append(gene_id)

            valid_rows += 1
            total_values += values.size

            negative_values += int(
                np.count_nonzero(
                    values < 0
                )
            )

            noninteger_values += int(
                np.count_nonzero(
                    np.abs(
                        values -
                        np.rint(values)
                    ) > 1e-8
                )
            )

            zero_values += int(
                np.count_nonzero(
                    values == 0
                )
            )

            min_value = min(
                min_value,
                float(values.min())
            )

            max_value = max(
                max_value,
                float(values.max())
            )

            sample_sums += values

            if valid_rows % 2500 == 0:
                print(
                    f"processed genes: {valid_rows:,}",
                    flush=True
                )

    result = {

        "file":
            str(path),

        "n_gene_rows":
            valid_rows,

        "n_samples":
            n_samples,

        "duplicate_gene_count":
            duplicate_gene_count,

        "blank_gene_count":
            blank_gene_count,

        "sample_duplicate_count":
            sample_duplicate_count,

        "sample_blank_count":
            sample_blank_count,

        "bad_width_rows":
            bad_width_rows,

        "nonnumeric_rows":
            nonnumeric_rows,

        "total_values":
            total_values,

        "negative_values":
            negative_values,

        "noninteger_values":
            noninteger_values,

        "noninteger_fraction":
            (
                noninteger_values
                / total_values
                if total_values
                else None
            ),

        "zero_values":
            zero_values,

        "zero_fraction":
            (
                zero_values
                / total_values
                if total_values
                else None
            ),

        "min_value":
            min_value,

        "max_value":
            max_value,

        "first_10_genes":
            first_genes,

        "first_10_samples":
            sample_ids[:10],

        "sample_sum_min":
            float(sample_sums.min()),

        "sample_sum_median":
            float(np.median(sample_sums)),

        "sample_sum_mean":
            float(sample_sums.mean()),

        "sample_sum_max":
            float(sample_sums.max()),

        "zero_sum_samples":
            int(
                np.count_nonzero(
                    sample_sums == 0
                )
            ),

        "_sample_ids":
            sample_ids,

        "_gene_ids":
            gene_ids
    }

    print()
    print("Genes              :", valid_rows)
    print("Samples            :", n_samples)
    print("Duplicate genes    :", duplicate_gene_count)
    print("Duplicate samples  :", sample_duplicate_count)
    print("Bad-width rows     :", bad_width_rows)
    print("Non-numeric rows   :", nonnumeric_rows)
    print("Negative values    :", negative_values)
    print("Non-integer values :", noninteger_values)

    print(
        "Non-integer frac   :",
        result["noninteger_fraction"]
    )

    print("Min                :", min_value)
    print("Max                :", max_value)

    print(
        "Zero fraction      :",
        result["zero_fraction"]
    )

    print(
        "First genes        :",
        first_genes
    )

    # Hard structural requirements.
    structure_pass = all([
        valid_rows > 0,
        n_samples > 0,
        duplicate_gene_count == 0,
        sample_duplicate_count == 0,
        blank_gene_count == 0,
        sample_blank_count == 0,
        bad_width_rows == 0,
        nonnumeric_rows == 0,
        negative_values == 0,
    ])

    print(
        "NONNEGATIVE COUNT-MATRIX STRUCTURE:",
        "PASS"
        if structure_pass
        else "FAIL"
    )

    return result


def main():

    raw = audit_matrix(RAW)
    adj = audit_matrix(ADJ)

    same_samples = (
        raw["_sample_ids"]
        == adj["_sample_ids"]
    )

    same_genes = (
        raw["_gene_ids"]
        == adj["_gene_ids"]
    )

    comparison = {
        "same_dimensions":
            (
                raw["n_gene_rows"]
                == adj["n_gene_rows"]
                and
                raw["n_samples"]
                == adj["n_samples"]
            ),

        "same_sample_ids_same_order":
            same_samples,

        "same_gene_ids_same_order":
            same_genes,

        "raw_all_nonnegative":
            raw["negative_values"] == 0,

        "raw_has_fractional_values":
            raw["noninteger_values"] > 0,

        "adj_has_fractional_values":
            adj["noninteger_values"] > 0,

        "adj_has_negative_values":
            adj["negative_values"] > 0,
    }

    raw.pop("_sample_ids")
    raw.pop("_gene_ids")

    adj.pop("_sample_ids")
    adj.pop("_gene_ids")

    output = {
        "raw_counts": raw,
        "adjusted_counts": adj,
        "cross_file_comparison":
            comparison
    }

    with open(
        JSON_OUT,
        "w",
        encoding="utf-8"
    ) as fh:
        json.dump(
            output,
            fh,
            indent=2,
            ensure_ascii=False
        )

    with open(
        TXT_OUT,
        "w",
        encoding="utf-8"
    ) as fh:

        fh.write(
            "MSCCR GSE193677 BULK AUDIT\n"
        )

        fh.write("=" * 60 + "\n\n")

        for label, x in (
            ("RAW COUNTS", raw),
            ("ADJCOUNTS", adj),
        ):

            fh.write(
                f"[{label}]\n"
            )

            for k, v in x.items():

                if isinstance(
                    v,
                    (str, int, float, bool)
                ) or v is None:

                    fh.write(
                        f"{k}: {v}\n"
                    )

            fh.write("\n")

        fh.write(
            "[CROSS-FILE COMPARISON]\n"
        )

        for k, v in comparison.items():

            fh.write(
                f"{k}: {v}\n"
            )

    print()
    print("=" * 72)
    print("CROSS-FILE COMPARISON")
    print("=" * 72)

    for k, v in comparison.items():
        print(f"{k}: {v}")

    # Hard requirements for the candidate
    # BayesPrism mixture.
    #
    # NOTE:
    # Fractional count values are recorded
    # but are NOT rounded and NOT treated
    # as an automatic failure.

    assert raw["n_gene_rows"] > 0
    assert raw["n_samples"] == 2490
    assert raw["duplicate_gene_count"] == 0
    assert raw["sample_duplicate_count"] == 0
    assert raw["bad_width_rows"] == 0
    assert raw["nonnumeric_rows"] == 0
    assert raw["negative_values"] == 0

    print()
    print("=" * 72)
    print(
        "MSCCR COUNT-MATRIX TECHNICAL AUDIT: PASS"
    )
    print("=" * 72)


if __name__ == "__main__":
    main()
