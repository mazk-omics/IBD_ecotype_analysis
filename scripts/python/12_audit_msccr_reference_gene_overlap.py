#!/usr/bin/env python3

"""
Audit gene-space compatibility between:

1. Unified scRNA reference
   14,431 gene symbols

2. GSE282122 gene metadata
   used as the audited gene_symbol <-> Ensembl bridge

3. MSCCR / GSE193677 raw bulk count matrix
   Ensembl gene IDs

This script does NOT modify or subset expression matrices.
"""

import gzip
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


# ============================================================================
# Paths
# ============================================================================

ROOT = Path("/home/mazekai/IBD_EcoTyper")

UNIFIED_GENE_FILE = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "07_unified_reference"
    / "unified_reference_gene_metadata.parquet"
)

GSE_GENE_FILE = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "06_raw_count_extraction"
    / "GSE282122"
    / "GSE282122_gene_metadata.parquet"
)

MSCCR_RAW = (
    ROOT
    / "01_raw_data"
    / "bulk"
    / "MSCCR_GSE193677"
    / "GSE193677_MSCCR_Biopsy_counts.txt.gz"
)

OUT_DIR = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "09_bayesprism_input_harmonization"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

MAP_PARQUET = (
    OUT_DIR
    / "unified_reference_to_MSCCR_gene_map.parquet"
)

MAP_TSV = (
    OUT_DIR
    / "unified_reference_to_MSCCR_gene_map.tsv"
)

AUDIT_JSON = (
    OUT_DIR
    / "MSCCR_reference_gene_overlap_audit.json"
)

AUDIT_TXT = (
    OUT_DIR
    / "MSCCR_reference_gene_overlap_audit.txt"
)


# ============================================================================
# Helpers
# ============================================================================

def read_parquet_columns(path):
    table = pq.read_table(path)
    return table.to_pydict()


def strip_ensembl_version(x):
    """
    ENSG00000123456.12 -> ENSG00000123456

    Used only for audit comparison.
    We do NOT silently replace identifiers.
    """
    if x is None:
        return None

    if "." in x:
        return x.split(".", 1)[0]

    return x


def read_msccr_gene_ids(path):

    genes = []

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="strict"
    ) as fh:

        # Header contains 2490 quoted sample IDs
        header = fh.readline()

        if not header:
            raise RuntimeError(
                f"Empty MSCCR file: {path}"
            )

        for line_number, line in enumerate(
            fh,
            start=2
        ):

            line = line.strip()

            if not line:
                continue

            pieces = line.split(
                maxsplit=1
            )

            if len(pieces) != 2:
                raise RuntimeError(
                    "Malformed MSCCR row at "
                    f"line {line_number}"
                )

            gene_id = pieces[0].strip('"')

            genes.append(gene_id)

    return genes


def duplicated_values(values):

    seen = set()
    dup = set()

    for x in values:
        if x in seen:
            dup.add(x)
        else:
            seen.add(x)

    return sorted(dup)


# ============================================================================
# Main
# ============================================================================

def main():

    print()
    print("=" * 72)
    print("MSCCR <-> UNIFIED REFERENCE GENE OVERLAP AUDIT")
    print("=" * 72)

    for path in (
        UNIFIED_GENE_FILE,
        GSE_GENE_FILE,
        MSCCR_RAW
    ):
        if not path.exists():
            raise FileNotFoundError(path)


    # ========================================================================
    # 1. Unified reference
    # ========================================================================

    print("\n[1] Loading unified reference genes...")

    unified = read_parquet_columns(
        UNIFIED_GENE_FILE
    )

    if "gene_symbol" not in unified:
        raise RuntimeError(
            "unified_reference_gene_metadata.parquet "
            "does not contain gene_symbol"
        )

    unified_symbols = unified[
        "gene_symbol"
    ]

    unified_missing = sum(
        x is None or x == ""
        for x in unified_symbols
    )

    unified_duplicates = duplicated_values(
        unified_symbols
    )

    print(
        "Unified genes       :",
        len(unified_symbols)
    )

    print(
        "Missing symbols     :",
        unified_missing
    )

    print(
        "Duplicated symbols  :",
        len(unified_duplicates)
    )

    assert len(unified_symbols) == 14431
    assert unified_missing == 0
    assert len(unified_duplicates) == 0

    print("Unified gene metadata: PASS")


    # ========================================================================
    # 2. GSE282122 symbol <-> Ensembl bridge
    # ========================================================================

    print(
        "\n[2] Loading GSE282122 gene annotation bridge..."
    )

    gse = read_parquet_columns(
        GSE_GENE_FILE
    )

    required = {
        "source_gene_id",
        "gene_id",
        "gene_symbol"
    }

    missing_cols = (
        required
        - set(gse.keys())
    )

    if missing_cols:
        raise RuntimeError(
            "Missing GSE282122 columns: "
            + ", ".join(
                sorted(missing_cols)
            )
        )

    gse_symbols = gse[
        "gene_symbol"
    ]

    gse_gene_ids = gse[
        "gene_id"
    ]

    gse_source_ids = gse[
        "source_gene_id"
    ]

    if not (
        len(gse_symbols)
        == len(gse_gene_ids)
        == len(gse_source_ids)
    ):
        raise RuntimeError(
            "GSE282122 gene metadata columns "
            "have inconsistent lengths"
        )

    gse_symbol_duplicates = duplicated_values(
        gse_symbols
    )

    gse_gene_id_duplicates = duplicated_values(
        gse_gene_ids
    )

    print(
        "GSE genes           :",
        len(gse_symbols)
    )

    print(
        "Duplicate symbols   :",
        len(gse_symbol_duplicates)
    )

    print(
        "Duplicate Ensembl   :",
        len(gse_gene_id_duplicates)
    )

    # Unified construction already required unique
    # symbols in the GSE reference gene space.
    assert len(gse_symbol_duplicates) == 0

    symbol_to_gse_index = {
        symbol: i
        for i, symbol
        in enumerate(gse_symbols)
    }


    # ========================================================================
    # 3. Map 14,431 unified symbols -> GSE Ensembl
    # ========================================================================

    print(
        "\n[3] Mapping unified symbols to Ensembl IDs..."
    )

    mapped_gene_ids = []
    mapped_source_ids = []

    missing_from_bridge = []

    for symbol in unified_symbols:

        idx = symbol_to_gse_index.get(
            symbol
        )

        if idx is None:

            mapped_gene_ids.append(None)
            mapped_source_ids.append(None)

            missing_from_bridge.append(
                symbol
            )

        else:

            mapped_gene_ids.append(
                gse_gene_ids[idx]
            )

            mapped_source_ids.append(
                gse_source_ids[idx]
            )

    print(
        "Mapped to Ensembl   :",
        len(unified_symbols)
        - len(missing_from_bridge)
    )

    print(
        "Missing from bridge :",
        len(missing_from_bridge)
    )

    assert len(missing_from_bridge) == 0

    mapped_ensembl_duplicates = duplicated_values(
        mapped_gene_ids
    )

    print(
        "Duplicate mapped IDs:",
        len(mapped_ensembl_duplicates)
    )

    assert len(mapped_ensembl_duplicates) == 0

    print(
        "Unified symbol -> Ensembl mapping: PASS"
    )


    # ========================================================================
    # 4. MSCCR Ensembl universe
    # ========================================================================

    print(
        "\n[4] Reading MSCCR raw-count gene IDs..."
    )

    msccr_gene_ids = read_msccr_gene_ids(
        MSCCR_RAW
    )

    msccr_duplicates = duplicated_values(
        msccr_gene_ids
    )

    print(
        "MSCCR genes         :",
        len(msccr_gene_ids)
    )

    print(
        "Duplicate MSCCR IDs :",
        len(msccr_duplicates)
    )

    assert len(msccr_gene_ids) == 56632
    assert len(msccr_duplicates) == 0

    msccr_set = set(
        msccr_gene_ids
    )

    msccr_index = {
        gene_id: i + 1
        for i, gene_id
        in enumerate(msccr_gene_ids)
    }

    print("MSCCR gene universe: PASS")


    # ========================================================================
    # 5. Exact overlap
    # ========================================================================

    print(
        "\n[5] Exact Ensembl overlap..."
    )

    exact_present = [
        gene_id in msccr_set
        for gene_id
        in mapped_gene_ids
    ]

    n_exact = sum(
        exact_present
    )

    n_exact_missing = (
        len(mapped_gene_ids)
        - n_exact
    )

    exact_fraction = (
        n_exact
        / len(mapped_gene_ids)
    )

    print(
        "Unified genes       :",
        len(mapped_gene_ids)
    )

    print(
        "Exact MSCCR overlap :",
        n_exact
    )

    print(
        "Missing from MSCCR  :",
        n_exact_missing
    )

    print(
        "Overlap fraction    :",
        f"{exact_fraction:.6f}"
    )


    # ========================================================================
    # 6. Version-stripped diagnostic only
    # ========================================================================

    print(
        "\n[6] Ensembl-version diagnostic..."
    )

    unified_stripped = [
        strip_ensembl_version(x)
        for x in mapped_gene_ids
    ]

    msccr_stripped = [
        strip_ensembl_version(x)
        for x in msccr_gene_ids
    ]

    msccr_stripped_duplicates = (
        duplicated_values(
            msccr_stripped
        )
    )

    msccr_stripped_set = set(
        msccr_stripped
    )

    stripped_present = [
        gene_id in msccr_stripped_set
        for gene_id
        in unified_stripped
    ]

    n_stripped = sum(
        stripped_present
    )

    print(
        "Exact overlap       :",
        n_exact
    )

    print(
        "Versionless overlap :",
        n_stripped
    )

    print(
        "MSCCR duplicate IDs "
        "after stripping     :",
        len(msccr_stripped_duplicates)
    )

    if n_stripped > n_exact:
        print(
            "NOTE: Ensembl version suffixes "
            "affect overlap."
        )
    else:
        print(
            "Ensembl version suffixes do not "
            "increase overlap."
        )


    # ========================================================================
    # 7. Construct formal mapping table
    # ========================================================================

    print(
        "\n[7] Writing formal gene map..."
    )

    msccr_row_index = [
        msccr_index.get(
            gene_id
        )
        for gene_id
        in mapped_gene_ids
    ]

    mapping = pa.table({

        "unified_gene_order":
            list(
                range(
                    1,
                    len(unified_symbols) + 1
                )
            ),

        "gene_symbol":
            unified_symbols,

        "ensembl_id":
            mapped_gene_ids,

        "gse282122_source_gene_id":
            mapped_source_ids,

        "present_in_msccr_exact":
            exact_present,

        "msccr_row_index":
            msccr_row_index,
    })

    pq.write_table(
        mapping,
        MAP_PARQUET
    )

    # TSV is small enough to keep as
    # a human-readable audit artifact.
    with open(
        MAP_TSV,
        "w",
        encoding="utf-8"
    ) as fh:

        fh.write(
            "unified_gene_order\t"
            "gene_symbol\t"
            "ensembl_id\t"
            "gse282122_source_gene_id\t"
            "present_in_msccr_exact\t"
            "msccr_row_index\n"
        )

        for i in range(
            len(unified_symbols)
        ):

            row_index = (
                ""
                if msccr_row_index[i] is None
                else str(
                    msccr_row_index[i]
                )
            )

            fh.write(
                f"{i + 1}\t"
                f"{unified_symbols[i]}\t"
                f"{mapped_gene_ids[i]}\t"
                f"{mapped_source_ids[i]}\t"
                f"{exact_present[i]}\t"
                f"{row_index}\n"
            )


    # ========================================================================
    # 8. Audit summary
    # ========================================================================

    audit = {

        "unified_reference_genes":
            len(unified_symbols),

        "gse282122_annotation_genes":
            len(gse_symbols),

        "msccr_raw_genes":
            len(msccr_gene_ids),

        "unified_missing_gene_symbol":
            unified_missing,

        "unified_duplicate_gene_symbol":
            len(unified_duplicates),

        "unified_missing_from_gse_bridge":
            len(missing_from_bridge),

        "mapped_ensembl_duplicates":
            len(mapped_ensembl_duplicates),

        "msccr_duplicate_gene_ids":
            len(msccr_duplicates),

        "exact_overlap_genes":
            n_exact,

        "exact_overlap_fraction":
            exact_fraction,

        "missing_from_msccr":
            n_exact_missing,

        "versionless_overlap_genes":
            n_stripped,

        "msccr_duplicates_after_version_strip":
            len(
                msccr_stripped_duplicates
            ),
    }

    with open(
        AUDIT_JSON,
        "w",
        encoding="utf-8"
    ) as fh:

        json.dump(
            audit,
            fh,
            indent=2
        )

    with open(
        AUDIT_TXT,
        "w",
        encoding="utf-8"
    ) as fh:

        fh.write(
            "MSCCR <-> UNIFIED REFERENCE "
            "GENE OVERLAP AUDIT\n"
        )

        fh.write(
            "=" * 60 + "\n\n"
        )

        for key, value in audit.items():

            fh.write(
                f"{key}: {value}\n"
            )


    # ========================================================================
    # Final technical gates
    # ========================================================================

    assert unified_missing == 0
    assert len(unified_duplicates) == 0
    assert len(missing_from_bridge) == 0
    assert len(mapped_ensembl_duplicates) == 0
    assert len(msccr_duplicates) == 0
    assert n_exact > 0

    print()
    print("=" * 72)
    print("GENE HARMONIZATION AUDIT COMPLETE")
    print("=" * 72)

    print(
        "Unified genes       :",
        len(unified_symbols)
    )

    print(
        "Exact overlap       :",
        n_exact
    )

    print(
        "Missing from MSCCR  :",
        n_exact_missing
    )

    print(
        "Overlap fraction    :",
        f"{exact_fraction:.6f}"
    )

    print(
        "Gene map            :",
        MAP_PARQUET
    )

    print()
    print(
        "STRUCTURAL GENE-MAPPING CHECKS: PASS"
    )


if __name__ == "__main__":
    main()
