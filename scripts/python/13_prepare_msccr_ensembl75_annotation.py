#!/usr/bin/env python3

"""
13_prepare_msccr_ensembl75_annotation.py

Purpose
-------
Prepare the MSCCR / GSE193677 raw bulk count matrix for downstream
gene-identifier harmonization.

The script:

1. Downloads the exact Ensembl release 75 / GRCh37 GTF used for the
   original MSCCR gene quantification.
2. Parses gene-level Ensembl annotations.
3. Maps MSCCR Ensembl gene IDs to Ensembl75 gene_name.
4. Audits unmapped IDs and duplicated gene symbols.
5. Creates a conservative gene-symbol count matrix containing only
   unambiguous 1:1 mappings.
6. Audits overlap with the unified scRNA reference.

IMPORTANT
---------
This script does NOT:
- normalize counts
- log-transform counts
- round fractional values
- perform CPM/TMM/TPM/FPKM/voom
- perform expression filtering
- perform differential expression
- modify the original MSCCR matrix

For retained genes, numeric count payloads are copied verbatim.
"""

import csv
import gzip
import hashlib
import json
import re
import shlex
import shutil
import urllib.request

from collections import Counter, defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


# ============================================================================
# Paths
# ============================================================================

ROOT = Path(
    "/home/mazekai/IBD_EcoTyper"
)

MSCCR_RAW = (
    ROOT
    / "01_raw_data"
    / "bulk"
    / "MSCCR_GSE193677"
    / "GSE193677_MSCCR_Biopsy_counts.txt.gz"
)

UNIFIED_GENE_FILE = (
    ROOT
    / "03_reference"
    / "combined_reference"
    / "output"
    / "07_unified_reference"
    / "unified_reference_gene_metadata.parquet"
)


ANNOT_DIR = (
    ROOT
    / "01_raw_data"
    / "annotation"
    / "Ensembl_GRCh37_release75"
)

ANNOT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


GTF_URL = (
    "https://ftp.ensembl.org/"
    "pub/release-75/gtf/homo_sapiens/"
    "Homo_sapiens.GRCh37.75.gtf.gz"
)

GTF_FILE = (
    ANNOT_DIR
    / "Homo_sapiens.GRCh37.75.gtf.gz"
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


ANNOT_PARQUET = (
    OUT_DIR
    / "Ensembl75_GRCh37_gene_annotation.parquet"
)

MSCCR_MAP_PARQUET = (
    OUT_DIR
    / "MSCCR_Ensembl75_gene_mapping.parquet"
)

MSCCR_MAP_TSV = (
    OUT_DIR
    / "MSCCR_Ensembl75_gene_mapping.tsv"
)

UNIFIED_MAP_PARQUET = (
    OUT_DIR
    / "unified_reference_to_MSCCR_Ensembl75_gene_map.parquet"
)

UNIFIED_MAP_TSV = (
    OUT_DIR
    / "unified_reference_to_MSCCR_Ensembl75_gene_map.tsv"
)

SYMBOL_COUNTS = (
    OUT_DIR
    / "GSE193677_MSCCR_Biopsy_counts_"
      "GRCh37_Ensembl75_gene_symbol_unique.txt.gz"
)

AUDIT_JSON = (
    OUT_DIR
    / "MSCCR_Ensembl75_annotation_audit.json"
)

AUDIT_TXT = (
    OUT_DIR
    / "MSCCR_Ensembl75_annotation_audit.txt"
)


# ============================================================================
# Utilities
# ============================================================================

def sha256_file(path):

    h = hashlib.sha256()

    with open(path, "rb") as fh:

        while True:

            block = fh.read(
                1024 * 1024
            )

            if not block:
                break

            h.update(block)

    return h.hexdigest()


def download_gtf():

    print()
    print("=" * 72)
    print("STEP 1: DOWNLOAD OFFICIAL ENSEMBL75 GTF")
    print("=" * 72)

    print("Source:")
    print(GTF_URL)

    print()
    print("Destination:")
    print(GTF_FILE)

    if GTF_FILE.exists():

        print()
        print(
            "GTF already exists; "
            "download skipped."
        )

    else:

        tmp = Path(
            str(GTF_FILE) + ".part"
        )

        if tmp.exists():
            tmp.unlink()

        request = urllib.request.Request(
            GTF_URL,
            headers={
                "User-Agent":
                    "IBD-EcoTyper-reference-audit/1.0"
            }
        )

        try:

            with urllib.request.urlopen(
                request,
                timeout=120
            ) as response:

                with open(
                    tmp,
                    "wb"
                ) as out:

                    shutil.copyfileobj(
                        response,
                        out,
                        length=1024 * 1024
                    )

            tmp.rename(
                GTF_FILE
            )

        except Exception:

            if tmp.exists():
                tmp.unlink()

            raise

        print(
            "Download complete."
        )

    if (
        not GTF_FILE.exists()
        or GTF_FILE.stat().st_size == 0
    ):
        raise RuntimeError(
            "Downloaded GTF is empty."
        )

    # gzip integrity / readability check
    with gzip.open(
        GTF_FILE,
        "rt",
        encoding="utf-8"
    ) as fh:

        first = fh.readline()

        if not first:
            raise RuntimeError(
                "Cannot read downloaded GTF."
            )

    checksum = sha256_file(
        GTF_FILE
    )

    print()
    print(
        "GTF size (bytes):",
        GTF_FILE.stat().st_size
    )

    print(
        "GTF SHA256      :",
        checksum
    )

    print(
        "GTF download/readability: PASS"
    )

    return checksum


# ============================================================================
# GTF parsing
# ============================================================================

ATTR_PATTERN = re.compile(
    r'(\S+)\s+"([^"]*)";'
)


def parse_attributes(text):

    return dict(
        ATTR_PATTERN.findall(text)
    )


def parse_gtf():

    print()
    print("=" * 72)
    print("STEP 2: PARSE ENSEMBL75 GENE ANNOTATION")
    print("=" * 72)

    records = {}

    duplicate_gene_rows = 0
    malformed_rows = 0
    missing_gene_id = 0
    missing_gene_name = 0

    with gzip.open(
        GTF_FILE,
        "rt",
        encoding="utf-8",
        errors="strict"
    ) as fh:

        for line_no, line in enumerate(
            fh,
            start=1
        ):

            if (
                not line
                or line.startswith("#")
            ):
                continue

            fields = line.rstrip(
                "\r\n"
            ).split("\t")

            if len(fields) != 9:

                malformed_rows += 1
                continue

            chromosome = fields[0]
            feature = fields[2]

            if feature != "gene":
                continue

            start = fields[3]
            end = fields[4]
            strand = fields[6]

            attrs = parse_attributes(
                fields[8]
            )

            gene_id = attrs.get(
                "gene_id",
                ""
            )

            gene_name = attrs.get(
                "gene_name",
                ""
            )

            gene_biotype = (
                attrs.get(
                    "gene_biotype"
                )
                or
                attrs.get(
                    "gene_type"
                )
                or
                ""
            )

            if gene_id == "":

                missing_gene_id += 1
                continue

            if gene_name == "":
                missing_gene_name += 1

            record = {
                "gene_id":
                    gene_id,

                "gene_name":
                    gene_name,

                "gene_biotype":
                    gene_biotype,

                "chromosome":
                    chromosome,

                "start":
                    int(start),

                "end":
                    int(end),

                "strand":
                    strand,
            }

            if gene_id in records:

                duplicate_gene_rows += 1

                # Same gene_id should not have
                # incompatible gene-level records.
                if (
                    records[gene_id]
                    != record
                ):

                    raise RuntimeError(
                        "Conflicting duplicate "
                        f"GTF gene_id: {gene_id}"
                    )

            else:

                records[
                    gene_id
                ] = record

    print(
        "Unique gene IDs       :",
        len(records)
    )

    print(
        "Duplicate gene rows   :",
        duplicate_gene_rows
    )

    print(
        "Malformed GTF rows    :",
        malformed_rows
    )

    print(
        "Missing gene_id       :",
        missing_gene_id
    )

    print(
        "Missing gene_name     :",
        missing_gene_name
    )

    if len(records) == 0:
        raise RuntimeError(
            "No gene records parsed from GTF."
        )

    if malformed_rows != 0:
        raise RuntimeError(
            "Malformed GTF rows detected."
        )

    if missing_gene_id != 0:
        raise RuntimeError(
            "Gene feature without gene_id."
        )

    print(
        "Ensembl75 gene annotation: PASS"
    )

    # Save exact parsed annotation.
    ordered = list(
        records.values()
    )

    table = pa.table({

        "gene_id":
            [
                x["gene_id"]
                for x in ordered
            ],

        "gene_name":
            [
                x["gene_name"]
                for x in ordered
            ],

        "gene_biotype":
            [
                x["gene_biotype"]
                for x in ordered
            ],

        "chromosome":
            [
                x["chromosome"]
                for x in ordered
            ],

        "start":
            [
                x["start"]
                for x in ordered
            ],

        "end":
            [
                x["end"]
                for x in ordered
            ],

        "strand":
            [
                x["strand"]
                for x in ordered
            ],
    })

    pq.write_table(
        table,
        ANNOT_PARQUET
    )

    return (
        records,
        duplicate_gene_rows,
        missing_gene_name
    )


# ============================================================================
# Read MSCCR schema / gene IDs
# ============================================================================

def read_msccr_gene_ids():

    print()
    print("=" * 72)
    print("STEP 3: READ MSCCR GENE IDENTIFIERS")
    print("=" * 72)

    if not MSCCR_RAW.exists():
        raise FileNotFoundError(
            MSCCR_RAW
        )

    gene_ids = []

    with gzip.open(
        MSCCR_RAW,
        "rt",
        encoding="utf-8",
        errors="strict"
    ) as fh:

        header_line = fh.readline()

        if not header_line:
            raise RuntimeError(
                "MSCCR matrix has no header."
            )

        sample_ids = shlex.split(
            header_line.strip()
        )

        for line_no, line in enumerate(
            fh,
            start=2
        ):

            line = line.rstrip(
                "\r\n"
            )

            if not line:
                continue

            parts = line.split(
                maxsplit=1
            )

            if len(parts) != 2:

                raise RuntimeError(
                    "Malformed MSCCR row "
                    f"at line {line_no}"
                )

            gene_id = (
                parts[0]
                .strip('"')
            )

            gene_ids.append(
                gene_id
            )

    duplicate_gene_ids = (
        len(gene_ids)
        -
        len(set(gene_ids))
    )

    duplicate_samples = (
        len(sample_ids)
        -
        len(set(sample_ids))
    )

    print(
        "MSCCR genes          :",
        len(gene_ids)
    )

    print(
        "MSCCR samples        :",
        len(sample_ids)
    )

    print(
        "Duplicate gene IDs   :",
        duplicate_gene_ids
    )

    print(
        "Duplicate sample IDs :",
        duplicate_samples
    )

    if len(gene_ids) != 56632:

        raise RuntimeError(
            "Unexpected MSCCR gene count: "
            f"{len(gene_ids)}"
        )

    if len(sample_ids) != 2490:

        raise RuntimeError(
            "Unexpected MSCCR sample count: "
            f"{len(sample_ids)}"
        )

    if duplicate_gene_ids != 0:
        raise RuntimeError(
            "Duplicated MSCCR Ensembl IDs."
        )

    if duplicate_samples != 0:
        raise RuntimeError(
            "Duplicated MSCCR sample IDs."
        )

    print(
        "MSCCR identifier schema: PASS"
    )

    return (
        gene_ids,
        sample_ids,
        header_line
    )


# ============================================================================
# MSCCR Ensembl -> Ensembl75 gene_name mapping
# ============================================================================

def build_msccr_mapping(
    msccr_gene_ids,
    gtf_records
):

    print()
    print("=" * 72)
    print("STEP 4: MAP MSCCR ENSEMBL IDs TO ENSEMBL75")
    print("=" * 72)

    mapped_symbols = []

    for gene_id in msccr_gene_ids:

        rec = gtf_records.get(
            gene_id
        )

        if rec is None:

            mapped_symbols.append(
                ""
            )

        else:

            mapped_symbols.append(
                rec["gene_name"]
            )

    symbol_counts = Counter(
        x
        for x in mapped_symbols
        if x != ""
    )

    statuses = []

    gene_names = []
    biotypes = []
    chromosomes = []

    n_unique = 0
    n_duplicate_symbol = 0
    n_unmapped = 0
    n_blank_symbol = 0

    for gene_id, gene_name in zip(
        msccr_gene_ids,
        mapped_symbols
    ):

        rec = gtf_records.get(
            gene_id
        )

        if rec is None:

            status = (
                "unmapped_ensembl_id"
            )

            gene_name = ""
            biotype = ""
            chromosome = ""

            n_unmapped += 1

        elif gene_name == "":

            status = (
                "blank_gene_name"
            )

            biotype = rec[
                "gene_biotype"
            ]

            chromosome = rec[
                "chromosome"
            ]

            n_blank_symbol += 1

        elif symbol_counts[
            gene_name
        ] > 1:

            status = (
                "duplicate_gene_symbol"
            )

            biotype = rec[
                "gene_biotype"
            ]

            chromosome = rec[
                "chromosome"
            ]

            n_duplicate_symbol += 1

        else:

            status = (
                "unique_gene_symbol"
            )

            biotype = rec[
                "gene_biotype"
            ]

            chromosome = rec[
                "chromosome"
            ]

            n_unique += 1

        statuses.append(
            status
        )

        gene_names.append(
            gene_name
        )

        biotypes.append(
            biotype
        )

        chromosomes.append(
            chromosome
        )

    print(
        "Unique symbol rows    :",
        n_unique
    )

    print(
        "Duplicate-symbol rows :",
        n_duplicate_symbol
    )

    print(
        "Blank gene_name rows  :",
        n_blank_symbol
    )

    print(
        "Unmapped Ensembl IDs  :",
        n_unmapped
    )

    mapping_table = pa.table({

        "msccr_ensembl_id":
            msccr_gene_ids,

        "ensembl75_gene_name":
            gene_names,

        "ensembl75_gene_biotype":
            biotypes,

        "chromosome":
            chromosomes,

        "mapping_status":
            statuses,
    })

    pq.write_table(
        mapping_table,
        MSCCR_MAP_PARQUET
    )

    with open(
        MSCCR_MAP_TSV,
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:

        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "msccr_ensembl_id",
            "ensembl75_gene_name",
            "ensembl75_gene_biotype",
            "chromosome",
            "mapping_status",
        ])

        for row in zip(
            msccr_gene_ids,
            gene_names,
            biotypes,
            chromosomes,
            statuses
        ):

            writer.writerow(
                row
            )

    print(
        "MSCCR annotation map written."
    )

    return {
        "gene_names":
            gene_names,

        "statuses":
            statuses,

        "symbol_counts":
            symbol_counts,

        "n_unique":
            n_unique,

        "n_duplicate_symbol":
            n_duplicate_symbol,

        "n_blank_symbol":
            n_blank_symbol,

        "n_unmapped":
            n_unmapped,
    }


# ============================================================================
# Unified reference overlap
# ============================================================================

def audit_unified_overlap(
    msccr_gene_ids,
    mapping
):

    print()
    print("=" * 72)
    print("STEP 5: UNIFIED REFERENCE OVERLAP")
    print("=" * 72)

    unified = (
        pq.read_table(
            UNIFIED_GENE_FILE
        )
        .to_pydict()
    )

    if (
        "gene_symbol"
        not in unified
    ):

        raise RuntimeError(
            "Unified gene metadata "
            "does not contain gene_symbol."
        )

    unified_symbols = unified[
        "gene_symbol"
    ]

    if len(unified_symbols) != 14431:

        raise RuntimeError(
            "Unexpected unified gene count."
        )

    unique_symbol_to_ensembl = {}

    symbol_to_msccr_ids = defaultdict(
        list
    )

    for gene_id, symbol, status in zip(
        msccr_gene_ids,
        mapping["gene_names"],
        mapping["statuses"]
    ):

        if symbol != "":

            symbol_to_msccr_ids[
                symbol
            ].append(
                gene_id
            )

        if (
            status
            == "unique_gene_symbol"
        ):

            unique_symbol_to_ensembl[
                symbol
            ] = gene_id

    # Ensembl75 complete symbol universe.
    # Used only to distinguish "not in annotation"
    # from "annotated but not represented in MSCCR".
    annotation_table = (
        pq.read_table(
            ANNOT_PARQUET
        )
        .to_pydict()
    )

    gtf_symbol_set = set(
        x
        for x in annotation_table[
            "gene_name"
        ]
        if x != ""
    )

    statuses = []
    matched_gene_ids = []

    n_unique_overlap = 0
    n_ambiguous = 0
    n_absent_msccr = 0
    n_not_v75 = 0

    for symbol in unified_symbols:

        msccr_ids = (
            symbol_to_msccr_ids.get(
                symbol,
                []
            )
        )

        if len(msccr_ids) == 1:

            status = (
                "unique_msccr_ensembl75_match"
            )

            matched_gene_id = (
                msccr_ids[0]
            )

            n_unique_overlap += 1

        elif len(msccr_ids) > 1:

            status = (
                "ambiguous_multiple_"
                "msccr_ensembl_ids"
            )

            matched_gene_id = ""

            n_ambiguous += 1

        elif symbol in gtf_symbol_set:

            status = (
                "present_in_ensembl75_"
                "but_absent_from_msccr"
            )

            matched_gene_id = ""

            n_absent_msccr += 1

        else:

            status = (
                "symbol_not_in_ensembl75"
            )

            matched_gene_id = ""

            n_not_v75 += 1

        statuses.append(
            status
        )

        matched_gene_ids.append(
            matched_gene_id
        )

    print(
        "Unified genes              :",
        len(unified_symbols)
    )

    print(
        "Unique Ensembl75 overlap   :",
        n_unique_overlap
    )

    print(
        "Ambiguous symbol mappings  :",
        n_ambiguous
    )

    print(
        "Annotated but absent MSCCR :",
        n_absent_msccr
    )

    print(
        "Symbol absent from v75     :",
        n_not_v75
    )

    print(
        "Usable overlap fraction    :",
        f"{n_unique_overlap / len(unified_symbols):.6f}"
    )

    overlap_table = pa.table({

        "unified_gene_order":
            list(
                range(
                    1,
                    len(unified_symbols) + 1
                )
            ),

        "gene_symbol":
            unified_symbols,

        "msccr_ensembl75_gene_id":
            matched_gene_ids,

        "mapping_status":
            statuses,
    })

    pq.write_table(
        overlap_table,
        UNIFIED_MAP_PARQUET
    )

    with open(
        UNIFIED_MAP_TSV,
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:

        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writerow([
            "unified_gene_order",
            "gene_symbol",
            "msccr_ensembl75_gene_id",
            "mapping_status",
        ])

        for i in range(
            len(unified_symbols)
        ):

            writer.writerow([
                i + 1,
                unified_symbols[i],
                matched_gene_ids[i],
                statuses[i],
            ])

    return {
        "unified_genes":
            len(unified_symbols),

        "unique_overlap":
            n_unique_overlap,

        "ambiguous":
            n_ambiguous,

        "annotated_absent_msccr":
            n_absent_msccr,

        "not_in_ensembl75":
            n_not_v75,
    }


# ============================================================================
# Create conservative gene-symbol raw count matrix
# ============================================================================

def write_gene_symbol_matrix(
    msccr_gene_ids,
    mapping
):

    print()
    print("=" * 72)
    print("STEP 6: WRITE GENE-SYMBOL RAW COUNT MATRIX")
    print("=" * 72)

    keep_map = {}

    for gene_id, symbol, status in zip(
        msccr_gene_ids,
        mapping["gene_names"],
        mapping["statuses"]
    ):

        if (
            status
            == "unique_gene_symbol"
        ):

            keep_map[
                gene_id
            ] = symbol

    print(
        "Rows retained:",
        len(keep_map)
    )

    source_payload_hash = (
        hashlib.sha256()
    )

    n_written = 0

    with gzip.open(
        MSCCR_RAW,
        "rt",
        encoding="utf-8",
        errors="strict"
    ) as src:

        with gzip.open(
            SYMBOL_COUNTS,
            "wt",
            encoding="utf-8",
            newline=""
        ) as dst:

            header = src.readline()

            # Preserve sample header exactly.
            dst.write(header)

            for line_no, line in enumerate(
                src,
                start=2
            ):

                stripped = line.rstrip(
                    "\r\n"
                )

                if not stripped:
                    continue

                parts = stripped.split(
                    maxsplit=1
                )

                if len(parts) != 2:

                    raise RuntimeError(
                        "Malformed source row "
                        f"at line {line_no}"
                    )

                gene_id = (
                    parts[0]
                    .strip('"')
                )

                payload = parts[1]

                symbol = keep_map.get(
                    gene_id
                )

                if symbol is None:
                    continue

                # Numeric payload is copied
                # verbatim. Only the row ID
                # changes from Ensembl ID to
                # gene symbol.
                source_payload_hash.update(
                    payload.encode(
                        "utf-8"
                    )
                )

                source_payload_hash.update(
                    b"\n"
                )

                dst.write(
                    f'"{symbol}" '
                )

                dst.write(
                    payload
                )

                dst.write(
                    "\n"
                )

                n_written += 1

    if n_written != len(
        keep_map
    ):

        raise RuntimeError(
            "Output row-count mismatch."
        )

    # Independent re-read of the output
    # numeric payload.
    output_payload_hash = (
        hashlib.sha256()
    )

    verify_rows = 0

    with gzip.open(
        SYMBOL_COUNTS,
        "rt",
        encoding="utf-8",
        errors="strict"
    ) as fh:

        fh.readline()

        seen_symbols = set()

        for line in fh:

            stripped = line.rstrip(
                "\r\n"
            )

            if not stripped:
                continue

            parts = stripped.split(
                maxsplit=1
            )

            if len(parts) != 2:
                raise RuntimeError(
                    "Malformed converted row."
                )

            symbol = (
                parts[0]
                .strip('"')
            )

            payload = parts[1]

            if symbol in seen_symbols:

                raise RuntimeError(
                    "Duplicate gene symbol "
                    "in converted matrix: "
                    f"{symbol}"
                )

            seen_symbols.add(
                symbol
            )

            output_payload_hash.update(
                payload.encode(
                    "utf-8"
                )
            )

            output_payload_hash.update(
                b"\n"
            )

            verify_rows += 1

    hash_match = (
        source_payload_hash.hexdigest()
        ==
        output_payload_hash.hexdigest()
    )

    print(
        "Written rows          :",
        n_written
    )

    print(
        "Verified rows         :",
        verify_rows
    )

    print(
        "Duplicate output rows :",
        0
    )

    print(
        "Numeric payload hash  :",
        "PASS"
        if hash_match
        else "FAIL"
    )

    if verify_rows != n_written:
        raise RuntimeError(
            "Converted matrix "
            "verification failed."
        )

    if not hash_match:
        raise RuntimeError(
            "Count payload changed "
            "during annotation conversion."
        )

    print()
    print(
        "COUNT VALUES PRESERVED EXACTLY: PASS"
    )

    return {
        "output_rows":
            n_written,

        "payload_hash":
            source_payload_hash.hexdigest(),

        "payload_concordance":
            hash_match,
    }


# ============================================================================
# Main
# ============================================================================

def main():

    if not UNIFIED_GENE_FILE.exists():
        raise FileNotFoundError(
            UNIFIED_GENE_FILE
        )

    gtf_sha256 = (
        download_gtf()
    )

    (
        gtf_records,
        gtf_duplicate_rows,
        gtf_missing_gene_name
    ) = parse_gtf()

    (
        msccr_gene_ids,
        sample_ids,
        _
    ) = read_msccr_gene_ids()

    mapping = build_msccr_mapping(
        msccr_gene_ids,
        gtf_records
    )

    overlap = audit_unified_overlap(
        msccr_gene_ids,
        mapping
    )

    converted = (
        write_gene_symbol_matrix(
            msccr_gene_ids,
            mapping
        )
    )

    audit = {

        "annotation_source":
            GTF_URL,

        "annotation_file":
            str(GTF_FILE),

        "annotation_sha256":
            gtf_sha256,

        "annotation_release":
            "Ensembl 75",

        "genome_build":
            "GRCh37",

        "gtf_unique_gene_ids":
            len(gtf_records),

        "gtf_duplicate_gene_rows":
            gtf_duplicate_rows,

        "gtf_missing_gene_name":
            gtf_missing_gene_name,

        "msccr_input_genes":
            len(msccr_gene_ids),

        "msccr_samples":
            len(sample_ids),

        "msccr_unique_gene_symbol_rows":
            mapping["n_unique"],

        "msccr_duplicate_symbol_rows":
            mapping[
                "n_duplicate_symbol"
            ],

        "msccr_blank_gene_name_rows":
            mapping[
                "n_blank_symbol"
            ],

        "msccr_unmapped_ensembl_ids":
            mapping[
                "n_unmapped"
            ],

        "unified_reference_genes":
            overlap[
                "unified_genes"
            ],

        "unified_unique_overlap":
            overlap[
                "unique_overlap"
            ],

        "unified_ambiguous_symbols":
            overlap[
                "ambiguous"
            ],

        "unified_annotated_but_absent_msccr":
            overlap[
                "annotated_absent_msccr"
            ],

        "unified_symbols_not_in_ensembl75":
            overlap[
                "not_in_ensembl75"
            ],

        "converted_matrix_rows":
            converted[
                "output_rows"
            ],

        "count_payload_sha256":
            converted[
                "payload_hash"
            ],

        "count_payload_concordance":
            converted[
                "payload_concordance"
            ],

        "normalization_performed":
            False,

        "rounding_performed":
            False,

        "expression_filtering_performed":
            False,

        "duplicate_symbol_aggregation_performed":
            False,
    }

    with open(
        AUDIT_JSON,
        "w",
        encoding="utf-8"
    ) as fh:

        json.dump(
            audit,
            fh,
            indent=2,
            ensure_ascii=False
        )

    with open(
        AUDIT_TXT,
        "w",
        encoding="utf-8"
    ) as fh:

        fh.write(
            "MSCCR GRCh37 / ENSEMBL75 "
            "ANNOTATION AUDIT\n"
        )

        fh.write(
            "=" * 60
            + "\n\n"
        )

        for key, value in (
            audit.items()
        ):

            fh.write(
                f"{key}: {value}\n"
            )

    print()
    print("=" * 72)
    print("MSCCR ENSEMBL75 ANNOTATION COMPLETE")
    print("=" * 72)

    print(
        "Original genes        :",
        len(msccr_gene_ids)
    )

    print(
        "Unique-symbol genes   :",
        mapping["n_unique"]
    )

    print(
        "Duplicate-symbol rows :",
        mapping[
            "n_duplicate_symbol"
        ]
    )

    print(
        "Unmapped Ensembl IDs  :",
        mapping[
            "n_unmapped"
        ]
    )

    print(
        "Unified reference     :",
        overlap[
            "unified_genes"
        ]
    )

    print(
        "Unique overlap        :",
        overlap[
            "unique_overlap"
        ]
    )

    print(
        "Ambiguous overlap     :",
        overlap[
            "ambiguous"
        ]
    )

    print(
        "Output raw counts     :",
        SYMBOL_COUNTS
    )

    print()
    print(
        "Annotation source     : "
        "Ensembl release 75 / GRCh37"
    )

    print(
        "Count transformation  : NONE"
    )

    print(
        "Count payload integrity:",
        "PASS"
    )

    print()
    print(
        "ANNOTATION / ID CONVERSION: PASS"
    )


if __name__ == "__main__":
    main()
