#!/usr/bin/env python3
import argparse
import csv
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

try:
    import pysam
except ImportError as exc:
    raise SystemExit(
        "pysam is required to run tests. Install via: pip install pysam"
    ) from exc


@dataclass
class ReadSignature:
    qname: str
    flag: int
    tid: int
    pos: int
    mapq: int
    cigar: Optional[str]
    next_tid: int
    next_pos: int
    tlen: int
    seq: Optional[str]
    qual: Optional[str]
    tags: List[Tuple[str, object, str]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate scbamop split output.")
    parser.add_argument("--input-bam", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--mode", choices=["split", "dedup"], required=True)
    parser.add_argument("--mapq-threshold", type=int, default=0)
    return parser.parse_args()


def sanitize_label(label: str) -> str:
    sanitized = list(label)
    for idx, char in enumerate(sanitized):
        if char in {"/", "\\", "~"}:
            sanitized[idx] = "_"
        elif not (char.isalnum() or char in {"_", "-", " ", "."}):
            sanitized[idx] = "_"

    if sanitized and sanitized[0] == ".":
        sanitized[0] = "_"

    sanitized_str = "".join(sanitized)
    while ".." in sanitized_str:
        sanitized_str = sanitized_str.replace("..", "__")
    return sanitized_str


def load_metadata(metadata_path: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open(metadata_path, newline="") as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            barcode = row[0].strip()
            label = row[1].strip()
            mapping[barcode] = sanitize_label(label)
    return mapping


def read_signature(read: pysam.AlignedSegment) -> ReadSignature:
    return ReadSignature(
        qname=read.query_name,
        flag=read.flag,
        tid=read.reference_id,
        pos=read.reference_start,
        mapq=read.mapping_quality,
        cigar=read.cigarstring,
        next_tid=read.next_reference_id,
        next_pos=read.next_reference_start,
        tlen=read.template_length,
        seq=read.query_sequence,
        qual=read.qual,
        tags=read.get_tags(with_value_type=True),
    )


def extract_tag(read: pysam.AlignedSegment, tag: str) -> Optional[str]:
    try:
        return read.get_tag(tag)
    except KeyError:
        return None


def passes_filters(
    read: pysam.AlignedSegment, cb_map: Dict[str, str], mapq_threshold: int, mode: str
) -> bool:
    cb = extract_tag(read, "CB")
    ub = extract_tag(read, "UB")

    if cb is None or ub is None:
        return False
    if cb not in cb_map:
        return False
    if read.mapping_quality < mapq_threshold:
        return False
    if mode == "dedup" and (read.is_secondary or read.is_supplementary):
        return False
    return True


def gather_expected_reads(
    input_bam: str, cb_map: Dict[str, str], mapq_threshold: int, mode: str
) -> Dict[str, List[ReadSignature]]:
    expected: Dict[str, List[ReadSignature]] = {
        label: [] for label in set(cb_map.values())
    }

    if mode == "split":
        with pysam.AlignmentFile(input_bam, "rb") as infile:
            for read in infile:
                if not passes_filters(read, cb_map, mapq_threshold, mode):
                    continue
                cb = extract_tag(read, "CB")
                if cb is None:
                    continue
                label = cb_map[cb]
                expected[label].append(read_signature(read))
        return expected

    best_by_key: Dict[Tuple[str, str, int, int, int], Tuple[int, int]] = {}
    read_index = 0
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        for read in infile:
            if not passes_filters(read, cb_map, mapq_threshold, mode):
                read_index += 1
                continue
            cb = extract_tag(read, "CB")
            ub = extract_tag(read, "UB")
            if cb is None or ub is None:
                read_index += 1
                continue
            strand = 1 if read.is_reverse else 0
            key = (cb, ub, read.reference_id, read.reference_start, strand)
            entry = best_by_key.get(key)
            if entry is None:
                best_by_key[key] = (read_index, read.mapping_quality)
            else:
                best_idx, best_mapq = entry
                if read.mapping_quality > best_mapq:
                    best_by_key[key] = (read_index, read.mapping_quality)
            read_index += 1

    read_index = 0
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        for read in infile:
            if not passes_filters(read, cb_map, mapq_threshold, mode):
                read_index += 1
                continue
            cb = extract_tag(read, "CB")
            ub = extract_tag(read, "UB")
            if cb is None or ub is None:
                read_index += 1
                continue
            strand = 1 if read.is_reverse else 0
            key = (cb, ub, read.reference_id, read.reference_start, strand)
            best_entry = best_by_key.get(key)
            if best_entry and best_entry[0] == read_index:
                label = cb_map[cb]
                expected[label].append(read_signature(read))
            read_index += 1

    return expected


def compare_headers(input_bam: str, output_bam: str, mode: str) -> None:
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        expected_header = infile.header.to_dict()

    with pysam.AlignmentFile(output_bam, "rb") as outfile:
        actual_header = outfile.header.to_dict()

    strip_sort_order(expected_header)
    strip_sort_order(actual_header)

    if actual_header != expected_header:
        raise AssertionError(f"Header mismatch for {output_bam}")


def strip_sort_order(header: Dict[str, object]) -> None:
    hd = header.get("HD")
    if isinstance(hd, dict) and "SO" in hd:
        hd = dict(hd)
        hd.pop("SO", None)
        if hd:
            header["HD"] = hd
        else:
            header.pop("HD", None)


def compare_reads(expected: List[ReadSignature], output_bam: str) -> None:
    with pysam.AlignmentFile(output_bam, "rb") as outfile:
        actual_reads = [read_signature(read) for read in outfile]

    if len(actual_reads) != len(expected):
        raise AssertionError(
            f"Read count mismatch for {output_bam}: expected {len(expected)}, got {len(actual_reads)}"
        )

    for idx, (exp, act) in enumerate(zip(expected, actual_reads)):
        if exp != act:
            raise AssertionError(
                f"Read mismatch at index {idx} in {output_bam}: expected {exp}, got {act}"
            )


def main() -> None:
    args = parse_args()
    cb_map = load_metadata(args.metadata)

    expected = gather_expected_reads(
        args.input_bam,
        cb_map,
        args.mapq_threshold,
        args.mode,
    )

    expected_labels = set(expected.keys())
    output_files = {
        name for name in os.listdir(args.output_dir) if name.endswith(".bam")
    }
    expected_files = {f"{label}.bam" for label in expected_labels}

    if output_files != expected_files:
        raise AssertionError(
            f"Output file mismatch: expected {sorted(expected_files)}, got {sorted(output_files)}"
        )

    for label, expected_reads in expected.items():
        output_bam = os.path.join(args.output_dir, f"{label}.bam")
        if not os.path.exists(output_bam):
            raise AssertionError(f"Missing output BAM: {output_bam}")
        compare_headers(args.input_bam, output_bam, args.mode)
        compare_reads(expected_reads, output_bam)


if __name__ == "__main__":
    main()
