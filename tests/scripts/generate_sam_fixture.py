#!/usr/bin/env python3
import argparse
import csv
import random
from dataclasses import dataclass
from typing import List, Optional, Tuple

try:
    import pysam
except ImportError as exc:
    raise SystemExit(
        "pysam is required to run tests. Install via: pip install pysam"
    ) from exc

BASES = "ACGT"


@dataclass
class Molecule:
    cb: str
    ub: str
    tid: int
    pos: int
    strand: int
    mapq: int


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def generate_barcodes(rng: random.Random, count: int) -> List[str]:
    return [f"{random_dna(rng, 16)}-1" for _ in range(count)]


def generate_umis(rng: random.Random, count: int) -> List[str]:
    return [random_dna(rng, 12) for _ in range(count)]


def generate_genes(rng: random.Random, count: int) -> List[Tuple[str, str]]:
    genes = []
    for idx in range(count):
        gene_id = f"ENSG{rng.randint(100000, 999999)}{idx:03d}"
        gene_name = f"GENE{idx:03d}"
        genes.append((gene_id, gene_name))
    return genes


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate SAM/BAM fixtures for scbamop tests."
    )
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output-bam", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--num-reads", type=int, default=200)
    parser.add_argument("--num-refs", type=int, default=3)
    parser.add_argument("--ref-length", type=int, default=100000)
    parser.add_argument("--read-length", type=int, default=90)
    parser.add_argument("--barcode-count", type=int, default=30)
    parser.add_argument("--metadata-count", type=int, default=20)
    parser.add_argument("--label-count", type=int, default=5)
    parser.add_argument("--dup-rate", type=float, default=0.2)
    parser.add_argument("--missing-cb-rate", type=float, default=0.05)
    parser.add_argument("--missing-ub-rate", type=float, default=0.05)
    parser.add_argument("--cb-unknown-rate", type=float, default=0.1)
    parser.add_argument("--secondary-rate", type=float, default=0.05)
    parser.add_argument("--supplementary-rate", type=float, default=0.05)
    parser.add_argument("--mapq-min", type=int, default=0)
    parser.add_argument("--mapq-max", type=int, default=60)
    parser.add_argument("--mapq-tie-rate", type=float, default=0.2)
    parser.add_argument("--unsafe-label", action="store_true")
    return parser.parse_args()


def build_header(num_refs: int, ref_length: int) -> dict:
    return {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": f"chr{idx + 1}", "LN": ref_length} for idx in range(num_refs)],
    }


def write_metadata(metadata_path: str, barcodes: List[str], labels: List[str]) -> None:
    with open(metadata_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile, lineterminator="\n")
        writer.writerow(["barcode", "label"])
        for barcode, label in zip(barcodes, labels):
            writer.writerow([barcode, label])


def choose_duplicate_mapq(
    rng: random.Random,
    base_mapq: int,
    mapq_min: int,
    mapq_max: int,
    mapq_tie_rate: float,
) -> int:
    if rng.random() < mapq_tie_rate:
        return base_mapq
    delta = rng.choice([-10, -5, 0, 5, 10])
    return max(mapq_min, min(mapq_max, base_mapq + delta))


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    if args.metadata_count > args.barcode_count:
        raise SystemExit("metadata-count cannot exceed barcode-count")

    header = build_header(args.num_refs, args.ref_length)
    barcodes = generate_barcodes(rng, args.barcode_count)
    metadata_barcodes = barcodes[: args.metadata_count]
    unknown_barcodes = generate_barcodes(rng, max(5, args.barcode_count // 3))
    umis = generate_umis(rng, max(30, args.barcode_count))
    genes = generate_genes(rng, max(10, args.label_count))

    labels = [f"cluster_{idx + 1}" for idx in range(args.label_count)]
    if args.unsafe_label and labels:
        labels[0] = "../bad/label"

    metadata_labels = [rng.choice(labels) for _ in metadata_barcodes]
    write_metadata(args.metadata, metadata_barcodes, metadata_labels)

    molecules: List[Molecule] = []

    with pysam.AlignmentFile(args.output_bam, "wb", header=header) as out_bam:
        for idx in range(args.num_reads):
            is_duplicate = rng.random() < args.dup_rate and molecules
            cb: Optional[str] = None
            ub: Optional[str] = None

            if is_duplicate:
                chosen = rng.choice(molecules)
                cb = chosen.cb
                ub = chosen.ub
                tid = chosen.tid
                pos = chosen.pos
                strand = chosen.strand
                mapq = choose_duplicate_mapq(
                    rng, chosen.mapq, args.mapq_min, args.mapq_max, args.mapq_tie_rate
                )
            else:
                if rng.random() >= args.missing_cb_rate:
                    if rng.random() < args.cb_unknown_rate:
                        cb = rng.choice(unknown_barcodes)
                    else:
                        cb = rng.choice(metadata_barcodes)

                if rng.random() >= args.missing_ub_rate:
                    ub = rng.choice(umis)

                tid = rng.randrange(args.num_refs)
                max_start = max(1, args.ref_length - args.read_length + 1)
                pos = rng.randrange(max_start)
                strand = 1 if rng.random() < 0.5 else 0
                mapq = rng.randint(args.mapq_min, args.mapq_max)

            read = pysam.AlignedSegment()
            read.query_name = f"read_{idx:06d}_{args.seed}"
            flag = 0
            if strand == 1:
                flag |= 0x10
            if rng.random() < args.secondary_rate:
                flag |= 0x100
            if rng.random() < args.supplementary_rate:
                flag |= 0x800
            read.flag = flag

            if cb and ub and cb in metadata_barcodes and not (flag & (0x100 | 0x800)):
                molecules.append(
                    Molecule(cb=cb, ub=ub, tid=tid, pos=pos, strand=strand, mapq=mapq)
                )
            read.reference_id = tid
            read.reference_start = pos
            read.mapping_quality = mapq
            # Keep a simple CIGAR for now; extend here if CIGAR handling expands later.
            read.cigar = [(0, args.read_length)]
            read.query_sequence = random_dna(rng, args.read_length)
            read.query_qualities = pysam.qualitystring_to_array("I" * args.read_length)

            gene_id, gene_name = rng.choice(genes)
            nm_tag = rng.randint(0, 3)
            as_tag = rng.randint(0, 100)
            tags: List[Tuple[str, object]] = [
                ("NM", nm_tag),
                ("AS", as_tag),
                ("GX", gene_id),
                ("GN", gene_name),
            ]
            if cb is not None:
                tags.append(("CB", cb))
            if ub is not None:
                tags.append(("UB", ub))

            read.set_tags(tags)
            out_bam.write(read)


if __name__ == "__main__":
    main()
