#!/usr/bin/env python3

from __future__ import annotations

import argparse
from itertools import zip_longest
from pathlib import Path
import shutil
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from msnaf.core import export_peptides


def compare_text_files(expected_path: Path, actual_path: Path, label: str) -> bool:
    if not actual_path.exists():
        print(f"[FAIL] {label}: missing output {actual_path}")
        return False

    with expected_path.open("r", encoding="utf-8") as expected_handle, actual_path.open(
        "r", encoding="utf-8"
    ) as actual_handle:
        for line_number, (expected_line, actual_line) in enumerate(
            zip_longest(expected_handle, actual_handle, fillvalue=None),
            start=1,
        ):
            if expected_line != actual_line:
                print(f"[FAIL] {label}: first difference at line {line_number}")
                print(f"  expected: {repr(expected_line)}")
                print(f"  actual:   {repr(actual_line)}")
                return False

    print(f"[OK] {label}")
    return True


def run_fixture_comparison(fixture_dir: Path, keep_output: bool = False) -> int:
    counts_path = fixture_dir / "counts.original.full.txt"
    refs_dir = fixture_dir / "snaf-data"
    expected_intermediates = fixture_dir / "snaf_intermediates.tsv"
    expected_stats = fixture_dir / "result" / "NeoJunction_statistics_maxmin.txt"
    genome_fasta = fixture_dir / "genome.fa"

    if not counts_path.exists():
        print(f"missing counts file: {counts_path}")
        return 2
    if not refs_dir.exists():
        print(f"missing refs dir: {refs_dir}")
        return 2
    if not expected_intermediates.exists():
        print(f"missing expected intermediates: {expected_intermediates}")
        return 2
    if not expected_stats.exists():
        print(f"missing expected stats: {expected_stats}")
        return 2
    if not genome_fasta.exists():
        print(f"missing genome fasta: {genome_fasta}")
        return 2

    tmpdir_obj = tempfile.TemporaryDirectory(prefix="msnaf_fixture_compare_")
    tmpdir = Path(tmpdir_obj.name)
    try:
        output_path = tmpdir / "result" / "msnaf_peptides.csv"
        export_peptides(
            counts_path=str(counts_path),
            refs_dir=str(refs_dir),
            genome_fasta_path=str(genome_fasta),
            output_path=str(output_path),
        )

        ok = True
        ok &= compare_text_files(
            expected_intermediates,
            tmpdir / "result" / "snaf_intermediates.tsv",
            "snaf_intermediates.tsv",
        )
        ok &= compare_text_files(
            expected_stats,
            tmpdir / "result" / "NeoJunction_statistics_maxmin.txt",
            "NeoJunction_statistics_maxmin.txt",
        )

        if keep_output:
            preserved = fixture_dir / "msnaf-test-output"
            if preserved.exists():
                shutil.rmtree(preserved)
            shutil.copytree(tmpdir / "result", preserved)
            print(f"saved outputs to {preserved}")

        return 0 if ok else 1
    finally:
        tmpdir_obj.cleanup()


def run_explicit_comparison(
    counts_path: Path,
    refs_dir: Path,
    genome_fasta: Path,
    expected_intermediates: Path,
    expected_stats: Path,
    keep_output_dir: Path | None = None,
) -> int:
    if not counts_path.exists():
        print(f"missing counts file: {counts_path}")
        return 2
    if not refs_dir.exists():
        print(f"missing refs dir: {refs_dir}")
        return 2
    if not expected_intermediates.exists():
        print(f"missing expected intermediates: {expected_intermediates}")
        return 2
    if not expected_stats.exists():
        print(f"missing expected stats: {expected_stats}")
        return 2
    if not genome_fasta.exists():
        print(f"missing genome fasta: {genome_fasta}")
        return 2

    tmpdir_obj = tempfile.TemporaryDirectory(prefix="msnaf_fixture_compare_")
    tmpdir = Path(tmpdir_obj.name)
    try:
        output_path = tmpdir / "result" / "msnaf_peptides.csv"
        export_peptides(
            counts_path=str(counts_path),
            refs_dir=str(refs_dir),
            genome_fasta_path=str(genome_fasta),
            output_path=str(output_path),
        )

        ok = True
        ok &= compare_text_files(
            expected_intermediates,
            tmpdir / "result" / "snaf_intermediates.tsv",
            "snaf_intermediates.tsv",
        )
        ok &= compare_text_files(
            expected_stats,
            tmpdir / "result" / "NeoJunction_statistics_maxmin.txt",
            "NeoJunction_statistics_maxmin.txt",
        )

        if keep_output_dir is not None:
            if keep_output_dir.exists():
                shutil.rmtree(keep_output_dir)
            shutil.copytree(tmpdir / "result", keep_output_dir)
            print(f"saved outputs to {keep_output_dir}")

        return 0 if ok else 1
    finally:
        tmpdir_obj.cleanup()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run msnaf on a fixture directory and compare legacy outputs against SNAF outputs."
    )
    parser.add_argument(
        "--fixture-dir",
        "--figure-dir",
        help="Directory containing counts.original.full.txt, snaf-data/, result/, and snaf_intermediates.tsv",
    )
    parser.add_argument(
        "--counts",
        "--count",
        help="Path to counts.original.full.txt",
    )
    parser.add_argument(
        "--refs",
        help="Path to the extracted SNAF reference directory",
    )
    parser.add_argument(
        "--genome-fasta",
        help="Path to the reference genome FASTA used for offline UTR sequence lookup",
    )
    parser.add_argument(
        "--expected-intermediates",
        help="Path to the expected legacy snaf_intermediates.tsv file",
    )
    parser.add_argument(
        "--expected-stats",
        help="Path to the expected NeoJunction_statistics_maxmin.txt file",
    )
    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Preserve the generated result directory for inspection.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.fixture_dir:
        return run_fixture_comparison(Path(args.fixture_dir), keep_output=args.keep_output)

    required = [
        args.counts,
        args.refs,
        args.genome_fasta,
        args.expected_intermediates,
        args.expected_stats,
    ]
    if any(value is None for value in required):
        print(
            "Either --fixture-dir or all of --counts, --refs, --genome-fasta, --expected-intermediates, and --expected-stats are required."
        )
        return 2

    keep_output_dir = None
    if args.keep_output:
        keep_output_dir = Path.cwd() / "msnaf-test-output"

    return run_explicit_comparison(
        counts_path=Path(args.counts),
        refs_dir=Path(args.refs),
        genome_fasta=Path(args.genome_fasta),
        expected_intermediates=Path(args.expected_intermediates),
        expected_stats=Path(args.expected_stats),
        keep_output_dir=keep_output_dir,
    )


if __name__ == "__main__":
    raise SystemExit(main())
