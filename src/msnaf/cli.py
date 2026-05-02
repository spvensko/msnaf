import argparse

from .core import export_peptides


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Filter junctions like SNAF, then export splice-derived peptides and coding "
            "sequences without any HLA or pMHC prediction."
        )
    )
    parser.add_argument("--counts", required=True, help="Path to counts.original.full.txt or equivalent.")
    parser.add_argument("--refs", required=True, help="Path to the extracted SNAF reference bundle.")
    parser.add_argument("--genome-fasta", required=True, help="Path to the reference genome FASTA used for offline UTR sequence lookup.")
    parser.add_argument("--output", required=True, help="CSV path to write.")
    parser.add_argument("--strict", action="store_true", help="Require start-codon support.")
    parser.add_argument(
        "--filter-mode",
        choices=["maxmin", "prevalance"],
        default="maxmin",
        help="Tumor-vs-control filter used by SNAF. Default: %(default)s.",
    )
    parser.add_argument(
        "--not-in-db",
        action="store_true",
        help="Also drop junctions already documented in Ensembl transcripts.",
    )
    parser.add_argument(
        "--skip-tcga-control",
        action="store_true",
        help="Do not auto-load refs/controls/tcga_matched_control_junction_count.h5ad.",
    )
    parser.add_argument("--t-min", type=int, default=20, help="SNAF maxmin tumor delta threshold.")
    parser.add_argument("--n-max", type=int, default=3, help="SNAF maxmin normal mean threshold.")
    parser.add_argument("--normal-cutoff", type=int, default=5, help="SNAF prevalence normal cutoff.")
    parser.add_argument("--tumor-cutoff", type=int, default=20, help="SNAF prevalence tumor cutoff.")
    parser.add_argument(
        "--normal-prevalance-cutoff",
        type=float,
        default=0.01,
        help="SNAF prevalence normal fraction cutoff.",
    )
    parser.add_argument(
        "--tumor-prevalance-cutoff",
        type=float,
        default=0.1,
        help="SNAF prevalence tumor fraction cutoff.",
    )
    parser.add_argument(
        "--species",
        choices=["human", "mouse", "auto"],
        default="auto",
        help="Species for reference database selection. 'auto' detects from junction IDs. Default: %(default)s.",
    )
    return parser


def main():
    args = build_parser().parse_args()
    export_peptides(
        counts_path=args.counts,
        refs_dir=args.refs,
        genome_fasta_path=args.genome_fasta,
        output_path=args.output,
        strict=args.strict,
        filter_mode=args.filter_mode,
        not_in_db=args.not_in_db,
        use_tcga_control=not args.skip_tcga_control,
        t_min=args.t_min,
        n_max=args.n_max,
        normal_cutoff=args.normal_cutoff,
        tumor_cutoff=args.tumor_cutoff,
        normal_prevalance_cutoff=args.normal_prevalance_cutoff,
        tumor_prevalance_cutoff=args.tumor_prevalance_cutoff,
        species=args.species,
    )


if __name__ == "__main__":
    main()
