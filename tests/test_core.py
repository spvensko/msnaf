import pathlib
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

try:
    import anndata as ad
except ModuleNotFoundError:
    class FakeAnnData:
        def __init__(self, X, obs=None, var=None):
            self.X = np.array(X)
            self.obs = obs if obs is not None else pd.DataFrame()
            self.var = var if var is not None else pd.DataFrame()
            self.obs_names = self.obs.index
            self.var_names = self.var.index
            self.shape = self.X.shape

        def to_df(self):
            return pd.DataFrame(self.X, index=self.obs.index, columns=self.var.index)

    ad_module = types.ModuleType("anndata")
    ad_module.AnnData = FakeAnnData
    sys.modules["anndata"] = ad_module
    import anndata as ad

from msnaf.core import GenomeReference, ReferenceData, collect_records, filter_junctions, iter_peptide_records, load_counts_matrix, retrieve_seq_from_genome, utr_junction


class CoreTests(unittest.TestCase):
    def test_load_counts_matrix_strips_trailing_coord_and_deduplicates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            counts_path = f"{tmpdir}/counts.original.full.txt"
            with open(counts_path, "w") as handle:
                handle.write("\tS1\n")
                handle.write("ENSG1:E1.1-E2.1=chr1:1-2(+)\t10\n")
                handle.write("ENSG1:E1.1-E2.1=chr1:1-2(+)\t99\n")
                handle.write("ENSG2:E2.1-E3.1=chr2:2-3(+)\t5\n")

            df = load_counts_matrix(counts_path)

        self.assertEqual(list(df.index), ["ENSG1:E1.1-E2.1", "ENSG2:E2.1-E3.1"])
        self.assertEqual(df.loc["ENSG1:E1.1-E2.1", "S1"], 10)

    def test_filter_junctions_maxmin_matches_snaf_threshold_logic(self):
        tumor = pd.DataFrame(
            {"S1": [30, 10, 2], "S2": [26, 12, 2]},
            index=["uid_a", "uid_b", "uid_c"],
        )
        gtex = ad.AnnData(
            X=np.array([[1, 1], [5, 5], [0, 0]]),
            obs=pd.DataFrame(index=["uid_a", "uid_b", "uid_c"]),
            var=pd.DataFrame({"tissue": ["t1", "t2"]}, index=["g1", "g2"]),
        )
        gtex.obs["mean"] = np.array([1.0, 5.0, 0.0])
        reference = ReferenceData(
            dict_exon_coords={},
            dict_exonlist={},
            dict_fa={},
            dict_start_codon={},
            phase_inferer_gtf_dict={},
            genome=GenomeReference(records={}),
            adata=gtex,
            adata_gtex=gtex,
            t_min=20,
            n_max=3,
            normal_cutoff=5,
            tumor_cutoff=20,
            normal_prevalance_cutoff=0.01,
            tumor_prevalance_cutoff=0.1,
        )

        query = filter_junctions(tumor, reference, add_control=None, filter_mode="maxmin", not_in_db=False)

        self.assertEqual(set(query.valid), {"uid_a"})
        self.assertEqual(set(query.invalid), {"uid_b", "uid_c"})
        self.assertTrue(bool(query.cond_df.loc["uid_a", "S1"]))
        self.assertFalse(bool(query.cond_df.loc["uid_b", "S1"]))
        self.assertEqual(query.stats_filename, "NeoJunction_statistics_maxmin.txt")
        self.assertEqual(list(query.stats_df.columns), ["max", "mean", "diff", "cond"])

    def test_iter_peptide_records_returns_cross_junction_peptide_and_cds(self):
        records = list(
            iter_peptide_records(
                de_facto_first="ATGAAA",
                second="CCCGGG",
                ks=[3],
                phase=0,
                evidences=(),
                coord="chr1:1-2(+)",
            )
        )

        peptides = {(record["peptide"], record["coding_sequence"]) for record in records}

        self.assertIn(("KPG", "AAACCCGGG"), peptides)

    def test_collect_records_flattens_translation_results(self):
        reference = ReferenceData(
            dict_exon_coords={},
            dict_exonlist={},
            dict_fa={},
            dict_start_codon={},
            phase_inferer_gtf_dict={},
            genome=GenomeReference(records={}),
            adata=None,
            adata_gtex=None,
            t_min=20,
            n_max=3,
            normal_cutoff=5,
            tumor_cutoff=20,
            normal_prevalance_cutoff=0.01,
            tumor_prevalance_cutoff=0.1,
        )

        def fake_translate(uid, reference, strict=False):
            del reference, strict
            return [{"uid": uid, "coord": "c", "peptide": "P", "coding_sequence": "ATG", "peptide_context": "XPY"}]

        with patch("msnaf.core.translate_uid", fake_translate):
            rows = collect_records(["u1", "u2"], reference, strict=True)

        self.assertEqual(
            rows,
            [
                {"uid": "u1", "coord": "c", "peptide": "P", "coding_sequence": "ATG", "peptide_context": "XPY"},
                {"uid": "u2", "coord": "c", "peptide": "P", "coding_sequence": "ATG", "peptide_context": "XPY"},
            ],
        )

    def test_export_peptides_writes_legacy_outputs(self):
        from msnaf import core

        reference = ReferenceData(
            dict_exon_coords={},
            dict_exonlist={},
            dict_fa={},
            dict_start_codon={},
            phase_inferer_gtf_dict={},
            genome=GenomeReference(records={}),
            adata=None,
            adata_gtex=None,
            t_min=20,
            n_max=3,
            normal_cutoff=5,
            tumor_cutoff=20,
            normal_prevalance_cutoff=0.01,
            tumor_prevalance_cutoff=0.1,
        )
        query = core.QueryResult(
            valid=["u1"],
            invalid=[],
            cond_df=pd.DataFrame(),
            stats_df=pd.DataFrame({"max": [10], "mean": [0], "diff": [10], "cond": [True]}, index=["u1"]),
            stats_filename="NeoJunction_statistics_maxmin.txt",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            counts_path = f"{tmpdir}/counts.txt"
            output_path = f"{tmpdir}/result/msnaf_peptides.csv"
            with open(counts_path, "w") as handle:
                handle.write("\tS1\nu1\t1\n")

            with patch("msnaf.core.load_reference_data", return_value=reference), \
                 patch("msnaf.core.filter_junctions", return_value=query), \
                 patch(
                     "msnaf.core.collect_records",
                     return_value=[
                         {
                             "uid": "u1",
                             "coord": "chr1:1-2(+)",
                             "peptide": "PEPTIDE",
                             "coding_sequence": "ATGAAA",
                             "peptide_context": "XXPEPTIDEYY",
                         }
                     ],
                 ):
                df = core.export_peptides(
                    counts_path=counts_path,
                    refs_dir=tmpdir,
                    genome_fasta_path=f"{tmpdir}/genome.fa",
                    output_path=output_path,
                )

            self.assertEqual(df.shape[0], 1)
            self.assertTrue(pathlib.Path(output_path).exists())
            self.assertTrue(pathlib.Path(f"{tmpdir}/result/snaf_intermediates.tsv").exists())
            self.assertTrue(pathlib.Path(f"{tmpdir}/result/NeoJunction_statistics_maxmin.txt").exists())

    def test_retrieve_seq_from_genome_returns_placeholder_for_missing_chromosome(self):
        reference = ReferenceData(
            dict_exon_coords={},
            dict_exonlist={},
            dict_fa={},
            dict_start_codon={},
            phase_inferer_gtf_dict={},
            genome=GenomeReference(records={}),
            adata=None,
            adata_gtex=None,
            t_min=20,
            n_max=3,
            normal_cutoff=5,
            tumor_cutoff=20,
            normal_prevalance_cutoff=0.01,
            tumor_prevalance_cutoff=0.1,
        )

        seq = retrieve_seq_from_genome(reference, "chr1", 1, 100)

        self.assertEqual(seq, "#" * 10)

    def test_utr_junction_reads_from_local_genome_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = pathlib.Path(tmpdir) / "genome.fa"
            fasta_path.write_text(">chr1\nACGTACGTACGT\n", encoding="utf-8")
            genome = GenomeReference.open(str(fasta_path))
            try:
                reference = ReferenceData(
                    dict_exon_coords={},
                    dict_exonlist={},
                    dict_fa={},
                    dict_start_codon={},
                    phase_inferer_gtf_dict={},
                    genome=genome,
                    adata=None,
                    adata_gtex=None,
                    t_min=20,
                    n_max=3,
                    normal_cutoff=5,
                    tumor_cutoff=20,
                    normal_prevalance_cutoff=0.01,
                    tumor_prevalance_cutoff=0.1,
                )

                seq = utr_junction("4", "ENSG1", "+", "chr1", "site1", reference, seq_len=4)
            finally:
                genome.close()

        self.assertEqual(seq, "ACGT")


if __name__ == "__main__":
    unittest.main()
