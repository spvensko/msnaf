import tempfile
import unittest
from unittest.mock import patch
import pathlib
import sys
import types

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

if "xmltodict" not in sys.modules:
    xmltodict_module = types.ModuleType("xmltodict")
    xmltodict_module.parse = lambda content: {}
    sys.modules["xmltodict"] = xmltodict_module

from msnaf.core import ReferenceData, collect_records, filter_junctions, iter_peptide_records, load_counts_matrix


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


if __name__ == "__main__":
    unittest.main()
