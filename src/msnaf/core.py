from __future__ import annotations

from ast import literal_eval
from dataclasses import dataclass
import os
import re
import sys
from typing import Iterable

import anndata as ad
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


@dataclass
class GenomeReference:
    records: object

    @classmethod
    def open(cls, fasta_path: str) -> "GenomeReference":
        return cls(records=SeqIO.index(fasta_path, "fasta"))

    def fetch(self, chrom: str, start: int, end: int) -> str:
        names = [chrom]
        if chrom.startswith("chr"):
            names.append(chrom[3:])
        else:
            names.append(f"chr{chrom}")
        for name in names:
            if name in self.records:
                return str(self.records[name].seq[start - 1 : end]).upper()
        raise KeyError(chrom)

    def close(self) -> None:
        close = getattr(self.records, "close", None)
        if close is not None:
            close()


@dataclass
class ReferenceData:
    dict_exon_coords: dict
    dict_exonlist: dict
    dict_fa: dict
    dict_start_codon: dict
    phase_inferer_gtf_dict: dict
    genome: GenomeReference
    adata: ad.AnnData
    adata_gtex: ad.AnnData
    t_min: int
    n_max: int
    normal_cutoff: int
    tumor_cutoff: int
    normal_prevalance_cutoff: float
    tumor_prevalance_cutoff: float


@dataclass
class QueryResult:
    valid: list[str]
    invalid: list[str]
    cond_df: pd.DataFrame
    stats_df: pd.DataFrame
    stats_filename: str


def log(message: str) -> None:
    print(f"[msnaf] {message}", file=sys.stderr, flush=True)


def load_counts_matrix(path: str, sep: str = "\t") -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index = [item.split("=")[0] for item in df.index]
    df = df.loc[~df.index.duplicated(), :]
    return df


def load_controls(refs_dir: str, use_tcga_control: bool = True) -> dict | None:
    if not use_tcga_control:
        return None
    tcga_path = os.path.join(refs_dir, "controls", "tcga_matched_control_junction_count.h5ad")
    if not os.path.exists(tcga_path):
        return None
    return {"tcga_control": ad.read_h5ad(tcga_path, backed="r")}


def subset_controls(add_control: dict | None, tested_junctions: set[str], filter_mode: str) -> dict | None:
    if add_control is None:
        return None
    if filter_mode == "maxmin":
        return add_control
    subsetted = {}
    for cohort_name, control in add_control.items():
        if isinstance(control, pd.DataFrame):
            control = control.loc[~control.index.duplicated(), :]
            control = control.loc[list(set(control.index).intersection(tested_junctions)), :]
        elif isinstance(control, ad.AnnData):
            selected_obs = [obs_name for obs_name in control.obs_names if obs_name in tested_junctions]
            control = control[selected_obs, :]
        else:
            raise TypeError("control must be either a pandas DataFrame or an AnnData object")
        subsetted[cohort_name] = control
    return subsetted


def export_peptides(
    counts_path: str,
    refs_dir: str,
    genome_fasta_path: str,
    output_path: str,
    strict: bool = False,
    filter_mode: str = "maxmin",
    not_in_db: bool = False,
    use_tcga_control: bool = True,
    t_min: int = 20,
    n_max: int = 3,
    normal_cutoff: int = 5,
    tumor_cutoff: int = 20,
    normal_prevalance_cutoff: float = 0.01,
    tumor_prevalance_cutoff: float = 0.1,
) -> pd.DataFrame:
    log(f"loading counts matrix from {counts_path}")
    df = load_counts_matrix(counts_path)
    log(f"loaded counts matrix with {df.shape[0]} junctions across {df.shape[1]} samples")
    log("loading additional control databases")
    add_control = subset_controls(
        load_controls(refs_dir, use_tcga_control=use_tcga_control),
        set(df.index),
        filter_mode,
    )
    if add_control is None:
        log("no additional controls loaded")
    else:
        log(f"loaded additional controls: {', '.join(sorted(add_control.keys()))}")
    log(f"loading reference bundle from {refs_dir}")
    reference = load_reference_data(
        df=df,
        refs_dir=refs_dir,
        genome_fasta_path=genome_fasta_path,
        filter_mode=filter_mode,
        t_min=t_min,
        n_max=n_max,
        normal_cutoff=normal_cutoff,
        tumor_cutoff=tumor_cutoff,
        normal_prevalance_cutoff=normal_prevalance_cutoff,
        tumor_prevalance_cutoff=tumor_prevalance_cutoff,
    )
    try:
        log(f"filtering junctions with mode={filter_mode}")
        query = filter_junctions(
            junction_count_matrix=df,
            reference=reference,
            add_control=add_control,
            filter_mode=filter_mode,
            not_in_db=not_in_db,
        )
        log(f"filtering complete: {len(query.valid)} valid, {len(query.invalid)} filtered out")
        log("translating retained junctions")
        records = collect_records(query.valid, reference, strict=strict)
        log(f"translation complete: {len(records)} peptide rows")
        result = pd.DataFrame.from_records(
            records,
            columns=["uid", "coord", "peptide", "coding_sequence", "peptide_context"],
        )
        if result.shape[0] != 0:
            result = result.drop_duplicates().sort_values(["uid", "peptide", "coding_sequence"])
        output_dir = os.path.dirname(os.path.abspath(output_path))
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        log(f"writing outputs to {output_dir}")
        result.to_csv(output_path, index=False)
        query.stats_df.to_csv(os.path.join(output_dir, query.stats_filename), sep="\t")
        result.loc[:, ["coord", "peptide", "coding_sequence", "peptide_context"]].to_csv(
            os.path.join(output_dir, "snaf_intermediates.tsv"),
            header=False,
            index=False,
        )
        log("run complete")
        return result
    finally:
        reference.genome.close()


def collect_records(valid_uids: list[str], reference: ReferenceData, strict: bool = False) -> list[dict]:
    records = []
    total = len(valid_uids)
    for index, uid in enumerate(valid_uids, start=1):
        if index == 1 or index % 250 == 0 or index == total:
            log(f"translating junction {index}/{total}")
        records.extend(translate_uid(uid, reference, strict=strict))
    return records


def load_reference_data(
    df: pd.DataFrame,
    refs_dir: str,
    genome_fasta_path: str,
    filter_mode: str,
    t_min: int = 20,
    n_max: int = 3,
    normal_cutoff: int = 5,
    tumor_cutoff: int = 20,
    normal_prevalance_cutoff: float = 0.01,
    tumor_prevalance_cutoff: float = 0.1,
) -> ReferenceData:
    exon_table = os.path.join(refs_dir, "Alt91_db", "Hs_Ensembl_exon_add_col.txt")
    transcript_db = os.path.join(refs_dir, "Alt91_db", "mRNA-ExonIDs.txt")
    fasta = os.path.join(refs_dir, "Alt91_db", "Hs_gene-seq-2000_flank.fa")
    gtf = os.path.join(refs_dir, "Alt91_db", "Homo_sapiens.GRCh38.91.gtf")
    start_codon_path = os.path.join(refs_dir, "Alt91_db", "df_start_codon.txt")
    gtex_db = os.path.join(refs_dir, "controls", "GTEx_junction_counts.h5ad")

    log("loading exon coordinate table")
    dict_exon_coords = exon_coords_to_dict(exon_table)
    log("loading transcript exon database")
    dict_exonlist = construct_dict_exonlist(transcript_db)
    log("loading splice flank FASTA")
    dict_fa = fasta_to_dict(fasta)
    log("loading GTF phase annotations")
    phase_inferer_gtf_dict = process_gtf(gtf)
    log(f"opening reference genome FASTA {genome_fasta_path}")
    genome = GenomeReference.open(genome_fasta_path)

    log("loading start codon table")
    df_start_codon = pd.read_csv(start_codon_path, sep="\t", index_col=0)
    df_start_codon["start_codon"] = [literal_eval(item) for item in df_start_codon["start_codon"]]
    dict_start_codon = df_start_codon["start_codon"].to_dict()

    log("loading GTEx control database")
    adata, adata_gtex = configure_controls(
        df=df,
        gtex_db=gtex_db,
        filter_mode=filter_mode,
        t_min=t_min,
        n_max=n_max,
        normal_cutoff=normal_cutoff,
        tumor_cutoff=tumor_cutoff,
        normal_prevalance_cutoff=normal_prevalance_cutoff,
        tumor_prevalance_cutoff=tumor_prevalance_cutoff,
    )
    return ReferenceData(
        dict_exon_coords=dict_exon_coords,
        dict_exonlist=dict_exonlist,
        dict_fa=dict_fa,
        dict_start_codon=dict_start_codon,
        phase_inferer_gtf_dict=phase_inferer_gtf_dict,
        genome=genome,
        adata=adata,
        adata_gtex=adata_gtex,
        t_min=t_min,
        n_max=n_max,
        normal_cutoff=normal_cutoff,
        tumor_cutoff=tumor_cutoff,
        normal_prevalance_cutoff=normal_prevalance_cutoff,
        tumor_prevalance_cutoff=tumor_prevalance_cutoff,
    )


def configure_controls(
    df: pd.DataFrame,
    gtex_db: str,
    filter_mode: str,
    t_min: int,
    n_max: int,
    normal_cutoff: int,
    tumor_cutoff: int,
    normal_prevalance_cutoff: float,
    tumor_prevalance_cutoff: float,
) -> tuple[ad.AnnData, ad.AnnData]:
    tested_junctions = set(df.index)
    adata = ad.read_h5ad(gtex_db, backed="r")
    if filter_mode != "maxmin":
        selected_obs = [obs_name for obs_name in adata.obs_names if obs_name in tested_junctions]
        adata = adata[selected_obs, :]
    if "mean" not in adata.obs.columns:
        adata.obs["mean"] = np.array(adata.X.mean(axis=1)).squeeze()
    if "total_count" not in adata.var.columns:
        adata.var["total_count"] = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
    adata_gtex = adata

    return adata, adata_gtex


def filter_junctions(
    junction_count_matrix: pd.DataFrame,
    reference: ReferenceData,
    add_control: dict | None = None,
    filter_mode: str = "maxmin",
    not_in_db: bool = False,
) -> QueryResult:
    if filter_mode == "prevalance":
        valid, invalid, cond_df, stats_df = _filter_prevalance(
            junction_count_matrix=junction_count_matrix,
            reference=reference,
            add_control=add_control,
            not_in_db=not_in_db,
        )
        stats_filename = "NeoJunction_statistics_prevalance.txt"
    elif filter_mode == "maxmin":
        valid, invalid, cond_df, stats_df = _filter_maxmin(
            junction_count_matrix=junction_count_matrix,
            reference=reference,
            add_control=add_control,
            not_in_db=not_in_db,
        )
        stats_filename = "NeoJunction_statistics_maxmin.txt"
    else:
        raise ValueError(f"unsupported filter mode: {filter_mode}")
    return QueryResult(
        valid=valid,
        invalid=invalid,
        cond_df=cond_df,
        stats_df=stats_df,
        stats_filename=stats_filename,
    )


def _filter_prevalance(
    junction_count_matrix: pd.DataFrame,
    reference: ReferenceData,
    add_control: dict | None,
    not_in_db: bool,
) -> tuple[list[str], list[str], pd.DataFrame]:
    df = pd.DataFrame(index=junction_count_matrix.index)
    prevalance_tumor = np.count_nonzero(
        (junction_count_matrix > reference.tumor_cutoff).values, axis=1
    ) / junction_count_matrix.shape[1]
    df["prevalance_tumor"] = prevalance_tumor

    prevalance_normal = np.count_nonzero(
        (reference.adata_gtex.X > reference.normal_cutoff).toarray(), axis=1
    ) / reference.adata_gtex.shape[1]
    normal_dict = {uid: value for uid, value in zip(reference.adata_gtex.obs_names, prevalance_normal)}
    df["prevalance_normal"] = df.index.map(normal_dict).fillna(value=0)
    df["cond"] = (
        (df["prevalance_tumor"] > reference.tumor_prevalance_cutoff)
        & (df["prevalance_normal"] < reference.normal_prevalance_cutoff)
    )
    valid = df.loc[df["cond"]].index.tolist()
    if not_in_db:
        valid = [uid for uid in valid if not uid_is_in_db(uid, reference.dict_exonlist)]

    if add_control is not None:
        for cohort_name, control in add_control.items():
            if isinstance(control, pd.DataFrame):
                prevalance_add = np.count_nonzero(
                    (control > reference.normal_cutoff).values, axis=1
                ) / control.shape[1]
                normal_add_dict = {uid: value for uid, value in zip(control.index, prevalance_add)}
            elif isinstance(control, ad.AnnData):
                prevalance_add = np.count_nonzero(
                    (control.X > reference.normal_cutoff).toarray(), axis=1
                ) / control.shape[1]
                normal_add_dict = {uid: value for uid, value in zip(control.obs_names, prevalance_add)}
            else:
                raise TypeError("control must be either a pandas DataFrame or an AnnData object")
            df["prevalance_normal_add"] = df.index.map(normal_add_dict).fillna(value=0)
            df["cond_add"] = (
                (df["prevalance_tumor"] > reference.tumor_prevalance_cutoff)
                & (df["prevalance_normal_add"] < reference.normal_prevalance_cutoff)
            )
            valid = list(set(valid).intersection(df.loc[df["cond_add"]].index.tolist()))
            df.rename(
                columns={
                    "prevalance_normal_add": f"prevalance_normal_{cohort_name}",
                    "cond_add": f"cond_{cohort_name}",
                },
                inplace=True,
            )

    invalid = list(set(junction_count_matrix.index).difference(set(valid)))
    valid_set = set(valid)
    placeholder = pd.DataFrame(
        index=junction_count_matrix.index,
        data={"placeholder": junction_count_matrix.index.map(lambda uid: uid in valid_set).values},
    )
    first_half = pd.concat([placeholder] * junction_count_matrix.shape[1], axis=1)
    first_half.columns = junction_count_matrix.columns
    cond_df = first_half & (junction_count_matrix > reference.tumor_cutoff)
    return valid, invalid, cond_df, df


def _filter_maxmin(
    junction_count_matrix: pd.DataFrame,
    reference: ReferenceData,
    add_control: dict | None,
    not_in_db: bool,
) -> tuple[list[str], list[str], pd.DataFrame, pd.DataFrame]:
    df = pd.DataFrame(index=junction_count_matrix.index, data={"max": junction_count_matrix.max(axis=1).values})
    df_to_write = [df.copy()]
    junction_to_mean = reference.adata_gtex.obs.loc[
        reference.adata_gtex.obs_names.isin(junction_count_matrix.index), "mean"
    ].to_dict()
    df["mean"] = df.index.map(junction_to_mean).fillna(value=0)
    df["diff"] = df["max"] - df["mean"]
    df["cond"] = (df["mean"] < reference.n_max) & (df["diff"] > reference.t_min)
    df_to_write[0] = df.copy()
    valid = df.loc[df["cond"]].index.tolist()
    if not_in_db:
        valid = [uid for uid in valid if not uid_is_in_db(uid, reference.dict_exonlist)]

    mean_add_list = []
    if add_control is not None:
        for cohort_name, control in add_control.items():
            if isinstance(control, pd.DataFrame):
                junction_to_mean_add = control.mean(axis=1).to_dict()
            elif isinstance(control, ad.AnnData):
                if "mean" in control.obs.columns:
                    junction_to_mean_add = control.obs.loc[
                        control.obs_names.isin(junction_count_matrix.index), "mean"
                    ].to_dict()
                else:
                    selected_obs = [obs_name for obs_name in control.obs_names if obs_name in junction_count_matrix.index]
                    junction_to_mean_add = np.asarray(control[selected_obs, :].X.mean(axis=1)).squeeze()
                    junction_to_mean_add = {
                        uid: value for uid, value in zip(selected_obs, junction_to_mean_add)
                    }
            else:
                raise TypeError("control must be either a pandas DataFrame or an AnnData object")
            df["mean_add"] = df.index.map(junction_to_mean_add).fillna(value=0)
            df["diff_add"] = df["max"] - df["mean_add"]
            df["cond_add"] = (df["mean_add"] < reference.n_max) & (df["diff_add"] > reference.t_min)
            mean_add_list.append(df["mean_add"])
            valid = list(set(valid).intersection(df.loc[df["cond_add"]].index.tolist()))
            tmp = df.copy()
            tmp.drop(columns=["mean", "diff", "cond"], inplace=True)
            tmp.rename(columns=lambda column: f"{column}_{cohort_name}", inplace=True)
            df_to_write.append(tmp)

    invalid = list(set(junction_count_matrix.index).difference(set(valid)))
    gtex_df = pd.concat([df["mean"]] * junction_count_matrix.shape[1], axis=1)
    gtex_df.columns = junction_count_matrix.columns
    diff_df = junction_count_matrix - gtex_df
    cond_df = (gtex_df < reference.n_max) & (diff_df > reference.t_min)
    for mean_add in mean_add_list:
        add_df = pd.concat([mean_add] * junction_count_matrix.shape[1], axis=1)
        add_df.columns = junction_count_matrix.columns
        diff_add = junction_count_matrix - add_df
        cond_df = cond_df & (add_df < reference.n_max) & (diff_add > reference.t_min)
    stats_df = pd.concat(df_to_write, axis=1)
    return valid, invalid, cond_df, stats_df


def uid_is_in_db(uid: str, dict_exonlist: dict) -> bool:
    ensg = uid.split(":")[0]
    exons = ":".join(uid.split(":")[1:])
    if "_" in exons or "U" in exons or "ENSG" in exons or "I" in exons:
        return False
    exonlist = dict_exonlist.get(ensg)
    if exonlist is None:
        return False
    exonstring = "|".join(exonlist)
    e1, e2 = exons.split("-")
    e1 = re.escape(e1)
    e2 = re.escape(e2)
    return bool(
        re.search(rf"^{e1}\|{e2}\|", exonstring)
        or re.search(rf"\|{e1}\|{e2}$", exonstring)
        or re.search(rf"\|{e1}\|{e2}\|", exonstring)
    )


def translate_uid(uid: str, reference: ReferenceData, strict: bool = False, ks: tuple[int, ...] = (9, 10)) -> list[dict]:
    event_type = detect_type(uid)
    if event_type == "invalid":
        return []
    junction = retrieve_junction_seq(uid, reference)
    coord = uid_to_coord(uid, reference.dict_exon_coords)
    if "$" in junction or "*" in junction or "#" in junction or "unknown" in coord:
        return []

    first, second = junction.split(",")
    ensg = uid.split(":")[0]
    strand = coord.split("(")[1].rstrip(")")
    if strand == "+":
        coord_first_exon_last_base = coord.split(":")[1].split("-")[0]
    else:
        coord_first_exon_last_base = coord.split(":")[1].split("(")[0].split("-")[1]

    support_phases_dict = {}
    for pssc in reference.dict_start_codon.get(ensg, []):
        supports = get_support_phase(
            phase_inferer_gtf_dict=reference.phase_inferer_gtf_dict,
            ensg=ensg,
            coord_first_exon_last_base=coord_first_exon_last_base,
            pssc=pssc,
            strand=strand,
            length_first=len(first),
        )
        for phase, pssc_value, enst, strand_value in supports:
            support_phases_dict.setdefault(phase, []).append((pssc_value, enst, strand_value))

    records = []
    for phase in [0, 1, 2]:
        evidences = tuple(support_phases_dict.get(phase, []))
        if strict and len(evidences) == 0:
            continue
        de_facto_first = first[phase:]
        for record in iter_peptide_records(
            de_facto_first=de_facto_first,
            second=second,
            ks=ks,
            phase=phase,
            evidences=evidences,
            coord=coord,
        ):
            records.append(
                {
                    "uid": uid,
                    "coord": record["coord"],
                    "peptide": record["peptide"],
                    "coding_sequence": record["coding_sequence"],
                    "peptide_context": record["peptide_context"],
                }
            )
    return records


def detect_type(uid: str) -> str:
    valid_pattern = re.compile(r"^ENSG\d+:.+?-.+")
    if not re.search(valid_pattern, uid):
        return "invalid"
    if len(re.findall("ENSG", uid)) == 2:
        return "trans_splicing"
    if "U" in uid:
        return "utr_event"
    if "_" in uid:
        subexon12 = uid.split(":")[1]
        subexon1, subexon2 = subexon12.split("-")
        if "I" in subexon12:
            return "novel_exon"
        if "_" in subexon1 and "_" in subexon2:
            return "alt5_alt3"
        if "_" in subexon1:
            return "alt5"
        if "_" in subexon2:
            return "alt3"
        return "invalid"
    if "I" in uid:
        return "intron_retention"
    if re.search(r"^ENSG\d+:E\d+\.\d+-E\d+\.\d+$", uid):
        e1, e2 = uid.split(":")[1].split("-")
        e1_int, e1_frac = int(e1.split(".")[0][1:]), int(e1.split(".")[1])
        e2_int, e2_frac = int(e2.split(".")[0][1:]), int(e2.split(".")[1])
        if e1 == e2:
            return "invalid"
        if e1_int > e2_int or (e1_int == e2_int and e1_frac > e2_frac):
            return "invalid"
        return "ordinary"
    return "invalid"


def retrieve_junction_seq(uid: str, reference: ReferenceData) -> str:
    ensid = uid.split(":")[0]
    subexon1, subexon2 = ":".join(uid.split(":")[1:]).split("-")
    code = is_consecutive(subexon1, subexon2)
    seq1 = subexon_tran(subexon1, ensid, "site1", code, reference)
    seq2 = subexon_tran(subexon2, ensid, "site2", code, reference)
    return ",".join([seq1, seq2])


def iter_peptide_records(
    de_facto_first: str,
    second: str,
    ks: Iterable[int],
    phase: int,
    evidences: tuple,
    coord: str,
):
    extra = len(de_facto_first) % 3
    num = len(de_facto_first) // 3
    aa_first = str(Seq(de_facto_first).translate(to_stop=True))
    if len(aa_first) != num:
        return
    if extra == 0:
        continue_second = second
    elif extra == 1:
        continue_second = de_facto_first[-1] + second
    else:
        continue_second = de_facto_first[-2:] + second
    aa_second = str(Seq(continue_second).translate(to_stop=True))
    if len(aa_second) == 0:
        return
    for k in ks:
        second_most = min(k, len(aa_second))
        first_most = len(aa_first)
        for n_from_second in range(second_most, 0, -1):
            n_from_first = k - n_from_second
            if n_from_first == 0 and extra == 0:
                continue
            if n_from_first > first_most:
                continue
            if n_from_first > 0:
                pep = aa_first[-n_from_first:] + aa_second[:n_from_second]
                pep_context = aa_first[-n_from_first - 5 :] + aa_second[: n_from_second + 10]
            else:
                pep = aa_second[:n_from_second]
                pep_context = aa_second[: n_from_second + 10]
            cds1 = str(Seq(de_facto_first))[len(de_facto_first) - len(pep) * 3 - n_from_first * 3 :]
            cds2 = str(Seq(continue_second))[: len(pep) * 3 + n_from_second * 3]
            raw_cds = cds1 + cds2
            actual_cds = ""
            for i in range(len(raw_cds) - (len(pep) * 3)):
                for j in [0, 1, 2]:
                    potential_cds = raw_cds[j + i : j + i + len(pep) * 3 + 1]
                    if str(Seq(potential_cds).translate(to_stop=True)) == pep:
                        actual_cds = potential_cds
            yield {
                "mer": k,
                "peptide": pep,
                "extra": extra,
                "n_from_first": n_from_first,
                "phase": phase,
                "evidences": evidences,
                "coord": coord,
                "coding_sequence": actual_cds,
                "peptide_context": pep_context,
            }


def fasta_to_dict(path: str) -> dict:
    result = {}
    with open(path, "r") as handle:
        for title, seq in SimpleFastaParser(handle):
            result[title.split("|")[0]] = [title.split("|")[1], title.split("|")[2], title.split("|")[3], seq]
    return result


def exon_coords_to_dict(path: str) -> dict:
    result = {}
    with open(path, "r") as handle:
        next(handle)
        for line in handle:
            items = line.split("\t")
            coords = (items[2], items[3], items[4], items[5], items[10].rstrip("\n"))
            result.setdefault(items[0], {})[items[1]] = coords
    return result


def construct_dict_exonlist(transcript_db: str) -> dict:
    df = pd.read_csv(
        transcript_db,
        sep="\t",
        header=None,
        names=["EnsGID", "EnsTID", "EnsPID", "Exons"],
    )
    return df.groupby(by="EnsGID")["Exons"].apply(lambda values: values.tolist()).to_dict()


def process_gtf(gtf: str) -> dict:
    gtf_dict = {}
    with open(gtf, "r") as handle:
        for line in handle:
            try:
                _, _, record_type, start, end, _, _, _, attrs = line.rstrip("\n").split("\t")
            except ValueError:
                continue
            if record_type == "gene":
                ensg = attrs.split(";")[0].split(" ")[1].strip('"')
                gtf_dict[ensg] = {}
            elif record_type == "transcript":
                enst = attrs.split(";")[2].split(" ")[2].strip('"')
                gtf_dict[ensg][enst] = []
            elif record_type == "exon":
                gtf_dict[ensg][enst].append((int(start), int(end)))
    return gtf_dict


def get_support_phase(
    phase_inferer_gtf_dict: dict,
    ensg: str,
    coord_first_exon_last_base: str,
    pssc: int,
    strand: str,
    length_first: int,
) -> list[tuple]:
    from bisect import bisect, bisect_left

    coord_first_exon_last_base = int(coord_first_exon_last_base)
    pssc = int(pssc)
    length_first = int(length_first)
    all_trans = phase_inferer_gtf_dict[ensg]
    supports = []
    if strand == "+":
        for enst, tran in all_trans.items():
            flat = []
            [flat.extend(list(exon)) for exon in tran]
            pssc_insert_pos = bisect(flat, pssc)
            coord_first_exon_first_base = coord_first_exon_last_base - length_first + 1
            junction_insert_pos = bisect(flat, coord_first_exon_first_base)
            if junction_insert_pos % 2 != 1 or pssc_insert_pos % 2 != 1:
                continue
            if junction_insert_pos == pssc_insert_pos:
                n_bases = coord_first_exon_first_base - pssc + 1
            elif junction_insert_pos > pssc_insert_pos:
                start_index = (pssc_insert_pos - 1) // 2
                end_index = (junction_insert_pos - 1) // 2
                n_bases = 0
                for i, exon in enumerate(tran):
                    if i == start_index:
                        n_bases += exon[1] - pssc + 1
                    elif i == end_index:
                        n_bases += coord_first_exon_first_base - exon[0] + 1
                    elif start_index < i < end_index:
                        n_bases += exon[1] - exon[0] + 1
            else:
                continue
            remainder = n_bases % 3
            phase = 0 if remainder == 1 else 2 if remainder == 2 else 1
            supports.append((phase, pssc, enst, strand))
    elif strand == "-":
        for enst, tran in all_trans.items():
            tran = sorted(tran, key=lambda exon: exon[0])
            flat = []
            [flat.extend(list(exon)) for exon in tran]
            pssc_insert_pos = bisect_left(flat, pssc)
            coord_first_exon_first_base = coord_first_exon_last_base + length_first - 1
            junction_insert_pos = bisect_left(flat, coord_first_exon_first_base)
            if junction_insert_pos % 2 != 1 or pssc_insert_pos % 2 != 1:
                continue
            if junction_insert_pos == pssc_insert_pos:
                n_bases = pssc - coord_first_exon_first_base + 1
            elif junction_insert_pos < pssc_insert_pos:
                start_index = (pssc_insert_pos - 1) // 2
                end_index = (junction_insert_pos - 1) // 2
                n_bases = 0
                for i, exon in enumerate(tran):
                    if i == start_index:
                        n_bases += pssc - exon[0] + 1
                    elif i == end_index:
                        n_bases += exon[1] - coord_first_exon_first_base + 1
                    elif end_index < i < start_index:
                        n_bases += exon[1] - exon[0] + 1
            else:
                continue
            remainder = n_bases % 3
            phase = 0 if remainder == 1 else 2 if remainder == 2 else 1
            supports.append((phase, pssc, enst, strand))
    return supports


def query_from_dict_fa(abs_start: int, abs_end: int, ensid: str, strand: str, dict_fa: dict) -> str:
    if strand == "+":
        start = int(dict_fa[ensid][1])
        seq = dict_fa[ensid][3]
        start_index = int(abs_start) - start + 2000
        end_index = int(abs_end) - start + 1 + 2000
        return seq[start_index:end_index]

    start = int(dict_fa[ensid][1])
    seq_reverse = dict_fa[ensid][3]
    seq_forward = str(Seq(seq_reverse).reverse_complement())
    start_index = int(abs_start) - start + 2000
    end_index = int(abs_end) - start + 1 + 2000
    exon_seq = seq_forward[start_index:end_index]
    return str(Seq(exon_seq).reverse_complement())


def utr_attrs(ensid: str, dict_exon_coords: dict) -> tuple[str, str]:
    attrs = next(iter(dict_exon_coords[ensid].values()))
    return attrs[0], attrs[1]


def retrieve_seq_from_genome(reference: ReferenceData, chr_: str, start: int, end: int) -> str:
    try:
        return reference.genome.fetch(chr_, start, end)
    except Exception:
        return "#" * 10


def utr_junction(site: str, ensid: str, strand: str, chr_: str, flag: str, reference: ReferenceData, seq_len: int = 100) -> str:
    del ensid
    if flag == "site1" and strand == "+":
        other_site = int(site) - seq_len + 1
        return retrieve_seq_from_genome(reference, chr_, int(other_site), int(site))
    if flag == "site1" and strand == "-":
        other_site = int(site) + seq_len - 1
        exon_seq = retrieve_seq_from_genome(reference, chr_, int(site), int(other_site))
        return str(Seq(exon_seq).reverse_complement())
    if flag == "site2" and strand == "+":
        other_site = int(site) + seq_len - 1
        return retrieve_seq_from_genome(reference, chr_, int(site), int(other_site))
    other_site = int(site) - seq_len + 1
    exon_seq = retrieve_seq_from_genome(reference, chr_, int(other_site), int(site))
    return str(Seq(exon_seq).reverse_complement())


def is_consecutive(subexon1: str, subexon2: str) -> int:
    code = 0
    match = re.search(r"^E(\d{1,3})\.(\d{1,3})$", subexon1)
    if match:
        pattern = re.compile(rf"^E{match.group(1)}\.{int(match.group(2)) + 1}$")
        if re.search(pattern, subexon2):
            code = 1
    return code


def subexon_tran(subexon: str, ensid: str, flag: str, code: int, reference: ReferenceData) -> str:
    try:
        attrs = reference.dict_exon_coords[ensid][subexon]
        if code == 1 and flag == "site2":
            if attrs[1] == "+":
                return query_from_dict_fa(int(attrs[2]) + 1, attrs[3], ensid, attrs[1], reference.dict_fa)
            return query_from_dict_fa(attrs[2], int(attrs[3]) - 1, ensid, attrs[1], reference.dict_fa)
        return query_from_dict_fa(attrs[2], attrs[3], ensid, attrs[1], reference.dict_fa)
    except KeyError:
        if ":" in subexon:
            fusion_ensid, fusion_exon = subexon.split(":")
            if "_" in fusion_exon:
                suffix = fusion_exon.split("_")[1]
                actual_exon = fusion_exon.split("_")[0]
                attrs = reference.dict_exon_coords[fusion_ensid][actual_exon]
                if attrs[1] == "+":
                    return query_from_dict_fa(suffix, attrs[3], fusion_ensid, attrs[1], reference.dict_fa)
                return query_from_dict_fa(attrs[2], suffix, fusion_ensid, attrs[1], reference.dict_fa)
            try:
                attrs = reference.dict_exon_coords[fusion_ensid][fusion_exon]
            except KeyError:
                return "*" * 10
            return query_from_dict_fa(attrs[2], attrs[3], fusion_ensid, attrs[1], reference.dict_fa)

        parts = subexon.split("_")
        if len(parts) == 1:
            return "*" * 10
        actual_exon, suffix = parts
        try:
            attrs = reference.dict_exon_coords[ensid][actual_exon]
        except KeyError:
            chr_utr, strand_utr = utr_attrs(ensid, reference.dict_exon_coords)
            return utr_junction(suffix, ensid, strand_utr, chr_utr, flag, reference)
        if flag == "site2":
            if attrs[1] == "+":
                return query_from_dict_fa(suffix, attrs[3], ensid, attrs[1], reference.dict_fa)
            return query_from_dict_fa(attrs[2], suffix, ensid, attrs[1], reference.dict_fa)
        if attrs[1] == "+":
            return query_from_dict_fa(attrs[2], suffix, ensid, attrs[1], reference.dict_fa)
        return query_from_dict_fa(suffix, attrs[3], ensid, attrs[1], reference.dict_fa)


def uid_to_coord(uid: str, dict_exon_coords: dict) -> str:
    parts = uid.split(":")
    if len(parts) == 2:
        ensg, exons = parts
    else:
        ensg = parts[0]
        exons = ":".join(parts[1:])
    first, second = exons.split("-")

    if "_" in first:
        actual_exon, trailing = first.split("_")
        try:
            attrs = dict_exon_coords[ensg][actual_exon]
        except KeyError:
            if "U" in actual_exon:
                attrs = dict_exon_coords[ensg][next(iter(dict_exon_coords[ensg].keys()))]
                chrom, strand, start_coord = attrs[0], attrs[1], trailing
            else:
                chrom, strand, start_coord = "unknown", "unknown", "unknown"
        else:
            chrom, strand, start_coord = attrs[0], attrs[1], trailing
    else:
        try:
            attrs = dict_exon_coords[ensg][first]
        except KeyError:
            chrom, strand, start_coord = "unknown", "unknown", "unknown"
        else:
            chrom, strand = attrs[0], attrs[1]
            start_coord = attrs[3] if strand == "+" else attrs[2]

    if "_" in second:
        actual_exon, trailing = second.split("_")
        try:
            dict_exon_coords[ensg][actual_exon]
        except KeyError:
            if "U" in actual_exon or "ENSG" in actual_exon:
                end_coord = trailing
            else:
                end_coord = "unknown"
        else:
            end_coord = trailing
    else:
        try:
            attrs = dict_exon_coords[ensg][second]
        except KeyError:
            if "ENSG" in second:
                ensg_second, actual_exon_second = second.split(":")
                try:
                    attrs = dict_exon_coords[ensg_second][actual_exon_second]
                except KeyError:
                    end_coord = "unknown"
                else:
                    end_coord = attrs[2] if strand == "+" else attrs[3]
            else:
                end_coord = "unknown"
        else:
            end_coord = attrs[2] if strand == "+" else attrs[3]

    if strand == "+":
        return f"{chrom}:{start_coord}-{end_coord}({strand})"
    return f"{chrom}:{end_coord}-{start_coord}({strand})"
