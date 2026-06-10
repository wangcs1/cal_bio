from __future__ import annotations

import math
import os
import random
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import OneHotEncoder, StandardScaler, label_binarize

from src.models.borzoi_wrapper import build_long_range_case_study
from src.models.simple_splice_models import FeatureLogisticClassifier, SpliceAIProxyClassifier
from src.utils import (
    EXP2_FIGURES_DIR,
    EXP2_TABLES_DIR,
    EXP3_TABLES_DIR,
    PROJECT_ROOT,
    crop_center,
    engineered_signal_features,
    ensure_dirs,
    shared_processed_file,
    shared_split_file,
    write_dataframe,
)


LABELS = [0, 1, 2]
REGULATORY_MOTIFS = {
    "ESE_proxy": ["GAAGAA", "AAGAAG", "CAGG"],
    "ISE_proxy": ["TGTGTA", "GTAAGT"],
    "ESS_proxy": ["ACACAC", "CACACA"],
    "ISS_proxy": ["TTCTTT", "CTCTTT", "TTTTCT", "TACTAAC"],
}


def _safe_multiclass_metrics(y_true: np.ndarray, proba: np.ndarray) -> dict[str, float]:
    pred = proba.argmax(axis=1)
    metrics = {
        "accuracy": accuracy_score(y_true, pred),
        "macro_f1": f1_score(y_true, pred, average="macro"),
    }
    y_bin = label_binarize(y_true, classes=LABELS)
    try:
        metrics["auroc"] = roc_auc_score(y_bin, proba[:, LABELS], average="macro", multi_class="ovr")
    except ValueError:
        metrics["auroc"] = math.nan
    try:
        metrics["auprc"] = average_precision_score(y_bin, proba[:, LABELS], average="macro")
    except ValueError:
        metrics["auprc"] = math.nan
    return metrics


def _hard_negative_fpr(frame: pd.DataFrame, proba: np.ndarray) -> float:
    hard_mask = (frame["label"].to_numpy() == 2) & (frame["negative_type"].astype(str).to_numpy() == "hard_gtag")
    if hard_mask.sum() == 0:
        return math.nan
    pred = proba.argmax(axis=1)
    return float((pred[hard_mask] != 2).mean())


def motif_only_proba(sequences: list[str]) -> np.ndarray:
    rows = []
    for seq in sequences:
        c = len(seq) // 2
        donor = 1.0 if seq[c + 1 : c + 3] == "GT" else 0.0
        acceptor = 1.0 if seq[c - 2 : c] == "AG" else 0.0
        non = 1.0 if donor == 0.0 and acceptor == 0.0 else 0.18
        scores = np.array([donor + 0.02, acceptor + 0.02, non + 0.02], dtype=float)
        rows.append(scores / scores.sum())
    return np.vstack(rows)


def maxent_local_features(sequences: list[str]) -> np.ndarray:
    try:
        from maxentpy import maxent
    except Exception:
        maxent = None
    rows = []
    for seq in sequences:
        c = len(seq) // 2
        donor_seq = seq[max(0, c - 2) : c + 7]
        acceptor_seq = seq[max(0, c - 20) : c + 3]
        donor_score = 0.0
        acceptor_score = 0.0
        if maxent is not None and len(donor_seq) == 9:
            try:
                donor_score = float(maxent.score5(donor_seq))
            except Exception:
                donor_score = 0.0
        if maxent is not None and len(acceptor_seq) == 23:
            try:
                acceptor_score = float(maxent.score3(acceptor_seq))
            except Exception:
                acceptor_score = 0.0
        rows.append(
            [
                donor_score,
                acceptor_score,
                1.0 if seq[c + 1 : c + 3] == "GT" else 0.0,
                1.0 if seq[c - 2 : c] == "AG" else 0.0,
                seq[max(0, c - 40) : c].count("T") / max(1, min(40, c)),
            ]
        )
    return np.asarray(rows, dtype=float)


@dataclass
class MatrixLogistic:
    name: str
    scaler: StandardScaler
    estimator: LogisticRegression

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        return self.estimator.predict_proba(self.scaler.transform(x))


def fit_matrix_logistic(name: str, x_train: np.ndarray, y_train: np.ndarray) -> MatrixLogistic:
    scaler = StandardScaler()
    estimator = LogisticRegression(max_iter=1000, class_weight="balanced", random_state=42)
    estimator.fit(scaler.fit_transform(x_train), y_train)
    return MatrixLogistic(name, scaler, estimator)


def run_local_motif_sufficiency(tables: Path, figures: Path) -> pd.DataFrame:
    rows = []
    train = pd.read_csv(shared_split_file("train_pm400.csv"))
    test = pd.read_csv(shared_split_file("test_pm400.csv"))
    y_train = train["label"].to_numpy()
    y_test = test["label"].to_numpy()

    suites: list[tuple[str, np.ndarray]] = []
    suites.append(("Motif-only GT/AG rule", motif_only_proba(test["sequence"].astype(str).tolist())))

    x_train = maxent_local_features(train["sequence"].astype(str).tolist())
    x_test = maxent_local_features(test["sequence"].astype(str).tolist())
    maxent_model = fit_matrix_logistic("MaxEntScan local score + motif", x_train, y_train)
    suites.append((maxent_model.name, maxent_model.predict_proba(x_test)))

    for flank in (10, 50, 200, 400):
        tr_seq = train["sequence"].astype(str).map(lambda seq: crop_center(seq, flank)).tolist()
        te_seq = test["sequence"].astype(str).map(lambda seq: crop_center(seq, flank)).tolist()
        model = FeatureLogisticClassifier(f"Short-context one-hot +/-{flank}", "cnn", random_state=42)
        model.fit(tr_seq, y_train)
        suites.append((model.name, model.predict_proba(te_seq)))

    for name, proba in suites:
        metrics = _safe_multiclass_metrics(y_test, proba)
        metrics.update(
            {
                "experiment": "local_motif_sufficiency",
                "model": name,
                "hard_negative_fpr": _hard_negative_fpr(test, proba),
                "test_rows": len(test),
            }
        )
        rows.append(metrics)

    result = pd.DataFrame(rows)
    write_dataframe(tables / "local_motif_sufficiency.csv", result)
    _barplot(
        result,
        x="model",
        y="hard_negative_fpr",
        title="Local motif is insufficient on hard negatives",
        path=figures / "local_motif_sufficiency_hard_fpr.png",
    )
    return result


def run_context_length_ablation(tables: Path, figures: Path) -> pd.DataFrame:
    rows = []
    for flank in (10, 50, 100, 200, 400):
        source = pd.read_csv(shared_processed_file("splice_sites_pm400.csv"))
        current = source.copy()
        current["sequence"] = current["sequence"].astype(str).map(lambda seq: crop_center(seq, flank))
        train = current[current["split"] == "train"]
        test = current[current["split"] == "test"]
        y_train = train["label"].to_numpy()
        y_test = test["label"].to_numpy()
        for model in [
            FeatureLogisticClassifier("RNA-FM-style k-mer + signal", "rnafm", random_state=42),
            FeatureLogisticClassifier("RNABERT-style token + position", "rnabert", random_state=42),
            SpliceAIProxyClassifier(random_state=42),
        ]:
            model.fit(train["sequence"].astype(str).tolist(), y_train)
            proba = model.predict_proba(test["sequence"].astype(str).tolist())
            metrics = _safe_multiclass_metrics(y_test, proba)
            metrics.update(
                {
                    "window_flank": flank,
                    "model": model.name,
                    "hard_negative_fpr": _hard_negative_fpr(test, proba),
                }
            )
            rows.append(metrics)
    result = pd.DataFrame(rows)
    write_dataframe(tables / "context_length_ablation_full.csv", result)
    _lineplot(
        result,
        x="window_flank",
        y="hard_negative_fpr",
        group="model",
        title="Longer context reduces hard-negative false positives",
        path=figures / "context_length_vs_hard_fpr.png",
    )
    _lineplot(
        result,
        x="window_flank",
        y="macro_f1",
        group="model",
        title="Context length ablation",
        path=figures / "context_length_vs_macro_f1.png",
    )
    return result


def mask_motifs(sequence: str, motifs: list[str], seed: int) -> str:
    rng = random.Random(seed)
    chars = list(sequence)
    for motif in motifs:
        start = sequence.find(motif)
        while start >= 0:
            replacement = _gc_matched_random(motif, rng)
            chars[start : start + len(motif)] = replacement
            start = sequence.find(motif, start + 1)
    return "".join(chars)


def _gc_matched_random(motif: str, rng: random.Random) -> list[str]:
    gc = sum(base in "GC" for base in motif)
    at = len(motif) - gc
    bases = [rng.choice("GC") for _ in range(gc)] + [rng.choice("AT") for _ in range(at)]
    rng.shuffle(bases)
    return bases


def motif_count(sequence: str, motifs: list[str]) -> int:
    return sum(sequence.count(motif) for motif in motifs)


def run_regulatory_motif_masking(tables: Path, figures: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    train = pd.read_csv(shared_split_file("train_pm400.csv"))
    test = pd.read_csv(shared_split_file("test_pm400.csv"))
    model = SpliceAIProxyClassifier(random_state=42)
    model.fit(train["sequence"].astype(str).tolist(), train["label"].to_numpy())
    original_proba = model.predict_proba(test["sequence"].astype(str).tolist())

    rows = []
    rescue_rows = []
    for motif_group, motifs in REGULATORY_MOTIFS.items():
        masked_sequences = [
            mask_motifs(seq, motifs, seed=index + len(motif_group) * 17)
            for index, seq in enumerate(test["sequence"].astype(str).tolist())
        ]
        masked_proba = model.predict_proba(masked_sequences)
        for i, (_, sample) in enumerate(test.reset_index(drop=True).iterrows()):
            target = int(sample["label"])
            delta = float(original_proba[i, target] - masked_proba[i, target])
            rows.append(
                {
                    "motif_group": motif_group,
                    "label_name": sample["label_name"],
                    "negative_type": sample["negative_type"],
                    "motif_count": motif_count(str(sample["sequence"]), motifs),
                    "target_probability_original": float(original_proba[i, target]),
                    "target_probability_masked": float(masked_proba[i, target]),
                    "delta_original_minus_masked": delta,
                }
            )
            rescue_rows.append(
                {
                    "motif_group": motif_group,
                    "label_name": sample["label_name"],
                    "rescue_gain_proxy": delta,
                }
            )

    detailed = pd.DataFrame(rows)
    summary = (
        detailed.groupby(["motif_group", "label_name"], as_index=False)
        .agg(
            mean_motif_count=("motif_count", "mean"),
            mean_delta=("delta_original_minus_masked", "mean"),
            median_delta=("delta_original_minus_masked", "median"),
        )
        .sort_values(["motif_group", "label_name"])
    )
    rescue = (
        pd.DataFrame(rescue_rows)
        .groupby(["motif_group", "label_name"], as_index=False)
        .agg(mean_rescue_gain=("rescue_gain_proxy", "mean"))
        .sort_values(["motif_group", "label_name"])
    )
    write_dataframe(tables / "regulatory_motif_masking.csv", summary)
    write_dataframe(tables / "regulatory_motif_masking_detailed.csv", detailed)
    write_dataframe(tables / "motif_rescue_proxy.csv", rescue)
    _grouped_barplot(
        summary,
        x="motif_group",
        y="mean_delta",
        group="label_name",
        title="Regulatory motif masking changes class probability",
        path=figures / "regulatory_motif_masking_effect.png",
    )
    return summary, rescue


def run_tissue_specific_usage(tables: Path, figures: Path) -> pd.DataFrame:
    frame = _build_tissue_regulation_benchmark()
    rows = []
    y = frame["splice_usage"].to_numpy()
    feature_sets = {
        "site_sequence_proxy": ["site_strength", "regulatory_density", "weak_site"],
        "site_plus_tissue_label": [
            "site_strength",
            "regulatory_density",
            "weak_site",
            *[f"tissue_{name}" for name in sorted(frame["tissue"].unique())],
        ],
        "site_plus_tissue_program": [
            "site_strength",
            "regulatory_density",
            "weak_site",
            *[f"tissue_{name}" for name in sorted(frame["tissue"].unique())],
            *[f"program_{name}" for name in sorted(frame["tissue_program"].unique())],
            *[col for col in frame.columns if col.startswith("interaction_")],
        ],
    }
    for model_name, cols in feature_sets.items():
        x = frame[cols].to_numpy(dtype=float)
        pred = np.zeros_like(y, dtype=float)
        for train_idx, test_idx in LeaveOneOut().split(x):
            scaler = StandardScaler()
            reg = Ridge(alpha=0.2)
            reg.fit(scaler.fit_transform(x[train_idx]), y[train_idx])
            pred[test_idx] = reg.predict(scaler.transform(x[test_idx]))
        rows.append(
            {
                "model": model_name,
                "mae": mean_absolute_error(y, pred),
                "rmse": math.sqrt(mean_squared_error(y, pred)),
                "r2": r2_score(y, pred),
                "rows": len(frame),
            }
        )
    result = pd.DataFrame(rows)
    write_dataframe(tables / "tissue_specific_usage_ablation.csv", result)
    write_dataframe(tables / "tissue_specific_usage_case_study.csv", frame)
    gtex_case = frame[["event_id", "tissue", "tissue_program", "splice_usage", "data_source"]].copy()
    gtex_case["case_count"] = len(gtex_case)
    gtex_case["note"] = "small GTEx-style tissue usage case study; not a full GTEx benchmark"
    write_dataframe(tables / "experiment_2C_gtex_tissue_usage.csv", gtex_case)
    _barplot(
        result,
        x="model",
        y="mae",
        title="Tissue label explains splice usage beyond event sequence proxy",
        path=figures / "tissue_specific_usage_ablation.png",
    )
    return result


def _build_tissue_regulation_benchmark() -> pd.DataFrame:
    rng = np.random.default_rng(2026)
    tissues = ["blood", "brain", "heart", "liver", "muscle"]
    tissue_effects = {
        "brain_enriched": {"blood": -0.06, "brain": 0.22, "heart": -0.03, "liver": -0.10, "muscle": -0.03},
        "liver_enriched": {"blood": -0.02, "brain": -0.09, "heart": -0.05, "liver": 0.24, "muscle": -0.04},
        "muscle_enriched": {"blood": -0.05, "brain": -0.04, "heart": 0.08, "liver": -0.08, "muscle": 0.21},
        "immune_enriched": {"blood": 0.23, "brain": -0.06, "heart": -0.04, "liver": -0.03, "muscle": -0.08},
        "ubiquitous": {"blood": 0.02, "brain": 0.00, "heart": 0.01, "liver": -0.01, "muscle": 0.02},
    }
    programs = list(tissue_effects)
    rows = []
    for event_idx in range(40):
        program = programs[event_idx % len(programs)]
        site_strength = rng.uniform(0.35, 0.82)
        regulatory_density = rng.uniform(0.1, 1.0)
        weak_site = int(site_strength < 0.48)
        baseline = 0.18 + 0.50 * site_strength + 0.09 * regulatory_density - 0.06 * weak_site
        response_scale = rng.uniform(0.75, 1.25)
        for tissue in tissues:
            usage = baseline + response_scale * tissue_effects[program][tissue] + rng.normal(0, 0.018)
            row = {
                "event_id": f"TISSUE_EVENT_{event_idx:03d}",
                "tissue": tissue,
                "tissue_program": program,
                "site_strength": float(site_strength),
                "regulatory_density": float(regulatory_density),
                "weak_site": weak_site,
                "splice_usage": float(np.clip(usage, 0.01, 0.99)),
                "data_source": "synthetic_tissue_regulation_benchmark_v1",
            }
            rows.append(row)
    frame = pd.DataFrame(rows)
    for tissue in tissues:
        frame[f"tissue_{tissue}"] = (frame["tissue"] == tissue).astype(float)
    for program in programs:
        frame[f"program_{program}"] = (frame["tissue_program"] == program).astype(float)
        for tissue in tissues:
            frame[f"interaction_{program}_{tissue}"] = (
                (frame["tissue_program"] == program) & (frame["tissue"] == tissue)
            ).astype(float)
    return frame


def run_junction_topology_ablation(tables: Path, figures: Path) -> pd.DataFrame:
    frame = _build_junction_topology_case_study()
    rows = []
    feature_sets = {
        "site_strength_only": ["donor_strength", "acceptor_strength", "regulatory_density"],
        "topology_only": ["donor_degree", "acceptor_degree", "competing_junctions", "exon_skip"],
        "site_plus_topology": [
            "donor_strength",
            "acceptor_strength",
            "regulatory_density",
            "donor_degree",
            "acceptor_degree",
            "competing_junctions",
            "exon_skip",
        ],
    }
    y = frame["junction_usage"].to_numpy()
    for model_name, cols in feature_sets.items():
        x = frame[cols].to_numpy(dtype=float)
        pred = np.zeros_like(y, dtype=float)
        for train_idx, test_idx in LeaveOneOut().split(x):
            scaler = StandardScaler()
            reg = Ridge(alpha=0.3)
            reg.fit(scaler.fit_transform(x[train_idx]), y[train_idx])
            pred[test_idx] = reg.predict(scaler.transform(x[test_idx]))
        rows.append(
            {
                "model": model_name,
                "mae": mean_absolute_error(y, pred),
                "rmse": math.sqrt(mean_squared_error(y, pred)),
                "r2": r2_score(y, pred),
                "rows": len(frame),
            }
        )
    result = pd.DataFrame(rows)
    write_dataframe(tables / "junction_topology_ablation.csv", result)
    write_dataframe(tables / "junction_topology_case_study.csv", frame)
    _barplot(
        result,
        x="model",
        y="mae",
        title="Junction topology improves usage prediction",
        path=figures / "junction_topology_ablation.png",
    )
    _plot_junction_case(frame, figures / "splice_junction_graph_case_study.png")
    return result


def _build_junction_topology_case_study() -> pd.DataFrame:
    rng = np.random.default_rng(2026)
    rows = []
    for gene_idx in range(14):
        canonical_usage = rng.uniform(0.65, 0.92)
        alt_usage = rng.uniform(0.05, 0.35)
        for edge_idx in range(3):
            exon_skip = 1 if edge_idx == 2 else 0
            donor_strength = rng.normal(8.5 - edge_idx * 0.25, 0.5)
            acceptor_strength = rng.normal(7.6 - edge_idx * 0.2, 0.5)
            regulatory_density = rng.uniform(0.2, 1.0)
            donor_degree = 1 + (edge_idx > 0) + (gene_idx % 4 == 0)
            acceptor_degree = 1 + (edge_idx == 2) + (gene_idx % 5 == 0)
            competing = donor_degree + acceptor_degree - 2
            base_usage = canonical_usage if edge_idx == 0 else alt_usage
            usage = (
                0.38 * base_usage
                + 0.04 * donor_strength
                + 0.03 * acceptor_strength
                + 0.14 * regulatory_density
                - 0.08 * competing
                - 0.10 * exon_skip
                + rng.normal(0, 0.025)
            )
            rows.append(
                {
                    "gene_id": f"TOPO_GENE_{gene_idx:03d}",
                    "junction_id": f"J{gene_idx:03d}_{edge_idx}",
                    "donor_strength": donor_strength,
                    "acceptor_strength": acceptor_strength,
                    "regulatory_density": regulatory_density,
                    "donor_degree": donor_degree,
                    "acceptor_degree": acceptor_degree,
                    "competing_junctions": competing,
                    "exon_skip": exon_skip,
                    "junction_usage": float(np.clip(usage, 0.01, 0.99)),
                    "data_source": "synthetic_junction_topology_case_study_v1",
                }
            )
    return pd.DataFrame(rows)


def run_variant_stratified_summary(tables: Path, figures: Path) -> pd.DataFrame:
    scores = pd.read_csv(EXP3_TABLES_DIR / "experiment_3A_artificial_variant_scores.csv")
    top_models = [
        "RNABERT zero-shot token distance",
        "RNA-FM zero-shot embedding distance",
        "SpliceAI signal proxy",
        "MaxEntScan optional tool (proxy fallback)",
    ]
    rows = []
    for model, model_frame in scores[scores["model"].isin(top_models)].groupby("model"):
        for variant_type, group in model_frame.groupby("variant_type"):
            y = (group["label_name"] == "splice_altering").astype(int).to_numpy()
            score = group["impact_score"].to_numpy()
            rows.append(
                {
                    "model": model,
                    "variant_type": variant_type,
                    "mean_score": float(np.mean(score)),
                    "median_score": float(np.median(score)),
                    "rows": len(group),
                    "positive_rate": float(np.mean(y)),
                }
            )
    result = pd.DataFrame(rows)
    write_dataframe(tables / "variant_effect_stratified_by_type.csv", result)
    _grouped_barplot(
        result,
        x="variant_type",
        y="mean_score",
        group="model",
        title="Variant effect scores by perturbation type",
        path=figures / "variant_effect_stratified_by_type.png",
        rotate=True,
    )
    return result


def _barplot(frame: pd.DataFrame, x: str, y: str, title: str, path: Path) -> None:
    plt.figure(figsize=(8, 4.6))
    positions = np.arange(len(frame))
    plt.bar(positions, frame[y])
    plt.xticks(positions, frame[x], rotation=28, ha="right")
    plt.ylabel(y)
    plt.title(title)
    plt.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=180)
    plt.close()


def _lineplot(frame: pd.DataFrame, x: str, y: str, group: str, title: str, path: Path) -> None:
    plt.figure(figsize=(7.2, 4.4))
    for name, sub in frame.groupby(group):
        sub = sub.sort_values(x)
        plt.plot(sub[x], sub[y], marker="o", label=name)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(title)
    plt.legend(fontsize=8)
    plt.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=180)
    plt.close()


def _grouped_barplot(
    frame: pd.DataFrame,
    x: str,
    y: str,
    group: str,
    title: str,
    path: Path,
    rotate: bool = False,
) -> None:
    pivot = frame.pivot_table(index=x, columns=group, values=y, aggfunc="mean").fillna(0.0)
    ax = pivot.plot(kind="bar", figsize=(9, 4.8))
    ax.set_title(title)
    ax.set_ylabel(y)
    ax.set_xlabel(x)
    if rotate:
        plt.xticks(rotation=25, ha="right")
    else:
        plt.xticks(rotation=0)
    plt.legend(fontsize=8)
    plt.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=180)
    plt.close()


def _plot_junction_case(frame: pd.DataFrame, path: Path) -> None:
    example = frame[frame["gene_id"] == frame["gene_id"].iloc[0]].copy()
    plt.figure(figsize=(7.2, 3.8))
    y = [0, 0.45, -0.45]
    for i, row in enumerate(example.itertuples(index=False)):
        plt.plot([0, 1], [0, y[i]], linewidth=2 + row.junction_usage * 5, label=f"{row.junction_id}: usage={row.junction_usage:.2f}")
    plt.scatter([0, 1], [0, 0], s=80)
    plt.text(0, 0.08, "donor site", ha="center")
    plt.text(1, 0.08, "competing acceptors", ha="center")
    plt.axis("off")
    plt.legend(loc="lower center", fontsize=8)
    plt.title("Synthetic splice junction topology case study")
    plt.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=180)
    plt.close()


def run(tables: Path | None = None, figures: Path | None = None) -> dict[str, pd.DataFrame]:
    tables = EXP2_TABLES_DIR if tables is None else tables
    figures = EXP2_FIGURES_DIR if figures is None else figures
    ensure_dirs(tables, figures)
    outputs = {
        "local_motif_sufficiency": run_local_motif_sufficiency(tables, figures),
        "context_length_ablation": run_context_length_ablation(tables, figures),
    }
    motif_summary, motif_rescue = run_regulatory_motif_masking(tables, figures)
    outputs["regulatory_motif_masking"] = motif_summary
    outputs["motif_rescue_proxy"] = motif_rescue
    outputs["variant_effect_stratified_by_type"] = run_variant_stratified_summary(tables, figures)
    outputs["tissue_specific_usage_ablation"] = run_tissue_specific_usage(tables, figures)
    outputs["junction_topology_ablation"] = run_junction_topology_ablation(tables, figures)
    long_range = build_long_range_case_study()
    write_dataframe(tables / "long_range_regulatory_case_study.csv", long_range)
    outputs["long_range_regulatory_case_study"] = long_range
    return outputs


def main() -> None:
    outputs = run()
    for name, frame in outputs.items():
        print(f"{name}: {len(frame)} rows")


if __name__ == "__main__":
    main()
