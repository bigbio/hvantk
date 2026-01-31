"""End-to-end ancestry inference pipeline.

This module provides the main pipeline orchestration for ancestry inference,
combining all steps: merging, filtering, PCA, and classification into a
single end-to-end workflow.

Classes
-------
PipelineConfig
    Configuration dataclass for all pipeline parameters.
AncestryInferenceResult
    Container for complete ancestry inference results.

Functions
---------
run_ancestry_inference
    Run the complete ancestry inference pipeline.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import hail as hl
import pandas as pd

from hvantk.ancestry.constants import (
    ANCESTRY_COL,
    ANCESTRY_PROB_COL,
    DEFAULT_HWE_P,
    DEFAULT_LD_R2,
    DEFAULT_LD_WINDOW,
    DEFAULT_MAX_AF,
    DEFAULT_MIN_AF,
    DEFAULT_MIN_CALL_RATE,
    DEFAULT_MIN_PROB,
    DEFAULT_N_CV_FOLDS,
    DEFAULT_N_ESTIMATORS,
    DEFAULT_N_PCS,
    DEFAULT_N_PCS_CLASSIFY,
    DEFAULT_RANDOM_SEED,
    KNOWN_ANCESTRY_COL,
    MIN_SAMPLES_PER_POP,
    PREDICTED_ANCESTRY_COL,
    SOURCE_COL,
)
from hvantk.ancestry.classify import (
    ClassificationResult,
    get_query_samples,
    get_training_samples,
    predict_ancestry,
    train_classifier,
)
from hvantk.ancestry.filter import (
    filter_variants_for_ancestry,
    ld_prune_variants,
)
from hvantk.ancestry.merge import merge_matrixtables
from hvantk.ancestry.pca import PCAResult, compute_pca

logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    """Configuration for ancestry inference pipeline.

    All parameters have sensible defaults based on standard practices
    in population genetics ancestry inference.

    Attributes
    ----------
    min_af : float
        Minimum allele frequency for variant inclusion. Default: 0.01.
    max_af : float
        Maximum allele frequency for variant inclusion. Default: 0.99.
    min_call_rate : float
        Minimum call rate for variant inclusion. Default: 0.98.
    apply_hwe_filter : bool
        Whether to apply Hardy-Weinberg equilibrium filter. Default: False.
    hwe_p_threshold : float
        HWE p-value threshold (variants below are removed). Default: 1e-6.
    ld_r2 : float
        LD r-squared threshold for pruning. Default: 0.2.
    ld_window : int
        LD window size in base pairs. Default: 500,000.
    skip_ld_pruning : bool
        Whether to skip LD pruning step. Default: False.
    n_pcs : int
        Number of principal components to compute. Default: 20.
    n_pcs_classify : int
        Number of PCs to use for classification. Default: 10.
    n_estimators : int
        Number of trees in Random Forest. Default: 100.
    min_prob : float
        Minimum probability for ancestry assignment. Default: 0.75.
    random_seed : int
        Random seed for reproducibility. Default: 42.
    validate_model : bool
        Whether to perform cross-validation. Default: True.
    n_cv_folds : int
        Number of cross-validation folds. Default: 5.
    min_samples_per_pop : int
        Minimum samples per population for training. Default: 10.
    min_shared_variants : int, optional
        Minimum shared variants required. Default: None (use constant).
    checkpoint_path : str, optional
        Path for intermediate checkpoints. Default: None.
    overwrite_checkpoints : bool
        Whether to overwrite existing checkpoints. Default: False.

    Example
    -------
    >>> config = PipelineConfig(
    ...     min_af=0.05,
    ...     min_prob=0.80,
    ...     n_pcs_classify=15,
    ... )
    >>> result = run_ancestry_inference(query_mt, ref_mt, config=config)
    """

    # Variant filtering
    min_af: float = DEFAULT_MIN_AF
    max_af: float = DEFAULT_MAX_AF
    min_call_rate: float = DEFAULT_MIN_CALL_RATE
    apply_hwe_filter: bool = False
    hwe_p_threshold: float = DEFAULT_HWE_P

    # LD pruning
    ld_r2: float = DEFAULT_LD_R2
    ld_window: int = DEFAULT_LD_WINDOW
    skip_ld_pruning: bool = False

    # PCA
    n_pcs: int = DEFAULT_N_PCS

    # Classification
    n_pcs_classify: int = DEFAULT_N_PCS_CLASSIFY
    n_estimators: int = DEFAULT_N_ESTIMATORS
    min_prob: float = DEFAULT_MIN_PROB
    random_seed: int = DEFAULT_RANDOM_SEED

    # Validation
    validate_model: bool = True
    n_cv_folds: int = DEFAULT_N_CV_FOLDS
    min_samples_per_pop: int = MIN_SAMPLES_PER_POP

    # Merge settings
    min_shared_variants: Optional[int] = None

    # Checkpointing
    checkpoint_path: Optional[str] = None
    overwrite_checkpoints: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary.

        Returns
        -------
        Dict[str, Any]
            Dictionary representation of all config parameters.
        """
        return {
            # Variant filtering
            "min_af": self.min_af,
            "max_af": self.max_af,
            "min_call_rate": self.min_call_rate,
            "apply_hwe_filter": self.apply_hwe_filter,
            "hwe_p_threshold": self.hwe_p_threshold,
            # LD pruning
            "ld_r2": self.ld_r2,
            "ld_window": self.ld_window,
            "skip_ld_pruning": self.skip_ld_pruning,
            # PCA
            "n_pcs": self.n_pcs,
            # Classification
            "n_pcs_classify": self.n_pcs_classify,
            "n_estimators": self.n_estimators,
            "min_prob": self.min_prob,
            "random_seed": self.random_seed,
            # Validation
            "validate_model": self.validate_model,
            "n_cv_folds": self.n_cv_folds,
            "min_samples_per_pop": self.min_samples_per_pop,
            # Merge settings
            "min_shared_variants": self.min_shared_variants,
            # Checkpointing
            "checkpoint_path": self.checkpoint_path,
            "overwrite_checkpoints": self.overwrite_checkpoints,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "PipelineConfig":
        """Create config from dictionary.

        Parameters
        ----------
        d : Dict[str, Any]
            Dictionary with config parameters.

        Returns
        -------
        PipelineConfig
            New config instance.
        """
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class AncestryInferenceResult:
    """Complete results from ancestry inference pipeline.

    This container holds all outputs from the ancestry inference pipeline,
    including predictions, PC scores, the trained model, and metadata.

    Attributes
    ----------
    predictions : hl.Table
        Hail Table with ancestry predictions for all samples.
    pc_scores : hl.Table
        Hail Table with PC scores for all samples.
    model : Any
        Trained Random Forest classifier.
    loadings : hl.Table
        PCA loadings for projection of new samples.
    eigenvalues : List[float]
        PCA eigenvalues.
    config : PipelineConfig
        Configuration used for the pipeline.
    pipeline_stats : Dict[str, Any]
        Statistics collected during pipeline execution.
    classification_result : ClassificationResult
        Full classification result including validation metrics.
    pca_result : PCAResult
        Full PCA result.

    Example
    -------
    >>> result = run_ancestry_inference(query_mt, ref_mt)
    >>> predictions_df = result.get_predictions_df()
    >>> print(predictions_df['predicted_ancestry'].value_counts())
    """

    # Core results
    predictions: hl.Table
    pc_scores: hl.Table
    model: Any  # RandomForestClassifier
    loadings: hl.Table
    eigenvalues: List[float]

    # Metadata
    config: PipelineConfig
    pipeline_stats: Dict[str, Any]

    # Full result objects
    classification_result: ClassificationResult
    pca_result: PCAResult

    # Cached DataFrames
    _predictions_df: Optional[pd.DataFrame] = field(default=None, repr=False)
    _scores_df: Optional[pd.DataFrame] = field(default=None, repr=False)

    def get_predictions_df(self) -> pd.DataFrame:
        """Get predictions as pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with all prediction columns.
        """
        if self._predictions_df is None:
            self._predictions_df = self.predictions.to_pandas()
        return self._predictions_df

    def get_scores_df(self) -> pd.DataFrame:
        """Get PC scores as pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with PC scores for all samples.
        """
        if self._scores_df is None:
            self._scores_df = self.pc_scores.to_pandas()
        return self._scores_df

    def get_query_predictions(self) -> pd.DataFrame:
        """Get predictions for query samples only.

        Returns
        -------
        pd.DataFrame
            DataFrame with predictions for query samples.
        """
        df = self.get_predictions_df()
        return df[df[SOURCE_COL] == "query"].copy()

    def get_reference_predictions(self) -> pd.DataFrame:
        """Get predictions for reference samples only.

        Returns
        -------
        pd.DataFrame
            DataFrame with predictions for reference samples.
        """
        df = self.get_predictions_df()
        return df[df[SOURCE_COL] == "reference"].copy()

    def get_accuracy(self) -> Optional[float]:
        """Get cross-validation accuracy.

        Returns
        -------
        float or None
            Accuracy if validation was performed.
        """
        return self.classification_result.get_accuracy()

    def variance_explained(self) -> List[float]:
        """Get proportion of variance explained by each PC.

        Returns
        -------
        List[float]
            Variance explained by each PC.
        """
        return self.pca_result.variance_explained()

    def prediction_summary(self) -> Dict[str, int]:
        """Get summary of predictions by ancestry.

        Returns
        -------
        Dict[str, int]
            Count of samples per predicted ancestry.
        """
        df = self.get_query_predictions()
        return df[PREDICTED_ANCESTRY_COL].value_counts().to_dict()

    def annotate_matrixtable(
        self,
        mt: hl.MatrixTable,
    ) -> hl.MatrixTable:
        """Annotate a MatrixTable with ancestry predictions.

        Adds column annotations for predicted ancestry, probability,
        and PC scores.

        Parameters
        ----------
        mt : hl.MatrixTable
            MatrixTable to annotate.

        Returns
        -------
        hl.MatrixTable
            Annotated MatrixTable.
        """
        # Create annotation table keyed by sample
        pred_ht = self.predictions

        # Join predictions to MT columns
        mt = mt.annotate_cols(
            **{
                PREDICTED_ANCESTRY_COL: pred_ht[mt.col_key][PREDICTED_ANCESTRY_COL],
                ANCESTRY_PROB_COL: pred_ht[mt.col_key][ANCESTRY_PROB_COL],
            }
        )

        # Add PC scores
        scores_ht = self.pc_scores
        n_pcs = self.config.n_pcs_classify
        pc_annotations = {
            f"PC{i + 1}": scores_ht[mt.col_key][f"PC{i + 1}"]
            for i in range(n_pcs)
        }
        mt = mt.annotate_cols(**pc_annotations)

        return mt

    def plot_pca(
        self,
        pc_x: int = 1,
        pc_y: int = 2,
        show_query_as_undefined: bool = False,
        **kwargs: Any,
    ) -> Any:
        """Plot PCA scatter plot.

        Parameters
        ----------
        pc_x : int, optional
            PC for x-axis. Default: 1.
        pc_y : int, optional
            PC for y-axis. Default: 2.
        show_query_as_undefined : bool, optional
            Show query samples as "Undefined" instead of their predicted ancestry.
            Default: False.
        **kwargs
            Additional arguments passed to plot_pca_scatter().

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        from hvantk.ancestry.plot import plot_pca_scatter

        return plot_pca_scatter(
            self.get_predictions_df(),
            pc_x=pc_x,
            pc_y=pc_y,
            show_query_as_undefined=show_query_as_undefined,
            **kwargs,
        )

    def plot_pca_panel(
        self,
        show_query_as_undefined: bool = False,
        **kwargs: Any,
    ) -> Any:
        """Plot two-panel PCA with full view and zoomed cluster.

        Parameters
        ----------
        show_query_as_undefined : bool, optional
            Show query samples as "Undefined" instead of their predicted ancestry.
            Default: False.
        **kwargs
            Additional arguments passed to plot_pca_panel().

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        from hvantk.ancestry.plot import plot_pca_panel

        return plot_pca_panel(
            self.get_predictions_df(),
            show_query_as_undefined=show_query_as_undefined,
            **kwargs,
        )

    def plot_ancestry_proportions(self, **kwargs: Any) -> Any:
        """Plot ancestry proportions bar chart.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_ancestry_proportions().

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        from hvantk.ancestry.plot import plot_ancestry_proportions

        return plot_ancestry_proportions(self.get_predictions_df(), **kwargs)

    def plot_probability_distribution(self, **kwargs: Any) -> Any:
        """Plot probability distribution histogram.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to plot_probability_distribution().

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        from hvantk.ancestry.plot import plot_probability_distribution

        return plot_probability_distribution(self.get_predictions_df(), **kwargs)

    def plot_variance_explained(self, n_pcs: int = 10, **kwargs: Any) -> Any:
        """Plot variance explained by principal components.

        Parameters
        ----------
        n_pcs : int, optional
            Number of PCs to show. Default: 10.
        **kwargs
            Additional arguments passed to plot_variance_explained().

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.
        """
        from hvantk.ancestry.plot import plot_variance_explained

        return plot_variance_explained(self.eigenvalues, n_pcs=n_pcs, **kwargs)

    def plot_confusion_matrix(self, normalize: bool = True, **kwargs: Any) -> Any:
        """Plot confusion matrix from cross-validation.

        Parameters
        ----------
        normalize : bool, optional
            Whether to normalize the matrix. Default: True.
        **kwargs
            Additional arguments passed to plot_confusion_matrix().

        Returns
        -------
        matplotlib.figure.Figure or None
            The generated figure, or None if validation was not performed.
        """
        from hvantk.ancestry.plot import plot_confusion_matrix

        labels = self.classification_result.confusion_matrix_labels
        if labels is None:
            logger.warning("No confusion matrix labels available (validation skipped?)")
            return None

        y_true, y_pred = labels
        return plot_confusion_matrix(
            y_true,
            y_pred,
            labels=self.classification_result.classes,
            normalize=normalize,
            **kwargs,
        )

    def generate_report(
        self,
        output_path: Union[str, Path],
        title: str = "Ancestry Inference Report",
        **kwargs: Any,
    ) -> Path:
        """Generate comprehensive HTML report.

        Parameters
        ----------
        output_path : str or Path
            Output file path for HTML report.
        title : str, optional
            Report title. Default: "Ancestry Inference Report".
        **kwargs
            Additional arguments passed to generate_ancestry_report().

        Returns
        -------
        Path
            Path to generated report.
        """
        from hvantk.ancestry.report import generate_ancestry_report

        return generate_ancestry_report(
            result=self,
            output_path=output_path,
            title=title,
            **kwargs,
        )

    def save(
        self,
        output_path: Union[str, Path],
        save_model: bool = True,
        save_loadings: bool = True,
    ) -> Dict[str, str]:
        """Save all results to disk.

        Parameters
        ----------
        output_path : str or Path
            Base output directory.
        save_model : bool, optional
            Whether to save the trained model. Default: True.
        save_loadings : bool, optional
            Whether to save PCA loadings. Default: True.

        Returns
        -------
        Dict[str, str]
            Dictionary of output paths.
        """
        import pickle

        output_dir = Path(output_path)
        output_dir.mkdir(parents=True, exist_ok=True)

        saved_paths = {}

        # Save predictions table
        predictions_path = str(output_dir / "predictions.ht")
        self.predictions.write(predictions_path, overwrite=True)
        saved_paths["predictions"] = predictions_path
        logger.info(f"Saved predictions: {predictions_path}")

        # Save PC scores table
        scores_path = str(output_dir / "pc_scores.ht")
        self.pc_scores.write(scores_path, overwrite=True)
        saved_paths["pc_scores"] = scores_path
        logger.info(f"Saved PC scores: {scores_path}")

        # Export predictions TSV
        tsv_path = output_dir / "predictions.tsv"
        self.get_predictions_df().to_csv(tsv_path, sep="\t", index=False)
        saved_paths["predictions_tsv"] = str(tsv_path)
        logger.info(f"Exported predictions TSV: {tsv_path}")

        # Save model
        if save_model:
            model_path = output_dir / "rf_model.pkl"
            with open(model_path, "wb") as f:
                pickle.dump(self.model, f)
            saved_paths["model"] = str(model_path)
            logger.info(f"Saved model: {model_path}")

        # Save loadings
        if save_loadings and self.loadings is not None:
            loadings_path = str(output_dir / "pca_loadings.ht")
            self.loadings.write(loadings_path, overwrite=True)
            saved_paths["loadings"] = loadings_path
            logger.info(f"Saved loadings: {loadings_path}")

        # Save pipeline stats as JSON
        import json
        stats_path = output_dir / "pipeline_stats.json"
        # Convert config to serializable format
        stats_to_save = {
            **self.pipeline_stats,
            "config": self.config.to_dict(),
            "eigenvalues": self.eigenvalues,
        }
        with open(stats_path, "w") as f:
            json.dump(stats_to_save, f, indent=2, default=str)
        saved_paths["stats"] = str(stats_path)
        logger.info(f"Saved pipeline stats: {stats_path}")

        return saved_paths


def _dataframe_to_hail_table(
    df: pd.DataFrame,
    key: str,
) -> hl.Table:
    """Convert pandas DataFrame to Hail Table via temporary file.

    This workaround avoids numpy compatibility issues with hl.Table.from_pandas()
    in certain Hail/numpy version combinations.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to convert.
    key : str
        Column to use as table key.

    Returns
    -------
    hl.Table
        Hail Table with the same data.
    """
    import tempfile
    import atexit

    # Create a temp directory that persists for the session
    tmp_dir = tempfile.mkdtemp(prefix="hail_ancestry_")
    tmp_path = str(Path(tmp_dir) / "data.tsv")

    # Write DataFrame to TSV
    df.to_csv(tmp_path, sep="\t", index=False)

    # Determine types for numeric columns
    types = {}
    for col in df.columns:
        if df[col].dtype in ("float64", "float32"):
            types[col] = hl.tfloat64
        elif df[col].dtype in ("int64", "int32"):
            types[col] = hl.tint64

    # Import table
    ht = hl.import_table(
        tmp_path,
        types=types,
        missing="",
    )
    ht = ht.key_by(key)

    # Register cleanup for when Python exits (temp files cleaned up by OS anyway)
    def cleanup():
        import shutil
        try:
            shutil.rmtree(tmp_dir, ignore_errors=True)
        except Exception:
            pass

    atexit.register(cleanup)

    return ht


def _checkpoint_or_load(
    mt: hl.MatrixTable,
    checkpoint_path: Optional[str],
    name: str,
    overwrite: bool,
) -> hl.MatrixTable:
    """Checkpoint MatrixTable or load from existing checkpoint.

    Parameters
    ----------
    mt : hl.MatrixTable
        MatrixTable to checkpoint.
    checkpoint_path : str, optional
        Base checkpoint directory.
    name : str
        Name for this checkpoint.
    overwrite : bool
        Whether to overwrite existing checkpoint.

    Returns
    -------
    hl.MatrixTable
        Checkpointed or loaded MatrixTable.
    """
    if checkpoint_path is None:
        return mt

    cp_path = str(Path(checkpoint_path) / f"{name}.mt")

    if Path(cp_path).exists() and not overwrite:
        logger.info(f"Loading from checkpoint: {cp_path}")
        return hl.read_matrix_table(cp_path)

    logger.info(f"Writing checkpoint: {cp_path}")
    return mt.checkpoint(cp_path, overwrite=overwrite)


def run_ancestry_inference(
    query_mt: hl.MatrixTable,
    reference_mt: hl.MatrixTable,
    ancestry_col: str = ANCESTRY_COL,
    config: Optional[PipelineConfig] = None,
    **kwargs: Any,
) -> AncestryInferenceResult:
    """Run complete ancestry inference pipeline.

    Performs the following steps:
    1. Merge query and reference MatrixTables
    2. Filter variants for ancestry analysis
    3. LD prune variants (optional)
    4. Compute PCA
    5. Train classifier on reference samples
    6. Predict ancestry for query samples
    7. Return complete results

    Parameters
    ----------
    query_mt : hl.MatrixTable
        Query cohort with unknown ancestry.
    reference_mt : hl.MatrixTable
        Reference panel with known ancestry labels.
    ancestry_col : str, optional
        Column name containing ancestry labels in reference_mt.
        Default: "ancestry".
    config : PipelineConfig, optional
        Pipeline configuration. If None, uses defaults.
    **kwargs
        Override individual config parameters.

    Returns
    -------
    AncestryInferenceResult
        Complete pipeline results.

    Raises
    ------
    ValueError
        If inputs are invalid or pipeline encounters errors.

    Example
    -------
    >>> result = run_ancestry_inference(
    ...     query_mt=hl.read_matrix_table("cohort.mt"),
    ...     reference_mt=hl.read_matrix_table("1kg.mt"),
    ...     ancestry_col="super_pop",
    ...     min_prob=0.80,
    ... )
    >>> predictions = result.get_predictions_df()
    >>> result.save("output/ancestry")
    """
    # Initialize config
    if config is None:
        config = PipelineConfig()

    # Apply any kwargs overrides
    if kwargs:
        config_dict = config.to_dict()
        config_dict.update(kwargs)
        config = PipelineConfig.from_dict(config_dict)

    # Validate n_pcs_classify <= n_pcs
    if config.n_pcs_classify > config.n_pcs:
        logger.warning(
            f"n_pcs_classify ({config.n_pcs_classify}) > n_pcs ({config.n_pcs}). "
            f"Reducing n_pcs_classify to {config.n_pcs}."
        )
        config = PipelineConfig.from_dict({
            **config.to_dict(),
            "n_pcs_classify": config.n_pcs,
        })

    # Initialize pipeline stats
    stats: Dict[str, Any] = {}

    logger.info("=" * 60)
    logger.info("ANCESTRY INFERENCE PIPELINE")
    logger.info("=" * 60)

    # Step 1: Merge MatrixTables
    logger.info("Step 1: Merging query and reference MatrixTables")
    stats["n_query_samples"] = query_mt.count_cols()
    stats["n_reference_samples"] = reference_mt.count_cols()
    stats["n_query_variants"] = query_mt.count_rows()
    stats["n_reference_variants"] = reference_mt.count_rows()

    merged_mt = merge_matrixtables(
        query_mt=query_mt,
        reference_mt=reference_mt,
        ancestry_col=ancestry_col,
        validate=True,
        min_shared_variants=config.min_shared_variants,
    )

    stats["n_merged_samples"] = merged_mt.count_cols()
    stats["n_shared_variants"] = merged_mt.count_rows()

    # Checkpoint merged MT
    merged_mt = _checkpoint_or_load(
        merged_mt,
        config.checkpoint_path,
        "merged",
        config.overwrite_checkpoints,
    )

    # Step 2: Filter variants
    logger.info("Step 2: Filtering variants for ancestry analysis")
    filtered_mt = filter_variants_for_ancestry(
        mt=merged_mt,
        min_af=config.min_af,
        max_af=config.max_af,
        min_call_rate=config.min_call_rate,
        apply_hwe_filter=config.apply_hwe_filter,
        hwe_p_threshold=config.hwe_p_threshold,
    )

    stats["n_variants_after_filter"] = filtered_mt.count_rows()

    # Checkpoint filtered MT
    filtered_mt = _checkpoint_or_load(
        filtered_mt,
        config.checkpoint_path,
        "filtered",
        config.overwrite_checkpoints,
    )

    # Step 3: LD pruning (optional)
    if not config.skip_ld_pruning:
        logger.info("Step 3: LD pruning variants")
        pruned_mt = ld_prune_variants(
            mt=filtered_mt,
            r2=config.ld_r2,
            bp_window_size=config.ld_window,
        )
        stats["n_variants_after_ld_prune"] = pruned_mt.count_rows()

        # Checkpoint pruned MT
        pruned_mt = _checkpoint_or_load(
            pruned_mt,
            config.checkpoint_path,
            "pruned",
            config.overwrite_checkpoints,
        )
    else:
        logger.info("Step 3: Skipping LD pruning")
        pruned_mt = filtered_mt
        stats["n_variants_after_ld_prune"] = stats["n_variants_after_filter"]

    # Step 4: Compute PCA
    logger.info("Step 4: Computing PCA")
    pca_result = compute_pca(
        mt=pruned_mt,
        n_pcs=config.n_pcs,
        compute_loadings=True,
    )

    stats["n_pcs_computed"] = pca_result.get_n_pcs()
    stats["variance_explained_pc1"] = pca_result.variance_explained()[0]
    stats["variance_explained_top5"] = sum(pca_result.variance_explained()[:5])

    # Step 5: Prepare training data
    logger.info("Step 5: Preparing training data")

    # Get scores DataFrame and merge with sample metadata
    scores_df = pca_result.get_scores_df()
    samples_df = merged_mt.cols().to_pandas()

    # Merge scores with sample metadata
    scores_with_meta = scores_df.merge(
        samples_df[["s", SOURCE_COL, KNOWN_ANCESTRY_COL]],
        on="s",
    )
    scores_with_meta = scores_with_meta.set_index("s")

    # Split into training and query
    train_scores, train_labels = get_training_samples(
        scores_with_meta,
        source_col=SOURCE_COL,
        ancestry_col=KNOWN_ANCESTRY_COL,
    )

    stats["n_training_samples"] = len(train_labels)
    stats["n_populations"] = len(train_labels.unique())
    stats["populations"] = list(train_labels.unique())

    # Step 6: Train classifier
    logger.info("Step 6: Training classifier")
    classification_result = train_classifier(
        scores_df=train_scores,
        labels=train_labels,
        n_pcs=config.n_pcs_classify,
        n_estimators=config.n_estimators,
        random_state=config.random_seed,
        validate=config.validate_model,
        n_cv_folds=config.n_cv_folds,
        min_samples_per_pop=config.min_samples_per_pop,
    )

    if classification_result.validation_metrics:
        stats["cv_accuracy"] = classification_result.get_accuracy()

    # Step 7: Predict ancestry for query samples
    logger.info("Step 7: Predicting ancestry for query samples")

    query_scores = get_query_samples(
        scores_with_meta,
        source_col=SOURCE_COL,
    )

    query_predictions = predict_ancestry(
        model=classification_result.model,
        scores_df=query_scores,
        n_pcs=config.n_pcs_classify,
        min_prob=config.min_prob,
    )

    # Also get predictions for reference samples (for validation/reporting)
    ref_predictions = predict_ancestry(
        model=classification_result.model,
        scores_df=train_scores,
        n_pcs=config.n_pcs_classify,
        min_prob=0.0,  # No threshold for reference samples
    )

    # Combine predictions
    all_predictions = pd.concat([query_predictions, ref_predictions])

    # Add source column
    all_predictions[SOURCE_COL] = scores_with_meta.loc[
        all_predictions.index, SOURCE_COL
    ]

    # Add known ancestry column
    all_predictions[KNOWN_ANCESTRY_COL] = scores_with_meta.loc[
        all_predictions.index, KNOWN_ANCESTRY_COL
    ]

    # Count predictions
    query_pred_counts = query_predictions[PREDICTED_ANCESTRY_COL].value_counts()
    stats["n_query_assigned"] = len(
        query_predictions[query_predictions[PREDICTED_ANCESTRY_COL] != "unassigned"]
    )
    stats["n_query_unassigned"] = len(
        query_predictions[query_predictions[PREDICTED_ANCESTRY_COL] == "unassigned"]
    )
    stats["query_predictions"] = query_pred_counts.to_dict()

    # Create Hail Tables from DataFrames
    logger.info("Step 8: Creating output Hail Tables")

    # Add PC columns to predictions
    n_pcs_to_include = config.n_pcs_classify
    pc_cols = [f"PC{i + 1}" for i in range(n_pcs_to_include)]
    for col in pc_cols:
        all_predictions[col] = scores_with_meta.loc[all_predictions.index, col]

    # Reset index to make sample ID a column
    all_predictions = all_predictions.reset_index()
    all_predictions = all_predictions.rename(columns={"index": "s"})

    # Convert to Hail Table via temporary file to avoid numpy compatibility issues
    predictions_ht = _dataframe_to_hail_table(all_predictions, key="s")

    # Create scores table with all PCs
    scores_for_ht = scores_with_meta.reset_index()
    scores_for_ht = scores_for_ht.rename(columns={"index": "s"})
    scores_ht = _dataframe_to_hail_table(scores_for_ht, key="s")

    # Log summary
    logger.info("=" * 60)
    logger.info("ANCESTRY INFERENCE COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Query samples: {stats['n_query_samples']}")
    logger.info(f"Assigned: {stats['n_query_assigned']}")
    logger.info(f"Unassigned: {stats['n_query_unassigned']}")
    if config.validate_model and "cv_accuracy" in stats:
        logger.info(f"Cross-validation accuracy: {stats['cv_accuracy']:.2%}")
    logger.info("Ancestry distribution:")
    for pop, count in query_pred_counts.items():
        pct = 100 * count / len(query_predictions)
        logger.info(f"  {pop}: {count} ({pct:.1f}%)")
    logger.info("=" * 60)

    return AncestryInferenceResult(
        predictions=predictions_ht,
        pc_scores=scores_ht,
        model=classification_result.model,
        loadings=pca_result.loadings,
        eigenvalues=pca_result.eigenvalues,
        config=config,
        pipeline_stats=stats,
        classification_result=classification_result,
        pca_result=pca_result,
    )
