"""
Unified Fluent API Wrapper for Single-Cell Downstream Analysis.
Enables rapid method chaining and intuitive inspection for Jupyter notebooks and Python scripts.
"""

from pathlib import Path
from typing import Optional, Union, List, Dict, Tuple, Any
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp

from . import io
from . import qc
from . import preprocess
from . import reduction
from . import clustering
from . import markers
from . import annotation
from . import pathway
from . import trajectory
from . import plotting

class SingleCellWorkbench:
    """
    High-level Single-Cell Downstream Workbench orchestrator.
    Encapsulates an AnnData object and provides seamless interoperability with
    Pandas, NumPy, and Scikit-Learn.
    """
    def __init__(self, adata: ad.AnnData):
        self.adata = adata

    @classmethod
    def from_matrix(
        cls,
        file_path: Union[str, Path],
        transpose: bool = True,
        sparse: bool = True
    ) -> "SingleCellWorkbench":
        """Initialize workbench by ingesting upstream pipeline count matrix."""
        adata = io.load_upstream_matrix(file_path, transpose=transpose, sparse=sparse)
        return cls(adata)

    @classmethod
    def from_df(
        cls,
        df: pd.DataFrame,
        obs: Optional[pd.DataFrame] = None,
        var: Optional[pd.DataFrame] = None,
        is_genes_by_cells: bool = False,
        sparse: bool = True
    ) -> "SingleCellWorkbench":
        """Initialize workbench directly from a Pandas DataFrame."""
        adata = io.from_dataframe(df, obs=obs, var=var, is_genes_by_cells=is_genes_by_cells, sparse=sparse)
        return cls(adata)

    @classmethod
    def from_numpy(
        cls,
        X: Union[np.ndarray, sp.spmatrix],
        cell_names: Optional[List[str]] = None,
        gene_names: Optional[List[str]] = None,
        obs: Optional[pd.DataFrame] = None,
        var: Optional[pd.DataFrame] = None,
        sparse: bool = True
    ) -> "SingleCellWorkbench":
        """Initialize workbench directly from a NumPy array or SciPy sparse matrix."""
        adata = io.from_numpy(X, cell_names=cell_names, gene_names=gene_names, obs=obs, var=var, sparse=sparse)
        return cls(adata)

    @classmethod
    def from_10x(cls, dir_path: Union[str, Path]) -> "SingleCellWorkbench":
        """Initialize workbench by ingesting 10x Genomics folder."""
        adata = io.load_10x_directory(dir_path)
        return cls(adata)

    @classmethod
    def from_h5ad(cls, file_path: Union[str, Path]) -> "SingleCellWorkbench":
        """Initialize workbench from existing .h5ad file."""
        import scanpy as sc
        adata = sc.read_h5ad(file_path)
        return cls(adata)

    # --------------------------------------------------------------------------
    # Pythonic Slicing & Direct Attribute Access (Pandas / NumPy Friendly)
    # --------------------------------------------------------------------------
    def __getitem__(self, index) -> "SingleCellWorkbench":
        """Slice the workbench using standard Pandas/NumPy masks, indices, or cell names."""
        sub_adata = self.adata[index].copy()
        return SingleCellWorkbench(sub_adata)

    @property
    def obs(self) -> pd.DataFrame:
        """Direct access to cell metadata Pandas DataFrame."""
        return self.adata.obs

    @obs.setter
    def obs(self, value: pd.DataFrame):
        self.adata.obs = value

    @property
    def var(self) -> pd.DataFrame:
        """Direct access to gene metadata Pandas DataFrame."""
        return self.adata.var

    @var.setter
    def var(self, value: pd.DataFrame):
        self.adata.var = value

    @property
    def X(self):
        """Direct access to expression data matrix (NumPy array or SciPy sparse matrix)."""
        return self.adata.X

    @X.setter
    def X(self, value):
        self.adata.X = value

    @property
    def n_cells(self) -> int:
        return self.adata.n_obs

    @property
    def n_genes(self) -> int:
        return self.adata.n_vars

    def add_obs(
        self,
        data: Union[pd.Series, pd.DataFrame, np.ndarray, list],
        col_name: Optional[str] = None
    ) -> "SingleCellWorkbench":
        """Add new cell metadata from a Pandas Series, DataFrame, or NumPy array."""
        if isinstance(data, pd.DataFrame):
            for col in data.columns:
                self.adata.obs[col] = io._align_metadata(data, self.adata.obs_names)[col]
        elif col_name is not None:
            self.adata.obs[col_name] = io._align_metadata(data.to_frame(name=col_name), self.adata.obs_names)[col_name] if isinstance(data, pd.Series) else data
        else:
            raise ValueError("col_name must be specified when adding a Series, list, or NumPy array.")
        return self

    def to_df(self, layer: Optional[str] = None, use_raw: bool = False) -> pd.DataFrame:
        """Export expression matrix as a Pandas DataFrame (Cells x Genes)."""
        return io.to_dataframe(self.adata, layer=layer, use_raw=use_raw)

    def to_numpy(self, layer: Optional[str] = None, use_raw: bool = False) -> np.ndarray:
        """Export expression matrix as a dense 2D NumPy array."""
        return io.to_numpy(self.adata, layer=layer, use_raw=use_raw)

    def get_embedding(self, name: str = "X_umap", as_df: bool = True) -> Union[pd.DataFrame, np.ndarray]:
        """Extract dimension reduction embedding (e.g. 'X_umap', 'X_pca') as DataFrame or NumPy array."""
        if as_df:
            return io.get_embedding_df(self.adata, embedding_key=name)
        return self.adata.obsm[name]

    def __repr__(self) -> str:
        return f"<SingleCellWorkbench: {self.n_cells} cells × {self.n_genes} genes>"

    # --------------------------------------------------------------------------
    # QC & Filtering
    # --------------------------------------------------------------------------
    def calculate_qc(
        self,
        mito_prefix: Union[str, Tuple[str, ...]] = ("MT-", "mt-"),
        ribo_prefix: Union[str, Tuple[str, ...]] = ("RPS", "RPL", "rps", "rpl")
    ) -> "SingleCellWorkbench":
        """Calculate per-cell and per-gene quality control metrics."""
        self.adata = qc.calculate_qc_metrics(self.adata, mito_prefix=mito_prefix, ribo_prefix=ribo_prefix)
        return self

    def detect_doublets(self, expected_doublet_rate: float = 0.06) -> "SingleCellWorkbench":
        """Detect cellular doublets using Scrublet."""
        self.adata = qc.detect_doublets(self.adata, expected_doublet_rate=expected_doublet_rate)
        return self

    def filter_cells(
        self,
        min_genes: int = 200,
        max_genes: Optional[int] = None,
        min_counts: int = 500,
        max_counts: Optional[int] = None,
        max_pct_mito: float = 20.0,
        filter_doublets: bool = False
    ) -> "SingleCellWorkbench":
        """Filter out low-quality cells based on QC thresholds."""
        self.adata = qc.filter_cells(
            self.adata,
            min_genes=min_genes,
            max_genes=max_genes,
            min_counts=min_counts,
            max_counts=max_counts,
            max_pct_mito=max_pct_mito,
            filter_doublets=filter_doublets
        )
        return self

    def filter_genes(self, min_cells: int = 3) -> "SingleCellWorkbench":
        """Filter out non-expressed genes."""
        self.adata = qc.filter_genes(self.adata, min_cells=min_cells)
        return self

    # --------------------------------------------------------------------------
    # Preprocessing & Normalization
    # --------------------------------------------------------------------------
    def normalize(self, target_sum: float = 1e4, save_raw: bool = True) -> "SingleCellWorkbench":
        """Normalize library depth and log-transform."""
        self.adata = preprocess.normalize_and_log(self.adata, target_sum=target_sum, save_raw=save_raw)
        return self

    def select_hvg(
        self,
        n_top_genes: int = 2000,
        flavor: str = "seurat",
        batch_key: Optional[str] = None
    ) -> "SingleCellWorkbench":
        """Select Highly Variable Genes."""
        self.adata = preprocess.select_hvg(self.adata, n_top_genes=n_top_genes, flavor=flavor, batch_key=batch_key)
        return self

    def scale(self, max_value: float = 10.0) -> "SingleCellWorkbench":
        """Z-score standardize features."""
        self.adata = preprocess.scale_features(self.adata, max_value=max_value)
        return self

    # --------------------------------------------------------------------------
    # Dimensionality Reduction & Integration
    # --------------------------------------------------------------------------
    def run_pca(self, n_comps: int = 50) -> "SingleCellWorkbench":
        """Compute PCA."""
        self.adata = reduction.run_pca(self.adata, n_comps=n_comps)
        return self

    def integrate_harmony(self, batch_key: str = "batch") -> "SingleCellWorkbench":
        """Run Harmony multi-batch alignment."""
        self.adata = reduction.run_harmony(self.adata, batch_key=batch_key)
        return self

    def compute_neighbors(self, n_neighbors: int = 15, n_pcs: int = 30, use_rep: Optional[str] = None) -> "SingleCellWorkbench":
        """Build cell-cell k-NN graph."""
        self.adata = reduction.compute_neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)
        return self

    def run_umap(self, min_dist: float = 0.3, spread: float = 1.0) -> "SingleCellWorkbench":
        """Compute 2D UMAP embedding."""
        self.adata = reduction.run_umap(self.adata, min_dist=min_dist, spread=spread)
        return self

    def run_tsne(self, perplexity: float = 30.0) -> "SingleCellWorkbench":
        """Compute 2D t-SNE embedding."""
        self.adata = reduction.run_tsne(self.adata, perplexity=perplexity)
        return self

    # --------------------------------------------------------------------------
    # Clustering (Single-Cell Graph + Scikit-Learn Algorithms)
    # --------------------------------------------------------------------------
    def cluster(
        self,
        resolution: float = 0.6,
        method: str = "leiden",
        n_clusters: Optional[int] = None,
        estimator: Optional[Any] = None,
        key_added: Optional[str] = None
    ) -> "SingleCellWorkbench":
        """
        Cluster cells via single-cell graph clustering (Leiden, Louvain)
        or standard Scikit-Learn clustering (KMeans, Hierarchical, or custom estimator).
        
        Parameters
        ----------
        resolution : float, default 0.6
            Leiden/Louvain resolution parameter.
        method : str, default 'leiden'
            Algorithm choice: 'leiden', 'louvain', 'kmeans', 'hierarchical'.
        n_clusters : int, optional
            Number of clusters when using 'kmeans' or 'hierarchical'.
        estimator : sklearn clustering estimator, optional
            Any custom Scikit-Learn clustering instance (e.g. DBSCAN(), SpectralClustering()).
        key_added : str, optional
            Column name in .obs to store cluster assignments.
        """
        key = key_added or (estimator.__class__.__name__.lower() if estimator else method)
        if estimator is not None:
            self.adata = clustering.cluster_sklearn(self.adata, estimator=estimator, key_added=key)
        elif method.lower() == "kmeans":
            k = n_clusters or 5
            self.adata = clustering.cluster_kmeans(self.adata, n_clusters=k, key_added=key)
        elif method.lower() in ["hierarchical", "agglomerative"]:
            k = n_clusters or 5
            self.adata = clustering.cluster_hierarchical(self.adata, n_clusters=k, key_added=key)
        elif method.lower() == "louvain":
            self.adata = clustering.cluster_louvain(self.adata, resolution=resolution, key_added=key)
        elif method.lower() == "leiden":
            self.adata = clustering.cluster_leiden(self.adata, resolution=resolution, key_added=key)
        else:
            raise ValueError(f"Unknown clustering method: {method}")
        return self

    # --------------------------------------------------------------------------
    # Scikit-Learn Transformers, Estimators & Supervised Learning
    # --------------------------------------------------------------------------
    def apply_sklearn(
        self,
        transformer_or_estimator: Any,
        input_source: str = "X_pca",
        key_added: Optional[str] = None,
        layer: Optional[str] = None
    ) -> "SingleCellWorkbench":
        """
        Apply any Scikit-Learn Transformer, Scaler, Decomposition, or Estimator.
        
        Parameters
        ----------
        transformer_or_estimator : sklearn Transformer or Estimator
            e.g. StandardScaler(), TruncatedSVD(), KMeans(), PCA(), RobustScaler()
        input_source : str, default 'X_pca'
            Input representation: 'X', a layer name, or a key in .obsm (e.g. 'X_pca', 'X_umap').
        key_added : str, optional
            Destination key in .obs, .obsm, or .layers.
        layer : str, optional
            Specific layer to use as input if input_source is not in .obsm.
        """
        est_name = transformer_or_estimator.__class__.__name__
        
        if layer is not None:
            X = self.adata.layers[layer]
        elif input_source == "X":
            X = self.adata.X
        elif input_source in self.adata.layers:
            X = self.adata.layers[input_source]
        else:
            X = self.adata.obsm[input_source]

        # 2. Execute Scikit-learn interface
        if hasattr(transformer_or_estimator, "fit_transform"):
            res = transformer_or_estimator.fit_transform(X)
        elif hasattr(transformer_or_estimator, "fit_predict"):
            res = transformer_or_estimator.fit_predict(X)
        else:
            transformer_or_estimator.fit(X)
            res = getattr(transformer_or_estimator, "labels_", None)
            if res is None:
                raise ValueError(f"Estimator {est_name} produced no labels_ after fit()!")

        # 3. Store result intelligently
        if res.ndim == 1:
            dest = key_added or est_name.lower()
            self.adata.obs[dest] = pd.Categorical([str(r) for r in res])
        elif res.ndim == 2:
            if (layer is not None or input_source == "X" or input_source in self.adata.layers) and res.shape[1] == self.n_genes:
                dest = key_added or est_name.lower()
                self.adata.layers[dest] = res
            else:
                dest = key_added or f"X_{est_name.lower()}"
                self.adata.obsm[dest] = res
        return self

    def train_classifier(
        self,
        classifier: Any,
        target_col: str = "leiden",
        use_rep: str = "X_pca",
        cv: int = 5
    ) -> Dict[str, Any]:
        """
        Train and evaluate any Scikit-Learn Classifier (e.g. RandomForestClassifier, LogisticRegression)
        to assess cluster separability or predict cell labels.
        """
        from sklearn.model_selection import cross_val_score, cross_val_predict, StratifiedKFold
        from sklearn.metrics import classification_report

        if use_rep in self.adata.obsm:
            X = self.adata.obsm[use_rep]
        elif use_rep == "X":
            X = self.adata.X
        else:
            raise KeyError(f"Unknown representation: {use_rep}")

        y = self.adata.obs[target_col].values
        splits = list(StratifiedKFold(n_splits=cv, shuffle=True, random_state=42).split(X, y))
        scores = cross_val_score(classifier, X, y, cv=splits)
        y_pred = cross_val_predict(classifier, X, y, cv=splits)

        classifier.fit(X, y)
        report = classification_report(y, y_pred, output_dict=True)

        return {
            "mean_cv_accuracy": float(scores.mean()),
            "cv_scores": scores.tolist(),
            "classifier": classifier,
            "classification_report": report
        }

    # --------------------------------------------------------------------------
    # Marker Genes & Biomarkers
    # --------------------------------------------------------------------------
    def find_markers(
        self,
        groupby: str = "leiden",
        method: str = "wilcoxon",
        n_genes: int = 25,
        pval_cutoff: Optional[float] = None,
        logfc_cutoff: Optional[float] = None
    ) -> pd.DataFrame:
        """Find cluster-specific marker genes and return tidy DataFrame."""
        self.adata = markers.find_markers(self.adata, groupby=groupby, method=method, n_genes=n_genes)
        return markers.get_markers_df(self.adata, pval_cutoff=pval_cutoff, logfc_cutoff=logfc_cutoff)

    # --------------------------------------------------------------------------
    # Annotation & Signatures
    # --------------------------------------------------------------------------
    def score_signatures(self, signatures: Dict[str, List[str]], prefix: str = "sig_") -> "SingleCellWorkbench":
        """Score cell signatures."""
        self.adata = annotation.score_gene_signatures(self.adata, signatures, prefix=prefix)
        return self

    def annotate_cell_types(
        self,
        cluster_to_celltype: Dict[Union[str, int], str],
        cluster_key: str = "leiden",
        new_key: str = "cell_type"
    ) -> "SingleCellWorkbench":
        """Assign biological cell types to clusters."""
        self.adata = annotation.assign_cell_types(
            self.adata,
            cluster_to_celltype,
            cluster_key=cluster_key,
            new_key=new_key
        )
        return self

    # --------------------------------------------------------------------------
    # Trajectory & Pseudotime
    # --------------------------------------------------------------------------
    def infer_trajectory(self, groups: str = "leiden", root_cluster: Optional[Union[str, int]] = None) -> "SingleCellWorkbench":
        """Run PAGA and Diffusion Pseudotime."""
        self.adata = trajectory.run_paga(self.adata, groups=groups)
        self.adata = trajectory.run_dpt(self.adata, root_cluster=root_cluster, cluster_key=groups)
        return self

    # --------------------------------------------------------------------------
    # Export & Persistence
    # --------------------------------------------------------------------------
    def save(self, output_path: Union[str, Path]):
        """Save AnnData object to compressed .h5ad."""
        io.save_h5ad(self.adata, output_path)

    def export_seurat(self, output_dir: Union[str, Path]):
        """Export matrix and metadata for R / Seurat."""
        io.export_for_seurat(self.adata, output_dir)
