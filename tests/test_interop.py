import numpy as np
import pandas as pd
import pytest
from sklearn.cluster import DBSCAN
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

from sc_workbench import SingleCellWorkbench


@pytest.fixture
def sample_pandas_df():
    np.random.seed(42)
    genes = [f"GENE_{i}" for i in range(25)]
    cells = [f"cell_{j}" for j in range(20)]
    data = np.random.poisson(lam=10, size=(20, 25))
    return pd.DataFrame(data, index=cells, columns=genes)


def test_dataframe_interop(sample_pandas_df):
    df = sample_pandas_df
    # Ingest from Pandas DataFrame
    wb = SingleCellWorkbench.from_df(df)
    assert wb.n_cells == 20
    assert wb.n_genes == 25

    # Export back to Pandas DataFrame
    exported_df = wb.to_df()
    assert isinstance(exported_df, pd.DataFrame)
    assert exported_df.shape == (20, 25)
    assert (exported_df.index == df.index).all()
    assert (exported_df.columns == df.columns).all()


def test_numpy_interop():
    np.random.seed(42)
    X = np.random.poisson(lam=15, size=(15, 30))
    wb = SingleCellWorkbench.from_numpy(X)
    assert wb.n_cells == 15
    assert wb.n_genes == 30

    # Export back to NumPy
    X_back = wb.to_numpy()
    assert isinstance(X_back, np.ndarray)
    assert X_back.shape == (15, 30)
    assert np.allclose(X, X_back)


def test_pandas_numpy_slicing_and_metadata(sample_pandas_df):
    wb = SingleCellWorkbench.from_df(sample_pandas_df)
    wb.calculate_qc()

    # Use pandas boolean series to filter/slice workbench
    mask = wb.obs["total_counts"] > np.median(wb.obs["total_counts"])
    sub_wb = wb[mask]
    assert sub_wb.n_cells < wb.n_cells
    assert sub_wb.n_cells == mask.sum()

    # Add new metadata column using pandas Series
    batch_series = pd.Series(
        ["batch_1" if i % 2 == 0 else "batch_2" for i in range(wb.n_cells)], index=wb.obs.index
    )
    wb.add_obs(batch_series, col_name="batch_id")
    assert "batch_id" in wb.obs.columns


def test_sklearn_clustering_and_transformers(sample_pandas_df):
    wb = SingleCellWorkbench.from_df(sample_pandas_df)
    wb.normalize().run_pca(n_comps=5)

    # 1. KMeans clustering via scikit-learn
    wb.cluster(method="kmeans", n_clusters=3, key_added="kmeans_labels")
    assert "kmeans_labels" in wb.obs.columns
    assert len(wb.obs["kmeans_labels"].unique()) <= 3

    # 2. Custom Scikit-learn estimator (DBSCAN)
    wb.cluster(estimator=DBSCAN(eps=2.0, min_samples=2), key_added="dbscan_labels")
    assert "dbscan_labels" in wb.obs.columns

    # 3. Apply arbitrary Scikit-Learn transformer (StandardScaler)
    wb.apply_sklearn(StandardScaler(), input_source="X_pca", key_added="X_pca_scaled")
    assert "X_pca_scaled" in wb.adata.obsm
    assert wb.adata.obsm["X_pca_scaled"].shape == (20, 5)


def test_sklearn_supervised_classifier(sample_pandas_df):
    wb = SingleCellWorkbench.from_df(sample_pandas_df)
    wb.normalize().run_pca(n_comps=5)
    wb.cluster(method="kmeans", n_clusters=2, key_added="cell_group")

    # Train Random Forest classifier using Scikit-Learn
    clf = RandomForestClassifier(n_estimators=10, random_state=42)
    res = wb.train_classifier(clf, target_col="cell_group", use_rep="X_pca", cv=2)
    assert "mean_cv_accuracy" in res
    assert 0.0 <= res["mean_cv_accuracy"] <= 1.0
    assert "classification_report" in res
