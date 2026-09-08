import numpy as np
import pandas as pd
import pytest
from sc_workbench import SingleCellWorkbench, io, qc, reduction


def test_orientation_and_metadata(tmp_path):
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=['a', 'b'], columns=['x', 'y', 'z'])
    path = tmp_path / 'counts.tsv'
    df.to_csv(path, sep='\t')
    assert io.load_upstream_matrix(path).shape == (3, 2)
    assert io.load_upstream_matrix(path, transpose=False).shape == (2, 3)
    wb = SingleCellWorkbench.from_df(df)
    wb.add_obs(pd.Series([20, 10], index=['b', 'a']), 'batch')
    assert wb.obs.batch.tolist() == [10, 20]
    with pytest.raises(ValueError):
        wb.add_obs(pd.Series([1], index=['missing']), 'batch')


def test_counts_and_qc_stable():
    wb = SingleCellWorkbench.from_df(pd.DataFrame([[3, 7], [5, 5]], columns=['MT-A', 'B']))
    wb.normalize()
    first = wb.to_numpy().copy()
    wb.normalize().calculate_qc()
    np.testing.assert_allclose(wb.to_numpy(), first)
    assert wb.obs.total_counts.tolist() == [10, 10]
    np.testing.assert_allclose(wb.obs.pct_counts_mito, [30, 50])
    with pytest.raises(KeyError):
        wb.to_numpy(layer='typo')
    with pytest.raises(ValueError):
        wb.cluster(method='typo')
    with pytest.raises(KeyError):
        from sklearn.preprocessing import StandardScaler
        wb.apply_sklearn(StandardScaler(), input_source='typo')


def test_pca_hvg_limit():
    wb = SingleCellWorkbench.from_numpy(np.random.default_rng(1).poisson(3, (10, 20)))
    wb.var['highly_variable'] = [True]*3 + [False]*17
    reduction.run_pca(wb.adata, n_comps=8)
    assert wb.adata.obsm['X_pca'].shape == (10, 2)


def test_bad_values_and_numpy_alignment():
    with pytest.raises(ValueError):
        io.from_dataframe(pd.DataFrame([[np.nan]]))
    with pytest.raises(ValueError):
        io.from_dataframe(pd.DataFrame([[1], [2]], index=['same', 'same']))
    a = io.from_numpy(np.ones((2, 2)), cell_names=['a', 'b'], obs=pd.DataFrame({'x':[2,1]}, index=['b','a']))
    assert a.obs.x.tolist() == [1, 2]


def test_gtf_annotation_for_mito_qc(tmp_path):
    path = tmp_path / 'genes.gtf'
    path.write_text('chrM\tsource\texon\t1\t10\t.\t+\t.\tgene_id "ENSG1.1"; gene_name "MT-A";\n')
    adata = io.from_dataframe(pd.DataFrame([[3, 7]], columns=['ENSG1.1', 'ENSG2.1']))
    io.annotate_from_gtf(adata, path)
    qc.calculate_qc_metrics(adata)
    np.testing.assert_allclose(adata.obs.pct_counts_mito, [30])
    assert adata.var.annotation_matched.tolist() == [True, False]


def test_trajectory_requires_root():
    from sc_workbench import trajectory
    with pytest.raises(ValueError, match='root'):
        trajectory.run_dpt(io.from_numpy(np.ones((5, 3))))
