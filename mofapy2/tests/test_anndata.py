import unittest
import pytest

import numpy as np
import pandas as pd

@pytest.mark.usefixtures("filepath_hdf5")
class TestAnnData:
    def test_multi_group(self, filepath_hdf5):
        pytest.importorskip("anndata")

        from anndata import AnnData
        from mofapy2.run.entry_point import mofa

        views_names = ["view1", "view2"]
        groups_names = ["groupA", "groupB"]

        # Set dimensions
        n_g1, n_g2 = 10, 20
        d_m = 30
        k = 5
        n = n_g1 + n_g2

        # Generate data
        np.random.seed(42)
        z1 = np.random.normal(size=(n_g1, k))
        z2 = np.random.normal(size=(n_g2, k))
        z = np.concatenate([z1, z2], axis=0)
        w = np.random.normal(size=(d_m, k))

        e1 = np.random.normal(size=(n_g1, d_m))
        e2 = np.random.normal(size=(n_g2, d_m))
        e = np.concatenate([e1, e2], axis=0)
        
        y = np.dot(z, w.T) + e
        
        # Make sample names
        samples_names = [f"sample{i}_group{g}" for g, g_size in {"A": n_g1, "B": n_g2}.items() for i in range(g_size)]
        np.random.shuffle(samples_names)
        samples_groups = [s.split("_")[1] for s in samples_names]

        adata = AnnData(X=y, obs=pd.DataFrame({"sample": samples_names, "group": samples_groups}, index=samples_names))
        
        mofa(adata, groups_label="group", outfile=filepath_hdf5, expectations=["W", "Z"])

        adata.obs['true_group'] = [s.split("_")[1] for s in adata.obs["sample"]]

        assert all(adata.obs.group.values == adata.obs.true_group.values)

        for sample, value in (("sample0_groupA", 0.209542), ("sample7_groupB", -0.276147)):
            si = np.where(adata.obs_names == sample)[0]
            assert adata.obsm["X_mofa"][si, 0] == pytest.approx(value)



if __name__ == "__main__":
    unittest.main()
