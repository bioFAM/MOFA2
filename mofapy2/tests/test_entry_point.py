import unittest
import pytest

import numpy as np
import pandas as pd

from mofapy2.run.entry_point import entry_point


class TestBuildModel:
    def test_build_basic(self):
        ent = entry_point()
        ent.set_data_options(scale_groups=False, scale_views=False)
        views_names = ["view1", "view2"]
        groups_names = ["groupA", "groupB"]

        # Set dimensions
        n_g1, n_g2 = 10, 20
        d_m1, d_m2 = 30, 40
        np.random.seed(42)
        ent.set_data_matrix(
            [
                [np.random.random((n_g1, d_m1)), np.random.random((n_g2, d_m1))],
                [np.random.random((n_g1, d_m2)), np.random.random((n_g2, d_m2))],
            ]
        )
        
        ent.set_model_options()
        ent.set_train_options()
        ent.build()


if __name__ == "__main__":
    unittest.main()
