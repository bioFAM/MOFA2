import pytest

@pytest.fixture(scope="module")
def filepath_hdf5(tmpdir_factory):
    yield str(tmpdir_factory.mktemp("tmp_mofa_dir").join("mofa_pytest.hdf5"))
