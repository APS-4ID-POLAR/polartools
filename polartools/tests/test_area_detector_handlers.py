# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools import area_detector_handlers
from os.path import join
from databroker import catalog
from polartools.manage_database import from_databroker_inplace


def test_lambda250k():
    path = join("polartools", "tests", "data_for_test", "lambda250k")
    if "lambda250k" not in list(catalog):
        from_databroker_inplace(path, "lambda250k", catalog)
        catalog.force_reload()
    cat = catalog["lambda250k"]
    # TODO: This is a temporary fix because of how I saved the data before!
    tmp = area_detector_handlers.LambdaHDF5Handler
    tmp.specs = {"AD_HDF5_lambda"}
    cat.register_handler("AD_HDF5_lambda", tmp, overwrite=True)
    assert (
        cat[276].primary.read()['lambda250k_image'].shape == (50, 1, 516, 516)
    )
