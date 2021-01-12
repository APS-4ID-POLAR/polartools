# Copyright (c) 2020, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools import load_data
from os.path import join
from databroker import catalog


def test_load_csv():
    path = join('polartools', 'tests', 'data_for_test', 'csv')
    table = load_data.load_csv(1049, folder=path)
    assert table.shape == (4, 21)


def test_load_databroker():
    db = catalog['test_data'].get()
    table = load_data.load_databroker(1049, db)

    # TODO: Check this. I got different results depending on versions.
    assert table.shape == (4, 19)
