# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

import pytest
from polartools import load_data
from os.path import join
from databroker import catalog
from polartools.manage_database import from_databroker_inplace
from spec2nexus.spec import SpecDataFile

CATALOG_NAME = 'data_2'


@pytest.fixture
def db():
    if CATALOG_NAME not in list(catalog):
        path = join('polartools', 'tests', 'data_for_test', 'databroker')
        from_databroker_inplace(path, CATALOG_NAME, catalog)
        catalog.force_reload()
    return catalog[CATALOG_NAME]


def test_load_csv():
    path = join('polartools', 'tests', 'data_for_test', 'csv')
    table = load_data.load_table(1049, 'csv', folder=path)
    assert table.shape == (4, 21)


def test_load_databroker(db):
    table = load_data.load_table(1049, db)
    assert table.shape == (4, 19)

    table = load_data.load_table(1049, db, use_db_v1=True)
    assert table.shape == (4, 20)


def test_spec():
    path = join('polartools', 'tests', 'data_for_test')
    table = load_data.load_table(1, 'pressure_calibration.dat', folder=path)
    assert table.shape == (51, 16)


def test_db_query(db):
    query = dict(since='2020-12-18', until='2020-12-19')
    search = load_data.db_query(db, query)
    assert len(list(search)) == 0


def test_is_Bluesky_specfile():

    folder = join('polartools', 'tests', 'data_for_test')

    # Tests with spec loading within polartools.
    result = load_data.is_Bluesky_specfile('pressure_calibration.dat',
                                           folder=folder)
    assert result is False

    result = load_data.is_Bluesky_specfile('bluesky_spec.dat',
                                           folder=folder)
    assert result is True

    result = load_data.is_Bluesky_specfile('absorption.dat',
                                           folder=folder)
    assert result is False

    # Tests loading the spec file here.
    spec_file = SpecDataFile(join(folder, 'pressure_calibration.dat'))
    result = load_data.is_Bluesky_specfile(spec_file)
    assert result is False

    spec_file = SpecDataFile(join(folder, 'bluesky_spec.dat'))
    result = load_data.is_Bluesky_specfile(spec_file)
    assert result is True

    spec_file = SpecDataFile(join(folder, 'absorption.dat'))
    result = load_data.is_Bluesky_specfile(spec_file)
    assert result is False
