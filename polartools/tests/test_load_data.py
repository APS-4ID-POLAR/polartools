# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools import load_data
from os.path import join
from databroker import catalog
from spec2nexus.spec import SpecDataFile


def test_load_csv():
    path = join('polartools', 'tests', 'data_for_test', 'csv')
    table = load_data.load_table(1049, 'csv', folder=path)
    assert table.shape == (4, 21)


def test_load_databroker():
    db = catalog['test_data']
    table = load_data.load_table(1049, db)
    assert table.shape == (4, 19)


def test_spec():
    path = join('polartools', 'tests', 'data_for_test')
    table = load_data.load_table(1, 'pressure_calibration.dat', folder=path)
    assert table.shape == (51, 16)


def test_db_query():
    query = dict(since='2020-12-18', until='2020-12-19')
    db = catalog['test_data']
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
