# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

import pytest
from polartools import load_data
from os.path import join
from databroker import catalog
from polartools.manage_database import from_databroker_inplace
from spec2nexus.spec import SpecDataFile

CATALOG_NAME = "data_2"


@pytest.fixture
def db():
    if CATALOG_NAME not in list(catalog):
        path = join("polartools", "tests", "data_for_test", "databroker")
        from_databroker_inplace(path, CATALOG_NAME, catalog)
        catalog.force_reload()
    return catalog[CATALOG_NAME]


def test_load_csv():
    path = join("polartools", "tests", "data_for_test", "csv")
    table = load_data.load_table(1049, "csv", folder=path)
    assert table.shape == (4, 21)


def test_load_databroker(db):
    table = load_data.load_table(1049, db, use_db_v1=False)
    assert table.shape == (4, 19)

    table = load_data.load_table(1049, db, use_db_v1=True)
    assert table.shape == (4, 20)


def test_spec():
    path = join("polartools", "tests", "data_for_test")
    table = load_data.load_table(1, "pressure_calibration.dat", folder=path)
    assert table.shape == (51, 16)


def test_is_Bluesky_specfile():

    folder = join("polartools", "tests", "data_for_test")

    # Tests with spec loading within polartools.
    result = load_data.is_Bluesky_specfile(
        "pressure_calibration.dat", folder=folder
    )
    assert result is False

    result = load_data.is_Bluesky_specfile("bluesky_spec.dat", folder=folder)
    assert result is True

    result = load_data.is_Bluesky_specfile("absorption.dat", folder=folder)
    assert result is False

    # Tests loading the spec file here.
    spec_file = SpecDataFile(join(folder, "pressure_calibration.dat"))
    result = load_data.is_Bluesky_specfile(spec_file)
    assert result is False

    spec_file = SpecDataFile(join(folder, "bluesky_spec.dat"))
    result = load_data.is_Bluesky_specfile(spec_file)
    assert result is True

    spec_file = SpecDataFile(join(folder, "absorption.dat"))
    result = load_data.is_Bluesky_specfile(spec_file)
    assert result is False


def test_db_query(db):
    query = dict(since="2020-12-18", until="2020-12-19")
    search = load_data.db_query(db, query)
    assert len(list(search)) == 0


def test_show_meta(capsys, db):
    load_data.show_meta(
        1049, meta_keys=["motors", "plan_name", "plan_pattern_args"], db=db
    )
    captured = capsys.readouterr()
    expected = "Scan #  motors plan_name init. pos. final pos.\n"
    expected += "1049    KBIC_x rel_scan  -20        20        \n\n"
    assert captured.out == expected


def test_collect_meta(db):
    meta = load_data.collect_meta([1049], ["scan_id", "uid", "time"], db=db)
    print(meta)
    assert meta[1049]["scan_id"] == [1049, None]
    assert meta[1049]["uid"] == [
        "befa48fc-3976-4812-8eee-688223770b46",
        "f155ab58-a681-4911-9e54-0bbcb7465b22",
    ]
    assert meta[1049]["time"] == [1608235923.9519994, 1608235933.5274835]
