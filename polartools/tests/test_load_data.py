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


def test_load_hdf5_data():
    path = join("polartools", "tests", "data_for_test")
    table = load_data.load_hdf5_data(25, path)
    assert table.shape == (2, 9)

    table2 = load_data.load_table(25, source="hdf", folder=path)
    assert table2.shape == (2, 9)


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
        1, meta_keys=["motors", "plan_name", "plan_pattern_args"], db=db
    )
    captured = capsys.readouterr()
    expected = "Scan #  motors plan_name init. pos. final pos.\n"
    expected += "1049    KBIC_x rel_scan  -20        20        \n\n"
    expected = "List of the 1 most recent scans:\n"
    expected += "motors plan_name init. pos. final pos.\n"
    expected += "KBIC_x rel_scan  -20        20        \n\n"
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


def test_load_spec():
    path = join("polartools", "tests", "data_for_test")
    table = load_data.load_spec(1, "pressure_calibration.dat", folder=path)
    assert table.shape == (51, 16)


def test_load_spec_with_specfile_object():
    path = join("polartools", "tests", "data_for_test")
    spec_file = SpecDataFile(join(path, "pressure_calibration.dat"))
    table = load_data.load_spec(1, spec_file)
    assert table.shape == (51, 16)


def test_load_hdf5_master():
    path = join("polartools", "tests", "data_for_test")
    f = load_data.load_hdf5_master(25, path)
    assert "entry" in f
    f.close()


def test_hdf5_to_dataframe():
    path = join("polartools", "tests", "data_for_test")
    f = load_data.load_hdf5_master(25, path)
    location = "entry/instrument/bluesky/streams/primary"
    df = load_data.hdf5_to_dataframe(f[location])
    assert df.shape == (2, 9)
    f.close()


def test_lookup_position(db, capsys):
    load_data.lookup_position(db, 1049, search_string="")
    captured = capsys.readouterr()
    # lookup_position always prints a "Positioner" header line
    assert "Positioner" in captured.out


def test__is_tiled(db):
    import types

    assert load_data._is_tiled(db) is False
    fake_tiled = types.SimpleNamespace()  # no .v2 attribute
    assert load_data._is_tiled(fake_tiled) is True


def test_db_query_tiled():
    pytest.importorskip("tiled")
    from unittest.mock import MagicMock
    from tiled.queries import TimeRange, Key  # noqa: F401

    mock_db = MagicMock(spec=[])  # no .v2 → _is_tiled returns True
    mock_db.search.return_value = mock_db
    query = {"since": "2024-01-01", "scan_id": 1049}
    load_data.db_query(mock_db, query)
    assert mock_db.search.call_count == 2


def test_load_catalog_tiled():
    pytest.importorskip("tiled")
    from unittest.mock import patch, MagicMock

    mock_cat = MagicMock(spec=[])  # no .v2 → _is_tiled returns True
    with patch(
        "polartools.load_data.getDatabase",
        side_effect=Exception("not found"),
    ), patch("tiled.client.from_profile", return_value=mock_cat) as mock_fp:
        result = load_data.load_catalog(name="test_profile")
    mock_fp.assert_called_once_with("test_profile")
    assert result is mock_cat
    mock_cat.register_handler.assert_not_called()


def test_load_databroker_tiled():
    pytest.importorskip("tiled")
    from unittest.mock import MagicMock
    import pandas as pd

    expected_df = pd.DataFrame({"x": [1, 2, 3]})
    mock_stream = MagicMock()
    mock_stream.read.return_value = expected_df
    mock_run = MagicMock(spec=["__contains__", "__getitem__"])
    mock_run.__contains__ = MagicMock(return_value=True)
    mock_run.__getitem__ = MagicMock(return_value=mock_stream)
    mock_db = MagicMock(spec=["__getitem__"])  # no .v2 → tiled
    mock_db.__getitem__ = MagicMock(return_value=mock_run)

    result = load_data.load_databroker(1049, db=mock_db)
    pd.testing.assert_frame_equal(result, expected_df)
