import pytest
from polartools import db_tools
from os.path import join
from databroker import catalog
from polartools.manage_database import from_databroker_inplace

CATALOG_NAME = 'data_2'


@pytest.fixture
def db():
    if CATALOG_NAME not in list(catalog):
        path = join('polartools', 'tests', 'data_for_test', 'databroker')
        from_databroker_inplace(path, CATALOG_NAME, catalog)
        catalog.force_reload()
    return catalog[CATALOG_NAME]


def test_db_query(db):
    query = dict(since='2020-12-18', until='2020-12-19')
    search = db_tools.db_query(db, query)
    assert len(list(search)) == 0


def test_show_meta(capsys, db):
    db_tools.show_meta(
        1049, db, meta_keys=["motors", "plan_name", "plan_pattern_args"]
    )
    captured = capsys.readouterr()
    expected = "Scan #  motors plan_name init. pos. final pos.\n"
    expected += "1049    KBIC_x rel_scan  -20        20        \n\n"
    assert captured.out == expected


def test_collect_meta(db):
    meta = db_tools.collect_meta(
        [1049], db, meta_keys=['scan_id', 'uid', 'time']
    )
    print(meta)
    assert meta[1049]['scan_id'] == [1049, None]
    assert meta[1049]['uid'] == ['befa48fc-3976-4812-8eee-688223770b46',
                                 'f155ab58-a681-4911-9e54-0bbcb7465b22']
    assert meta[1049]['time'] == [1608235923.9519994, 1608235933.5274835]
