# Copyright (c) 2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools import manage_database
from os.path import join
from glob import glob
from databroker import catalog

CATALOG_NAME = 'data_1'


def test_manage_databroker_database(tmpdir):

    # load databroker
    path = join('polartools', 'tests', 'data_for_test', 'databroker')
    manage_database.from_databroker_inplace(
        path, CATALOG_NAME, catalog
        )
    catalog.force_reload()
    assert CATALOG_NAME in list(catalog)

    # export databroker
    db = catalog[CATALOG_NAME]
    path = str(tmpdir.mkdir('test_db'))
    manage_database.to_databroker(db, path, query=dict(scan_id=1049))
    files = glob(join(path, '*.*')) + glob(join(path, '*', '*'))
    assert len(files) == 3

    # export csv
    path = str(tmpdir.mkdir('test_csv'))
    manage_database.to_csv_json(db, path, query=dict(scan_id=1049),
                                overwrite=True)
    files = glob(join(path, '*.*'))
    assert len(files) == 3

    # remove db
    manage_database.remove_catalog(CATALOG_NAME, catalog)
    catalog.force_reload()
    assert CATALOG_NAME not in list(catalog)
