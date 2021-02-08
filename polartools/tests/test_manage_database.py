# Copyright (c) 2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools import manage_database
from os.path import join
from glob import glob
from databroker import catalog


def test_manage_databroker_database():

    manage_database.remove_catalog('pytest_data_jfhworu1287')

    # load databroker
    path = join('polartools', 'tests', 'data_for_test', 'databroker')
    manage_database.from_databroker_inplace(path, 'pytest_data_jfhworu1287')
    catalog()  # Updates catalog list
    assert 'pytest_data_jfhworu1287' in list(catalog)

    # export databroker
    db = catalog['pytest_data_jfhworu1287']
    path = join('polartools', 'tests', 'tmp', 'test_db')
    manage_database.to_databroker(db, path, query=dict(scan_id=1049))
    files = glob(join(path, '*.*')) + glob(join(path, '*', '*'))
    assert len(files) == 3

    # export csv
    path = join('polartools', 'tests', 'tmp', 'test_csv')
    manage_database.to_csv_json(db, path, query=dict(scan_id=1049),
                                overwrite=True)
    files = glob(join(path, '*.*'))
    assert len(files) == 3

    # remove db
    manage_database.remove_catalog('pytest_data_jfhworu1287')
    catalog()  # Updates catalog list
    assert 'pytest_data_jfhworu1287' not in list(catalog)
