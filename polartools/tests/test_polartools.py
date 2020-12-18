from polartools import load_data
from os.path import join
from os import system
from databroker import catalog


def test_load_csv():
    path = join('polartools', 'tests', 'data_for_test', 'csv')
    table = load_data.load_csv(1049, folder=path)
    assert table.shape == (4, 21)


def test_load_databroker():
    path = join('polartools', 'tests', 'data_for_test', 'databroker')
    system(f'databroker-unpack inplace {path} test_data')
    db = catalog['test_data'].get()
    table = load_data.load_databroker(1049, db)
    assert table.shape == (4, 81)
