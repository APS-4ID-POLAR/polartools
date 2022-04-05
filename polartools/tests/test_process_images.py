
# Copyright (c) 2022, UChicago Argonne, LLC.
# See LICENSE file for details.

from polartools.load_data import load_catalog
from polartools.manage_database import from_databroker_inplace
import polartools.process_images as process_images
import pytest
from databroker import catalog
from os.path import join
from numpy import nanmax, allclose, mean

CATALOG_NAME = "lambda250k"


@pytest.fixture
def cat():
    if CATALOG_NAME not in list(catalog):
        path = join("polartools", "tests", "data_for_test", "lambda250k")
        from_databroker_inplace(path, CATALOG_NAME, catalog)
        catalog.force_reload()
    return load_catalog(CATALOG_NAME)


@pytest.fixture
def im(cat):
    return process_images.load_images(
        [276], cat, "lambda250k_image", cleanup=dict(threshold=(100,)),
    )


@pytest.fixture
def ims(cat):
    return process_images.load_images(
        [276],
        cat,
        "lambda250k_image",
        cleanup=dict(threshold=(100,)),
        positioner="rxes_motors_arot"
    )


def test_load_image(cat):

    # with threshold cleanup and normalization, no positioner
    image = process_images.load_images(
        [276],
        cat,
        "lambda250k_image",
        cleanup=dict(threshold=(100,),),
        normalize='Ion Ch 4'
    )

    assert image.shape == (516, 516)
    assert allclose(nanmax(image).compute(), 3.845249593597569e-05)

    # Only positioner
    images, positioner = process_images.load_images(
        [276],
        cat,
        "lambda250k_image",
        positioner="rxes_motors_arot"
    )

    assert images.shape == (50, 516, 516)
    assert allclose(nanmax(images).compute(), 4095.0)
    assert positioner.shape == (50,)
    assert allclose(nanmax(positioner), 52.550000000000125)

    def custom_clean(images, value):
        new_images = images.copy()
        new_images[new_images > value] = 100
        return new_images

    # Test custom function
    image = process_images.load_images(
        [276],
        cat,
        'lambda250k_image',
        cleanup=dict(function=(custom_clean, (100,)))
    )

    assert allclose(nanmax(image), 100)


def test_get_curvature(im):
    curvature = process_images.get_curvature(im, binx=64)

    assert allclose(
        curvature, [-8.28334263e-05, 1.33928571e-02, 1.64000000e+02]
    )


def test_get_spectrum(im):
    spectrum = process_images.get_spectrum(
        im, [-8.28334263e-05, 1.33928571e-02, 1.64000000e+02], biny=2
    )

    assert spectrum.shape == (264, 2)
    assert allclose(spectrum.mean(), 132.42162878787877)


def test_get_spectra(ims):
    spectra = process_images.get_spectra(
        ims[0], [-8.28334263e-05, 1.33928571e-02, 1.64000000e+02], biny=2
    )

    assert spectra.shape == (50, 264, 2)
    assert allclose(mean(spectra), 132.42162878787877)


def test_process_rxes(cat):

    # No positioner
    spectrum = process_images.process_rxes(
        [276],
        cat,
        "lambda250k_image",
        [-8.28334263e-05, 1.33928571e-02, 1.64000000e+02],
        positioner=None,
    )
    assert spectrum.shape == (529, 2)

    # With positioner
    spectra, positioner = process_images.process_rxes(
        [276],
        cat,
        "lambda250k_image",
        [-8.28334263e-05, 1.33928571e-02, 1.64000000e+02],
        positioner="rxes_motors_arot",
    )
    assert spectra.shape == (50, 529, 2)
    assert positioner.shape == (50,)
