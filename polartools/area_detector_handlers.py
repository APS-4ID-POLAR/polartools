# Copyright (c) 2020-2021, UChicago Argonne, LLC.
# See LICENSE file for details.

from area_detector_handlers.handlers import AreaDetectorHDF5SingleHandler
from area_detector_handlers.eiger import EigerHandler
from dask.array import from_array
from os.path import getsize
from glob import glob
from pathlib import Path
from h5py import File


class LambdaHDF5Handler(AreaDetectorHDF5SingleHandler):
    specs = {"AD_HDF5_lambda"} | AreaDetectorHDF5SingleHandler.specs
    """
    Handler for the Lambda detector HDF5 files.
    """
    def __call__(self, point_number):
        return from_array(super().__call__(point_number))


class MyEigerHandler(EigerHandler):
    """
    Modified Eiger handler -> APS seems to use a different file naming.
    """
    def __call__(self, image_num):
        '''
        This returns data contained in the file.
        Parameters
        ----------
        image_num int
            Image number as read from eiger.cam.num_images_counter
        Returns
        -------
            A dask array
        '''

        fpath = Path(
            f'{self._file_prefix}_data_'
            f'{1 + (image_num // self._images_per_file):06d}.h5'
        ).absolute()

        try:
            file = self._files[fpath]
        except KeyError:
            file = File(fpath, 'r')
            self._files[fpath] = file

        da = from_array(file['entry/data/data'])[
            image_num % self._images_per_file
        ]

        return da.reshape((1,) + da.shape)

    def get_file_list(self):
        '''
        Get the file list.
        '''
        return glob(f'{self._file_prefix}_*')

    def get_file_sizes(self):
        '''
        Get the file size.
        Returns size in bytes.
        '''
        sizes = []
        file_name = self.get_file_list()
        for file in file_name:
            sizes.append(getsize(file))

        return sizes
