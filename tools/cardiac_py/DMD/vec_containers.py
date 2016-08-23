import numpy as np
import modred as MR
from copy import deepcopy
import pio.iter_read


class CustomVector(MR.Vector):
    """
    This class is an interface which allows modred to operate on the dataset
    """
    def __init__(self, data_array):
        self.data_array = data_array

    def __add__(self, other):
        """Return a new object that is the sum of self and other"""
        sum_vec = deepcopy(self)
        sum_vec.data_array = self.data_array + other.data_array
        return sum_vec

    def __mul__(self, scalar):
        """Return a new object that is ``self * scalar`` """
        mult_vec = deepcopy(self)
        mult_vec.data_array = mult_vec.data_array * scalar
        return mult_vec

    def inner_product(self, other):
        return np.dot(self.data_array, other.data_array)


class CustomVecHandle(MR.VecHandle):
    '''
    This is a 'out of memory' interface which allows modred to load data
    progressively. It also allows mod_red to write data. In order to produce a
    dataset compliant with read_results
    '''
    def __init__(self, vec_path, write_stem='DMD_mode', base_handle=None,
                 scale=None):
        MR.VecHandle.__init__(self, base_handle, scale)
        self.vec_path = vec_path

        ''' this should be changed'''
        self.write_directory = './'
        self.write_stem = write_stem

    def _get(self):
        # lets presume this
        reader = pio.iter_read.Reader(self.vec_path)
        n_points = int(reader.header_vars["nrecords"])
        data = np.zeros((n_points))
        ii = 0
        for line in reader:
            data[ii] = float(line.rstrip().split(" ")[-1])
            ii = ii + 1
        return CustomVector(data)

    def _put(self, vec):
        # Here we will use the template to write the data
        # get vector data array
        writer = pio.append_complex.Piofy(self.vec_path, self.write_stem,
                                           vec.data_array)
        writer.Write()


def inner_product(v1, v2):
    '''
    Inner product wrapper. Here we create a wrapper of an inner product.
    Not sure whether
    '''
    return v1.inner_product(v2)
