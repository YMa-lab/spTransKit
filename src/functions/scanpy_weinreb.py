import scanpy.preprocessing._deprecated
import numpy
import pandas
import time
import tracemalloc

class ScanpyWeinreb:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        self.name = "scanpy Weinreb"
        self.runtime = 0
        self.memory = 0
        self.transform()
        self.hvgenes = self.hvg(num_hvg)

    def transform(self):
        start = time.time()
        tracemalloc.start()
        
        data = numpy.array(self.matrix)
        scanpy.preprocessing.log1p(data)
        scanpy.preprocessing._deprecated.normalize_per_cell_weinreb16_deprecated(data)
        self.matrix = pandas.DataFrame(data, index = self.matrix.index, columns = self.matrix.columns)

        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()

    def hvg(self, num_hvg : int) -> list:
        cv = self.matrix.std(0).div(self.matrix.mean(0))
        return cv.sort_values(ascending = False).index.to_list()[0:num_hvg]
