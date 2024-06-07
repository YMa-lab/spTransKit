import pandas
import numpy
import time
import tracemalloc

class CPMShiftedLog:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        self.name = "log(CPM + 1)"
        self.runtime = 0
        self.memory = 0
        self.transform()
        self.hvgenes = self.hvg(num_hvg)

    def cpm(self) -> pandas.Series:
        return self.matrix.sum(1) / 1000000
    
    def transform(self):
        start = time.time()
        tracemalloc.start()

        self.matrix = self.matrix.div(self.cpm(), axis = 0) + 1
        self.matrix = pandas.DataFrame(numpy.log10(self.matrix))
        
        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()

    def hvg(self, num_hvg : int) -> list:
        cv = self.matrix.std(0).div(self.matrix.mean(0))
        return cv.sort_values(ascending = False).index.to_list()[0:num_hvg]
