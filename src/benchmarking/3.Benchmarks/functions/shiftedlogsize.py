import pandas
import numpy
import time
import tracemalloc

class ShiftedLogSize:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        self.name = "log(y/s + 1)/u"
        self.runtime = 0
        self.memory = 0
        self.transform()
        self.hvgenes = self.hvg(num_hvg)
    
    def size_factor(self) -> pandas.Series:
        return self.matrix.sum(1) / (self.matrix.sum().sum() / len(self.matrix.index))
    
    def transform(self):
        start = time.time()
        tracemalloc.start()
        
        self.matrix = self.matrix.div(self.size_factor(), axis = 0) + 1
        self.matrix = pandas.DataFrame(numpy.log(self.matrix))
        u = self.size_factor()
        self.matrix = self.matrix.div(u, axis = 0)

        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()
    
    def hvg(self, num_hvg : int) -> list:
        var = self.matrix.var(0)
        return var.sort_values(ascending = False).index.to_list()[0:num_hvg]
