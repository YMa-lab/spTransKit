import pandas
import numpy
import time
import tracemalloc

class Acosh:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        self.name = "acosh(2αy/s + 1)"
        self.runtime = 0
        self.memory = 0
        self.transform()
        self.hvgenes = self.hvg(num_hvg)
    
    def size_factor(self) -> pandas.Series:
        return self.matrix.sum(1) / (self.matrix.sum().sum() / len(self.matrix.index))

    def overdispersion(self) -> pandas.Series:
        return pandas.Series(0.05, index = self.matrix.sum(0).index)

    def transform(self):
        start = time.time()
        tracemalloc.start()
        
        self.matrix = self.matrix.div(self.size_factor(), axis = 0)
        self.matrix = (self.matrix.mul(self.overdispersion(), axis = 1) * 2) + 1
        self.matrix = pandas.DataFrame((1 / numpy.sqrt(0.05)) * numpy.arccosh(self.matrix))
        
        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()
    
    def hvg(self, num_hvg : int) -> list:
        var = self.matrix.var(0)
        return var.sort_values(ascending = False).index.to_list()[0:num_hvg]
