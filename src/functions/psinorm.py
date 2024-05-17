import pandas
import numpy
import time
import tracemalloc

class PsiNorm:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        self.name = "PsiNorm"
        self.runtime = 0
        self.memory = 0
        self.transform()
        self.hvgenes = self.hvg(num_hvg)

    def transform(self):
        start = time.time()
        tracemalloc.start()
        
        n = len(self.matrix.index)
        m = (self.matrix + 1).min(0)
        log_counts = pandas.DataFrame(numpy.log(self.matrix + 1))
        log_m = pandas.Series(numpy.log(m))
        sf = log_counts.sub(log_m).sum(1) / n
        self.matrix = self.matrix.div(sf, axis = 0)

        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()
    
    def hvg(self, num_hvg : int) -> list:
        cv = self.matrix.std(0).div(self.matrix.mean(0))
        return cv.sort_values(ascending = False).index.to_list()[0:num_hvg]
        