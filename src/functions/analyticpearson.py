import pandas
import time
import tracemalloc

class PearsonResidual:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        self.name = "Analytic Pearson"
        self.runtime = 0
        self.memory = 0
        self.transform()
        self.hvgenes = self.hvg(num_hvg)

    def n(self) -> pandas.Series:
        return self.matrix.sum(1)
    
    def p(self) -> pandas.Series:
        return self.matrix.sum(0) / self.n().sum()
    
    def overdispersion(self) -> pandas.Series:
        return pandas.Series(0.05, index = self.matrix.sum(0).index)
    
    def transform(self):
        start = time.time()
        tracemalloc.start()

        nb_matrix = pandas.DataFrame(0, index = self.matrix.index, columns = self.matrix.columns)
        nb_matrix = nb_matrix.add(self.n(), axis = 0)
        nb_matrix = nb_matrix.multiply(self.p(), axis = 1)
        var_matrix = pandas.DataFrame(data = nb_matrix, copy = True)
        var_matrix = var_matrix.pow(2)
        var_matrix = var_matrix.multiply(self.overdispersion(), axis = 1)
        var_matrix = var_matrix.add(nb_matrix)
        var_matrix = var_matrix.pow(0.5)
        self.matrix = self.matrix.sub(nb_matrix)
        self.matrix = self.matrix.div(var_matrix)
        
        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()
    
    def hvg(self, num_hvg : int) -> list:
        cv = self.matrix.std(0).div(self.matrix.mean(0))
        return cv.sort_values(ascending = False).index.to_list()[0:num_hvg]
