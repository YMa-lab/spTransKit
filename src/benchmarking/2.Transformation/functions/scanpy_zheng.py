import scanpy
import pandas
import time
import tracemalloc

class ScanpyZheng:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        self.name = "scanpy Zheng"
        self.runtime = 0
        self.memory = 0
        self.transform()
        self.hvgenes = self.hvg(num_hvg)

    def transform(self):
        start = time.time()
        tracemalloc.start()
        
        data = scanpy.AnnData(self.matrix)
        scanpy.pp.normalize_total(data, key_added = "n_counts_all")
        scanpy.pp.log1p(data)
        scanpy.pp.scale(data)
        self.matrix = data.to_df()

        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()

    def hvg(self, num_hvg : int) -> list:
        var = self.matrix.var(0)
        return var.sort_values(ascending = False).index.to_list()[0:num_hvg]
