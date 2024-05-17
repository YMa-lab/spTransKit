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
        scanpy.pp.recipe_zheng17(data, n_top_genes = len(self.matrix.columns))
        self.matrix = data.to_df()

        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()

    def hvg(self, num_hvg : int) -> list:
        cv = self.matrix.std(0).div(self.matrix.mean(0))
        return cv.sort_values(ascending = False).index.to_list()[0:num_hvg]
