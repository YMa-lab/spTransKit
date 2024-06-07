import pydeseq2.preprocessing
import pandas
import numpy
import time
import tracemalloc

class DESeq2:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        self.name = "DESeq2"
        self.runtime = 0
        self.memory = 0
        self.transform()
        self.hvgenes = self.hvg(num_hvg)

    def transform(self):
        start = time.time()
        tracemalloc.start()

        # log_counts = pandas.DataFrame(numpy.log(self.matrix), index = self.matrix.index, columns = self.matrix.columns)
        # log_means = log_counts.mean(0)
        # filtered_genes = log_means.loc[log_means != float("-inf")].index.to_list()
        # log_ratios = log_counts.loc[:, filtered_genes].sub(log_means.loc[filtered_genes], axis = 1)

        log_counts = pandas.DataFrame(numpy.log(self.matrix + 1), index = self.matrix.index, columns = self.matrix.columns)
        log_means = log_counts.mean(0)
        filtered_genes = log_means.loc[log_means > 1.0].index.to_list()
        log_ratios = log_counts.loc[:, filtered_genes].sub(log_means.loc[filtered_genes], axis = 1)

        log_medians = log_ratios.median(axis = 1)
        size_factors = pandas.Series(numpy.exp(log_medians), index = log_medians.index)
        self.matrix = self.matrix.div(size_factors, axis = 0)

        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()

    def hvg(self, num_hvg : int) -> list:
        cv = self.matrix.std(0).div(self.matrix.mean(0))
        return cv.sort_values(ascending = False).index.to_list()[0:num_hvg]
