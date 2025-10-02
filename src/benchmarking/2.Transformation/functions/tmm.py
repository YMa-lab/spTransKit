import pandas
import numpy
import time
import tracemalloc
import scanpy

class TMM:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        self.name = "TMM"
        self.runtime = 0
        self.memory = 0
        self.transform()
        self.hvgenes = self.hvg(num_hvg)

    def transform(self):
        start = time.time()
        tracemalloc.start()
        
        f = self.matrix.quantile(q = 0.75, axis = 1).div(self.matrix.sum(1))
        ref_index = (abs(f - f.mean())).argmin()
        ref_row = self.matrix.iloc[ref_index, :]

        norm_factors = pandas.Series(1.0, index = self.matrix.index)

        for i in range(0, len(norm_factors.index)):
            current_row = self.matrix.iloc[i, :]
        
            numpy.seterr(divide = "ignore")
            log_ratios = numpy.log2((current_row / current_row.sum()) / (ref_row / ref_row.sum()))
            abs_exp = (numpy.log2(current_row / current_row.sum()) + numpy.log2(ref_row / ref_row.sum())) / 2
            weights = ((current_row.sum() - current_row) / current_row.sum() / current_row) + ((ref_row.sum() - ref_row) / ref_row.sum() / ref_row)

            keep = numpy.logical_and(numpy.logical_and(numpy.isfinite(log_ratios), numpy.isfinite(abs_exp)), (abs_exp > -1e10))
            log_ratios = log_ratios[keep]
            abs_exp = abs_exp[keep]
            weights = weights[keep]

            if log_ratios.max() < 1e-6:
                continue

            n = len(log_ratios)
            low_l = numpy.floor(n * 0.3) + 1
            high_l = n + 1 - low_l
            low_s = numpy.floor(n * 0.05) + 1
            high_s = n + 1 - low_s

            keep_l = numpy.logical_and(log_ratios.rank() >= low_l, log_ratios.rank() <= high_l)
            keep_s = numpy.logical_and(abs_exp.rank() >= low_s, abs_exp.rank() <= high_s)
            keep = numpy.logical_and(keep_l, keep_s)

            log_ratios = log_ratios[keep]
            weights = weights[keep]

            sf = (log_ratios / weights).sum(skipna = True) / (1 / weights).sum(skipna = True)

            if str(sf) == str(numpy.nan):
                continue
            else:
                norm_factors.iloc[i] = 2 ** sf

        norm_factors = norm_factors / numpy.exp(numpy.log(norm_factors).mean())
        self.matrix = self.matrix.div(norm_factors, axis = 0)

        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()

    def hvg(self, num_hvg : int) -> list:
        var = self.matrix.var(0)
        return var.sort_values(ascending = False).index.to_list()[0:num_hvg]
