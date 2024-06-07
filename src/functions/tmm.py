import rnanorm
import conorm
import pandas
import numpy
import time
import tracemalloc

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
        
        ref_row = self.matrix.median(axis = 0).argmin()
        ref_counts = self.matrix.iloc[ref_row, :]
        ref_counts = numpy.array(ref_counts / ref_counts.sum())
        ref_counts[ref_counts == 0.0] = numpy.nan
        norm_factors = [1] * len(self.matrix.index.to_list())
        
        for i in range(0, self.matrix.shape[0]):
            current_counts = self.matrix.iloc[i, :]
            current_counts = numpy.array(current_counts / current_counts.sum())
            current_counts[current_counts == 0] = numpy.nan
            
            log_ratios = numpy.log2((current_counts) / (ref_counts))
            abs_log_ratios = numpy.abs(log_ratios)
            weights = (1 / (current_counts)) + (1 / (ref_counts))
            
            finite = numpy.isfinite(log_ratios)
            log_ratios = log_ratios[finite]
            abs_log_ratios = abs_log_ratios[finite]
            weights = weights[finite]
            
            trim_fraction = 0.2
            lower_trim = int(numpy.floor(trim_fraction * len(log_ratios)))
            upper_trim = int(numpy.ceil((1 - trim_fraction) * len(log_ratios)))
            sorted_indices = numpy.argsort(abs_log_ratios)
            trimmed_indices = sorted_indices[lower_trim:upper_trim]
            trimmed_log_ratios = log_ratios[trimmed_indices]
            trimmed_weights = weights[trimmed_indices]
            weighted_mean_log_ratio = numpy.sum(trimmed_log_ratios * trimmed_weights) / numpy.sum(trimmed_weights)
            norm_factors[i] = 2 ** weighted_mean_log_ratio
        
        norm_factors = pandas.Series(norm_factors, index = self.matrix.index)
        self.matrix = self.matrix.div(norm_factors, axis = 0)

        #self.matrix = conorm.tmm(self.matrix.T).T

        #tmm = rnanorm.TMM()
        #self.matrix = pandas.DataFrame(tmm.fit_transform(X = self.matrix), index = self.matrix.index, columns = self.matrix.columns)

        self.runtime = time.time() - start
        self.memory = tracemalloc.get_tracemalloc_memory() / 1000000
        tracemalloc.stop()

    def hvg(self, num_hvg : int) -> list:
        cv = self.matrix.std(0).div(self.matrix.mean(0))
        return cv.sort_values(ascending = False).index.to_list()[0:num_hvg]
