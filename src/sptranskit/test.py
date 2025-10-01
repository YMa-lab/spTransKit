import transformations as sp
import filter as f
import helpers as h
import scanpy as sc
import pandas as pd
import numpy as np

if __name__ == "__main__":
    x, coord = h.get_unfiltered_dlpfc_data("151673")
    x, coord = f.filter_counts(x, coord)
    x = sp.analytic_pearson_clip(x)
