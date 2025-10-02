import pandas

class Raw:

    def __init__(self, count_matrix : pandas.DataFrame, num_hvg : int, name : str = None):
        self.matrix = pandas.DataFrame(data = count_matrix, copy = True)
        if name:
            self.name = name
        else:
            self.name = "y"
        self.hvgenes = self.hvg(num_hvg)
        self.runtime = 0
        self.memory = 0
    
    def hvg(self, num_hvg : int) -> list:
        var = self.matrix.var(0)
        return var.sort_values(ascending = False).index.to_list()[0:num_hvg]
    