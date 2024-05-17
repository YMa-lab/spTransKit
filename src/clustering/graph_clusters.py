import pandas
import tkinter

class GraphCluster:

    def __init__(self, coordinate_matrix : pandas.DataFrame, labels : list[pandas.Series], names : list[str]):
        self.graph_cluster(coordinate_matrix, labels, names)

    def graph_cluster(self, coordinate_matrix : pandas.DataFrame, labels : list[pandas.Series], names : list[str]):
        # Generates window
        window = tkinter.Tk()
        window.title("Clusters")
        window.geometry("1400x800")
        window.config(bg = "white")

        frame = tkinter.Frame(window, width = 1400, height = 800)
        frame.pack(expand = True)
        canvas = tkinter.Canvas(frame, background = "white", width = 1400, height = 800, scrollregion = (0, 0, 3000, 3000))
        hbar = tkinter.Scrollbar(frame, orient = tkinter.HORIZONTAL)
        hbar.pack(side = tkinter.BOTTOM,fill = tkinter.X)
        hbar.config(command = canvas.xview)
        vbar = tkinter.Scrollbar(frame, orient = tkinter.VERTICAL)
        vbar.pack(side = tkinter.RIGHT, fill = tkinter.Y)
        vbar.config(command = canvas.yview)
        canvas.config(xscrollcommand = hbar.set, yscrollcommand = vbar.set)
        canvas.pack(expand = True)
        
        # Spatially orients the cell spot positions
        len_box = self.box_length(max(coordinate_matrix.iloc[:, 0].to_list()) + 1, max(coordinate_matrix.iloc[:, 1].to_list()) + 1)

        colors = ["red", "green", "blue", "purple", "orange", "brown", "pink"]

        i = 0
        for label in labels:
            self.title(canvas, names[i], i)
            self.graph_each(coordinate_matrix, len_box, i, canvas, colors, label)
            i = i + 1

        canvas.pack()
        window.mainloop()

    def graph_each(self, coordinate_matrix : pandas.DataFrame, len_box : float, i : int, canvas : tkinter.Canvas, colors : list, label : pandas.Series):
        x = i % 4
        y = int(i / 4)
        
        for cell in coordinate_matrix.index:
            x0 = 100 + (650 * x) + ((coordinate_matrix.loc[cell, "array_col"]) * len_box)
            y0 = 100 + (700 * y) + ((coordinate_matrix.loc[cell, "array_row"]) * len_box)
            x1 = 100 + (650 * x) + ((coordinate_matrix.loc[cell, "array_col"] + 1) * len_box)
            y1 = 100 + (700 * y) + ((coordinate_matrix.loc[cell, "array_row"] + 1) * len_box)
            canvas.create_oval(x0, y0, x1, y1, fill = colors[int(label.loc[cell])], width = 0)
        
    def box_length(self, max_x, max_y) -> float:
        len_box = 0
        if max_x < max_y:
            len_box = 600 / max_y
        else:
            len_box = 600 / max_x
        return len_box
    
    def title(self, canvas : tkinter.Canvas, name : str, i : int):
        x = i % 4
        y = int(i / 4)
        
        canvas.create_text(100 + 300 + (650 * x), 70 + (700 * y), text = name, fill = "black", font = ("Helvetica", "20"))