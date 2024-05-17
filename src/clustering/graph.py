import tkinter
import pandas

class Grapher:
    
    def __init__(self, count_matrices : list, coordinate_matrix : pandas.DataFrame, gene : str, names : list[str]):
        self.graph(count_matrices, coordinate_matrix, gene, names)


    def graph(self, count_matrices : list[pandas.DataFrame], coordinate_matrix : pandas.DataFrame, gene : str, names : list[str]):
        # Generates window
        window = tkinter.Tk()
        window.title(gene)
        window.geometry("1375x800")
        window.config(bg = "white")

        # Generates canvas
        canvas = tkinter.Canvas(background = "white", width = 1375, height = 800)
        
        # Spatially orients the cell spot positions
        len_box = self.box_length(max(coordinate_matrix.loc[:, "x"].to_list()), max(coordinate_matrix.loc[:, "y"].to_list()))
        
        # Determines the bounds for colorification
        i = 0
        for count_matrix in count_matrices:
            max_color = self.color_max(max(count_matrix.loc[:, gene].to_list()), min(count_matrix.loc[:, gene].to_list()))
            self.title(canvas, names[i], i)
            self.graph_each(canvas, count_matrix, coordinate_matrix, gene, len_box, max_color, i)
            self.key(canvas, max_color, i)
            i = i + 1
        
        canvas.pack()
        window.mainloop()

    def title(self, canvas : tkinter.Canvas, name : str, i : int):
        x = 0
        y = 0
        if i < 4:
            x = i
            y = 0
        else:
            x = i - 4
            y = 1
        
        canvas.create_text(50 + 150 + (325 * x), 40 + (375 * y), text = name, fill = "black", font = ("Helvetica", "10"))

    def graph_each(self, canvas : tkinter.Canvas, count_matrix : pandas.DataFrame, coordinate_matrix : pandas.DataFrame, gene : str, len_box : float, max_color : float, i : int):
        x = 0
        y = 0
        if i < 4:
            x = i
            y = 0
        else:
            x = i - 4
            y = 1
        
        for cell in count_matrix.index:
            x0 = 50 + (325 * x) + ((coordinate_matrix.at[cell, "x"] - 1) * len_box)
            y0 = 50 + (375 * y) + ((coordinate_matrix.at[cell, "y"] - 1) * len_box)
            x1 = 50 + (325 * x) + ((coordinate_matrix.at[cell, "x"]) * len_box)
            y1 = 50 + (375 * y) + ((coordinate_matrix.at[cell, "y"]) * len_box)
            fill = self.color_val(count_matrix.at[cell, gene], max_color)
            canvas.create_rectangle(x0, y0, x1, y1, fill = fill)

    def box_length(self, max_x, max_y) -> float:
        len_box = 0
        if max_x < max_y:
            len_box = 300 / max_y
        else:
            len_box = 300 / max_x
        return len_box
    
    def color_max(self, max_count : float, min_count : float) -> float:
        max_color = 0
        if abs(max_count) > abs(min_count):
            max_color = abs(max_count)
        else:
            max_color = abs(min_count)
        return max_color
    
    def color_val(self, val : float, max_color : float) -> str:
        red = 0
        green = 0
        blue = 0
        
        if val > 0:
            red = 255
            green = 255 - ((val / max_color) * 255)
            blue = 255 - ((val / max_color) * 255)
        elif val == 0:
            red = 255
            green = 255
            blue = 255
        elif val < 0:
            red = 255 + ((val / max_color) * 255)
            green = 255 + ((val / max_color) * 255)
            blue = 255
            
        return "#%02x%02x%02x" % (int(red), int(green), int(blue))

    def key(self, canvas : tkinter.Canvas, max_color : float, i : int):
        x = 0
        y = 0
        if i < 4:
            x = i
            y = 0
        else:
            x = i - 4
            y = 1

        canvas.create_rectangle(50 + 143.5 + (325 * x), 50 + 320 + (375 * y), 50 + 156.5 + (325 * x), 50 + 333 + (375 * y), fill = "white")
        canvas.create_text(50 + 150 + (325 * x), 50 + 315 + (375 * y), text = "0.0", fill = "black", font = ("Helvetica", "6"))
        for interval in range(1, 11):
            # Red = upregulation
            red_key = 255
            green_key = 255 - (255 * interval / 10)
            blue_key = 255 - (255 * interval / 10)
            color_key = "#%02x%02x%02x" % (int(red_key), int(green_key), int(blue_key))
            canvas.create_rectangle((50 + 143.5 + (325 * x) - (interval * 13)), 50 + 320 + (375 * y), (50 + 143.5 + (325 * x) - ((interval - 1) * 13)), 50 + 333 + (375 * y), fill = color_key)
            canvas.create_text(50 + 150 + (325 * x) - (interval * 13), 50 + 315 + (375 * y), text = str(round(max_color / 10 * interval, 1)), fill = "black", font = ("Helvetica", "6"))
        for interval in range(1, 11):
            # Blue = downregulation
            red_key = 255 - (255 * interval / 10)
            green_key = 255 - (255 * interval / 10)
            blue_key = 255
            color_key = "#%02x%02x%02x" % (int(red_key), int(green_key), int(blue_key))
            canvas.create_rectangle((50 + 156.5 + (325 * x) + ((interval - 1) * 13)), 50 + 320 + (375 * y), (50 + 156.5 + (325 * x) + (interval * 13)), 50 + 333 + (375 * y), fill = color_key)
            canvas.create_text(50 + 150 + (325 * x) + (interval * 13), 50 + 315 + (375 * y), text = str(round(-max_color / 10 * interval, 1)), fill = "black", font = ("Helvetica", "6"))