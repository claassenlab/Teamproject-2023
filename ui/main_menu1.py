import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog
from PIL import ImageTk, Image
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import *
import os
import PIL.Image


class UIClaasen:
    def __init__(self, master=None):
        # build ui
        toplevel1 = tk.Tk() if master is None else tk.Toplevel(master)
        toplevel1.configure(height=600, width=800, bg='#d1ccc4')
        toplevel1.resizable(True, True)
        self.frame_main = tk.Frame(toplevel1)
        self.frame_main.configure(height=600, width=800, bg='#d1ccc4')
        self.Unspliced_File = tk.Button(self.frame_main)
        self.Unspliced_File.configure(text='Upload Unspliced Data File', bg='#baa874')
        self.Unspliced_File.pack(pady=5, side="top")
        self.Unspliced_File.configure(command=self.unspliced_file_upload)
        self.Spliced_File = tk.Button(self.frame_main)
        self.Spliced_File.configure(text='Upload Spliced Date File', bg='#baa874')
        self.Spliced_File.pack(pady=5, side="top")
        self.Spliced_File.configure(command=self.spliced_file_upload)
        self.Analyse = tk.Menubutton(self.frame_main)
        self.Analyse.configure(text='Analyse', bg='#baa874')
        menu1 = tk.Menu(self.Analyse)
        self.spliced_unspliced_proportions = tk.Menu(menu1)
        menu1.add(
            tk.CASCADE,
            menu=self.spliced_unspliced_proportions,
            label='spliced/unspliced proportions')
        self.number_of_genes_which_are_usable_for_velocity_analysis = tk.Menu(
            menu1)
        menu1.add(
            tk.CASCADE,
            menu=self.number_of_genes_which_are_usable_for_velocity_analysis,
            label='number of genes which are usable for velocity analysis')
        self.how_well_the_model_assumptions_hold = tk.Menu(menu1)
        menu1.add(
            tk.CASCADE,
            menu=self.how_well_the_model_assumptions_hold,
            label='how well the model assumptions hold')
        self.statistical_evaluation = tk.Menu(menu1)
        menu1.add(
            tk.CASCADE,
            menu=self.statistical_evaluation,
            label='statistical evaluation')
        self.Analyse.configure(menu=menu1)
        self.Analyse.pack(pady=5, side="top")
        self.plotting_Menu = tk.Menubutton(self.frame_main)
        self.plotting_Menu.configure(text='Plotting Menu', bg='#baa874')
        menu2 = tk.Menu(self.plotting_Menu)
        self.phase_plot = tk.Menu(menu2)
        menu2.add(tk.CASCADE, menu=self.phase_plot, label='phase plot')
        self.scatter_plot = tk.Menu(menu2)
        menu2.add(tk.CASCADE, menu=self.scatter_plot, label='scatter_plot')
        self.stream_plot = tk.Menu(menu2)
        menu2.add(tk.CASCADE, menu=self.stream_plot, label='stream_plot')
        self.table = tk.Menu(menu2)
        menu2.add(tk.CASCADE, menu=self.table, label='table')
        self.plotting_Menu.configure(menu=menu2)
        self.plotting_Menu.pack(pady=5, side="top")
        self.safe_as = tk.Menubutton(self.frame_main)
        self.safe_as.configure(text='safe as:', bg='#baa874')
        menu3 = tk.Menu(self.safe_as)
        self.safe_as_pdf = tk.Menu(menu3)
        menu3.add(tk.CASCADE, menu=self.safe_as_pdf, label='pdf')
        self.safe_as_png = tk.Menu(menu3)
        menu3.add(tk.CASCADE, menu=self.safe_as_png, label='png')
        self.safe_as.configure(menu=menu3)
        self.safe_as.pack(pady=5, side="top")
        self.frame_main.pack(side="left")
        frame4 = tk.Frame(toplevel1)
        frame4.configure(height=200, width=200, bg='#d1ccc4')
        self.runbutton = tk.Button(frame4)
        self.runbutton.configure(text='▶ Run ', bg='#baa874')
        self.runbutton.configure(command=self.clickRun)
        self.runbutton.pack(pady=10, side="right")
        frame4.pack(side="top")

    
        self.canvas2 = tk.Canvas(toplevel1)
        self.canvas2.configure(height=200, width=500)
        self.canvas2.pack(padx=20,pady=20, side="bottom")
        self.canvas2.update()
        coord = 10, 10, self.canvas2.winfo_width(), self.canvas2.winfo_height()
        self.arc1 = self.canvas2.create_arc(coord, start=0, extent=150, fill="red")
        self.arc2 = self.canvas2.create_arc(coord, start=150, extent=215, fill="green")
        self.canvas2.pack(fill=BOTH,expand=TRUE)
        """ canvas2.bind('<Configure>', lambda event: resize_arc(canvas2, arc1, coord)) """

        image_Claasen = PIL.Image.open("/home/ulyana/Schreibtisch/UI_Teamprojekt/Teamproject-2023/ui/images/claassen_lab_logo.png")
        image_Claasen = image_Claasen.resize((200, 200), PIL.Image.ANTIALIAS)
        image_Claasen = ImageTk.PhotoImage(image_Claasen)
        panel = Label(self.frame_main, image=image_Claasen, bg='#d1ccc4')
        panel.image = image_Claasen
        panel.pack(pady=10, side="top")

        imageUniKlinikum = PIL.Image.open('/home/ulyana/Schreibtisch/UI_Teamprojekt/Teamproject-2023/ui/images/UniklinikumTübingen.png')
        imageUniKlinikum = imageUniKlinikum.resize((137*2, 46*2), PIL.Image.ANTIALIAS)
        imageUniKlinikum = ImageTk.PhotoImage(imageUniKlinikum)
        panel = Label(self.frame_main, image=imageUniKlinikum, bg='#d1ccc4')
        panel.image = imageUniKlinikum
        panel.pack(side="left", anchor="sw")


        # Main widget
        self.mainwindow = toplevel1
    def resize_arc(canvas, arc, coord):
        canvas.coords(arc, coord)

    def run(self):
        self.mainwindow.mainloop()

    def unspliced_file_upload(self):
        self.filename = filedialog.askopenfilename(initialdir="", title="Select A File", filetypes=(("BAM files", "*.bam"),("FASTQ Files", "*.fastq")))

    def spliced_file_upload(self):
        self.filename = filedialog.askopenfilename(initialdir="", title="Select A File", filetypes=(("BAM files", "*.bam"),("FASTQ Files", "*.fastq")))

    def clickRun(self):
        self.canvas2.delete(self.arc1)


if __name__ == "__main__":
    app = UIClaasen()
    app.run()
