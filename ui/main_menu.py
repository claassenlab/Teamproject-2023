import tkinter as tk
from tkinter import *
#from tkinter.ttk import *
from PIL import ImageTk
import os
import PIL.Image
from ui.colors import Colors
from data.file_handler import FileHandler

main_bg_color = Colors.ukt_blue_dark_1
visualization_bg_color = Colors.ukt_blue_dark_1
main_button_color = Colors.ukt_gold

sidebar_width = 300
bottombar_height = 100
claasenlab_image_height = 67


class UI:

    def __init__(self, master=None):
        self.fh = FileHandler(self)
        self.file_name = self.fh.file_name

        # top UI level widget
        top_level = tk.Tk() if master is None else tk.Toplevel(master)
        top_level.title("RNA velocity")
        top_level.state("zoomed")
        top_level.configure(height=1080, width=1920, bg=main_bg_color)
        top_level.resizable(True, True)
        self.mainwindow = top_level

        # sidebar frame (left)
        self.sidebar_frame = tk.Frame(top_level)
        self.sidebar_frame.configure(
            height=1080, width=sidebar_width, bg=Colors.ukt_black)
        self.sidebar_frame.pack(side="left")
        self.sidebar_frame.pack_propagate(False)

        # analysis frame (right)
        self.analysis_frame = tk.Frame(top_level)
        self.analysis_frame.configure(
            height=1080, width=1920-sidebar_width, bg=main_bg_color)
        self.analysis_frame.pack(side="right")
        self.analysis_frame.pack_propagate(False)

        self.create_file_name_label()
        self.create_file_dnd_field()
        self.create_buttons()
        self.place_images()
        self.create_visualization_canvas()

        self.mainwindow.mainloop()

    def create_file_name_label(self):
        # text indicating the currently loaded file
        self.file_name_label = tk.Label(
            self.sidebar_frame, text=self.fh.file_name)
        self.file_name_label.pack(pady=10, side="top")

    def create_file_dnd_field(self):
        # the field you can drop files at
        self.dnd_field = tk.Frame(self.sidebar_frame)
        self.dnd_field.configure(width=200, height=200, bg=main_button_color)
        self.dnd_field.pack(pady=20, side="top")
        self.dnd_field.pack_propagate(False)
        label = tk.Label(
            self.dnd_field, text="Drop a BAM/FASTQ file here", bg=main_button_color)
        label.place(relx=0.5, rely=0.5, anchor=CENTER)

    def create_buttons(self):
        # button to browse the file explorer
        file_upload_button = tk.Button(self.sidebar_frame)
        file_upload_button.configure(
            text="Browse files", bg=main_button_color)
        file_upload_button.pack(pady=20, side="top")
        file_upload_button.configure(command=self.fh.open_file_by_explorer)

        # button to configure the analysis
        analysis_menu_button = tk.Button(self.sidebar_frame)
        analysis_menu_button.configure(
            text="Analysis menu", bg=main_button_color)
        analysis_menu_button.pack(pady=20, side="top")

        # button to run the analysis
        frame4 = tk.Frame(self.analysis_frame)
        frame4.configure(height=200, width=200, bg=main_bg_color)
        self.runbutton = tk.Button(frame4)
        self.runbutton.configure(text="▶ Run Analysis", bg=main_button_color)
        self.runbutton.configure(command=self.clickRun)
        self.runbutton.pack(padx=10, pady=10, side="right")
        frame4.pack(side="top")

        # button to save the results
        save_result_button = tk.Button(self.sidebar_frame)
        save_result_button.configure(
            text="Save results", bg=main_button_color)
        save_result_button.pack(pady=20, side="top")

    def place_images(self):
        bottom_bar = tk.Frame(self.sidebar_frame)
        bottom_bar.configure(width=sidebar_width,
                             height=bottombar_height, bg=main_bg_color)
        bottom_bar.pack(side="bottom")
        bottom_bar.pack_propagate(False)

        image_claasenlab = PIL.Image.open(
            "ui/images/claassen_lab_logo.png")
        image_claasenlab = image_claasenlab.resize(
            (claasenlab_image_height, claasenlab_image_height), PIL.Image.ANTIALIAS)
        image_claasenlab = ImageTk.PhotoImage(image_claasenlab)
        panel = Label(bottom_bar, image=image_claasenlab,
                      bg=Colors.ukt_greige)
        panel.image = image_claasenlab
        panel.pack(padx=5, side="left")

        image_uni_klinikum = PIL.Image.open(
            "ui/images/UniklinikumTübingen.png")
        # both should have the same height
        scale = image_uni_klinikum.height / claasenlab_image_height
        image_uni_klinikum = image_uni_klinikum. resize(
            (int(image_uni_klinikum.width / scale), int(image_uni_klinikum.height / scale)), PIL.Image.ANTIALIAS)
        image_uni_klinikum = ImageTk.PhotoImage(image_uni_klinikum)
        panel = Label(bottom_bar, image=image_uni_klinikum,
                      bg=Colors.ukt_greige)
        panel.image = image_uni_klinikum
        panel.pack(padx=5, side="right")

    def create_visualization_canvas(self):
        self.canvas2 = tk.Canvas(self.analysis_frame)
        self.canvas2.configure(height=200, width=500,
                               bg=visualization_bg_color)
        self.canvas2.pack(side="bottom")
        self.canvas2.update()
        coord = 10, 10, self.canvas2.winfo_width(), self.canvas2.winfo_height()
        self.arc1 = self.canvas2.create_arc(
            coord, start=0, extent=150, fill="red")
        self.arc2 = self.canvas2.create_arc(
            coord, start=150, extent=215, fill="green")
        self.canvas2.pack(fill=BOTH, expand=TRUE)
        """ canvas2.bind("<Configure>", lambda event: resize_arc(canvas2, arc1, coord)) """

    def resize_arc(canvas, arc, coord):
        canvas.coords(arc, coord)

    # def run(self):
        # self.mainwindow.mainloop()

    def clickRun(self):
        self.canvas2.delete(self.arc1)

    def updateUI(self):
        """
        This method is called after any event that should cause a UI change
        """

        # update the file name label
        self.file_name_label.config(text=self.fh.file_name)
