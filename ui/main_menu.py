import tkinter as tk
from tkinter import *
from PIL import ImageTk
import PIL.Image
from ui.colors import Colors
import data.file_handler as fh
import tkinterDnD

main_bg_color = Colors.ukt_blue_dark_1
sidebar_color = Colors.ukt_black
visualization_bg_color = Colors.ukt_blue_dark_1
main_button_color = Colors.ukt_gold

sidebar_width = 300
bottombar_height = 100
claasenlab_image_height = 67


class UI:

    def __init__(self):
        # add file handler
        self.fh = fh.FileHandler(self)
        self.file_name = self.fh.file_name

        # top UI level widget
        top_level = tkinterDnD.Tk()
        top_level.title("RNA velocity")
        # "fullscreen" when starting the app
        top_level.state("zoomed")
        # window size when minimizing
        top_level.geometry("1200x700")
        # minimum window size
        top_level.minsize(width=900, height=700)
        top_level.configure(width=1920, height=1080, bg=main_bg_color)
        top_level.resizable(True, True)
        self.mainwindow = top_level

        # sidebar frame (left)
        self.sidebar_frame = tk.Frame(top_level)
        self.sidebar_frame.configure(
            width=sidebar_width, height=1080, bg=sidebar_color)
        self.sidebar_frame.pack(side="left")
        self.sidebar_frame.pack_propagate(False)

        # analysis frame (right)
        self.analysis_frame = tk.Frame(top_level)
        self.analysis_frame.configure(
            width=1920-sidebar_width, height=1080, bg=main_bg_color)
        self.analysis_frame.pack(side="right")
        self.analysis_frame.pack_propagate(False)

        # load button images here due to garbage collection
        self.dnd_field_image = PhotoImage(
            file="ui/images/buttons/dnd_field_image.png", width=208, height=208)
        self.browse_files_button_image = PhotoImage(
            file="ui/images/buttons/browse_files_button_image.png", width=208, height=58)
        self.analysis_menu_button_image = PhotoImage(
            file="ui/images/buttons/analysis_menu_button_image.png", width=208, height=58)
        self.save_results_button_image = PhotoImage(
            file="ui/images/buttons/save_results_button_image.png", width=208, height=58)
        self.run_button_image = PhotoImage(
            file="ui/images/buttons/run_button_image.png", width=208, height=58)

        # create the UI elements
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
        self.file_name_label.pack(pady=20, side="top")

    def create_file_dnd_field(self):
        # the field you can drop files into
        self.dnd_field = tk.Label(self.sidebar_frame)
        self.dnd_field.configure(
            image=self.dnd_field_image, width=200, height=200, borderwidth=0)
        self.dnd_field.pack(pady=20, side="top")
        self.dnd_field.pack_propagate(False)
        self.dnd_field.register_drop_target("*")
        self.dnd_field.bind("<<Drop>>", self.drop)

    def drop(self, event):
        # This function is called when stuff is dropped into the dnd_field
        self.fh.open_file_by_dnd(event.data)

    def create_buttons(self):
        # button to browse the file explorer
        file_upload_button = tk.Button(self.sidebar_frame)
        file_upload_button.configure(
            image=self.browse_files_button_image, width=200, height=50, borderwidth=0)
        file_upload_button.pack(pady=20, side="top")
        file_upload_button.configure(command=self.fh.open_file_by_explorer)

        # button to configure the analysis
        analysis_menu_button = tk.Button(self.sidebar_frame)
        analysis_menu_button.configure(
            image=self.analysis_menu_button_image, width=200, height=50, borderwidth=0)
        analysis_menu_button.pack(pady=20, side="top")

        # button to save the results
        save_result_button = tk.Button(self.sidebar_frame)
        save_result_button.configure(
            image=self.save_results_button_image, width=200, height=50, borderwidth=0)
        save_result_button.pack(pady=20, side="top")

        # button to run the analysis
        run_frame = tk.Frame(self.analysis_frame)
        run_frame.configure(width=200, height=200, bg=main_bg_color)
        self.runbutton = tk.Button(run_frame)
        self.runbutton.configure(
            image=self.run_button_image, width=200, height=50, borderwidth=0)
        self.runbutton.configure(command=self.clickRun)
        self.runbutton.pack(padx=10, pady=10, side="right")
        run_frame.pack(side="top")

    def place_images(self):
        bottom_bar = tk.Frame(self.sidebar_frame)
        bottom_bar.configure(width=sidebar_width,
                             height=bottombar_height, bg=main_bg_color)
        bottom_bar.pack(side="bottom")
        bottom_bar.pack_propagate(False)

        image_claasenlab = PIL.Image.open(
            "ui/images/logos/claassen_lab_logo.png")
        image_claasenlab = image_claasenlab.resize(
            (claasenlab_image_height, claasenlab_image_height), PIL.Image.ANTIALIAS)
        image_claasenlab = ImageTk.PhotoImage(image_claasenlab)
        panel = Label(bottom_bar, image=image_claasenlab,
                      bg=Colors.ukt_greige)
        panel.image = image_claasenlab
        panel.pack(padx=5, side="left")

        image_uni_klinikum = PIL.Image.open(
            "ui/images/logos/UniklinikumTübingen.png")
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
        self.canvas2.configure(width=500, height=200,
                               bg=visualization_bg_color)
        self.canvas2.pack(side="bottom")
        self.canvas2.update()
        coord = 10, 10, self.canvas2.winfo_width(), self.canvas2.winfo_height()
        self.arc1 = self.canvas2.create_arc(
            coord, start=0, extent=150, fill="red")
        self.arc2 = self.canvas2.create_arc(
            coord, start=150, extent=215, fill="green")
        self.canvas2.pack(fill=BOTH, expand=TRUE)

    def resize_arc(canvas, arc, coord):
        canvas.coords(arc, coord)

    def clickRun(self):
        self.canvas2.delete(self.arc1)

    def updateUI(self):
        """
        This method is called after any event that should cause a UI change
        """

        # update the file name label
        self.file_name_label.config(text=self.fh.file_name)