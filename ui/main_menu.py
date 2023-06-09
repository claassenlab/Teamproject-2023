from tkinter import *
import tkinter.messagebox as messagebox
from PIL import ImageTk
import PIL.Image
from ui.colors import Colors
import data.file_handler as fh
from tkinterdnd2 import *
import ui.tooltip as tooltip
from data.analysis import Analysis

main_bg_color = Colors.ukt_black
sidebar_color = Colors.ukt_black
visualization_bg_color = Colors.ukt_blue_dark_1
main_button_color = Colors.ukt_gold
default_label_color = Colors.ukt_white
default_font_color = Colors.ukt_black
analysis_menu_bg_color = Colors.ukt_black
analysis_menu_button_color = Colors.ukt_white
analysis_menu_text_color = Colors.ukt_black

sidebar_width = 300
bottombar_height = 100
claassenlab_image_height = 78
min_height = 700
min_width_window = 900
min_width = 1200
base_height = 1080
base_width = 1920
button_height = 58
click_field_height = 50
button_width = 208
click_field_width = 200
min_window_dimensions = "1200x700"
data_overview_font_size = 15
do_pre_window_text = "Enter the name of the text file to be created.\nLeave blank to not create one."
do_pre_window_w = 400
do_pre_window_h = 120
do_pre_window_font_size = 15
standard_font = "Calibri"
analysis_menu_x_offset = 140
analysis_menu_width = 1000
analysis_menu_height = 600
analysis_menu_font = ("Calibri", 15)
analysis_menu_section_height = 40

# error messages
no_file_message_type = "warning"
no_file_title = "No file loaded!"
no_file_message = "A file must be loaded for this! Please load a file and try again!"
no_analysis_selected_type = "warning"
no_analysis_selected_title = "No analysis selected!"
no_analysis_selected_message = "You have not selected an analysis option!\nTo do this, open the analysis menu!"
data_overview_error_type = "error"
data_overview_error_title = "Data overview error!"
data_overview_error_message = "An error occurred while performing the data overview!\nPlease check if the file and the name you entered are valid!"
umap_error_type = "error"
umap_error_title = "UMAP error!"
umap_error_message = "The UMAP projection could not be created!\nPlease check if your file and the\nselected UMAP options are valid!"


class UI:

    def __init__(self):
        # add file handler
        self.fh = fh.FileHandler(self)
        self.file_name = self.fh.file_name

        # add analysis object
        self.analysis = Analysis()

        # top UI level widget
        top_level = TkinterDnD.Tk()
        top_level.title("RNA velocity")
        # "fullscreen" when starting the app
        top_level.state("zoomed")
        # window size when minimizing
        top_level.geometry(min_window_dimensions)
        # minimum window size
        top_level.minsize(width=min_width_window, height=min_height)
        top_level.configure(
            width=base_width, height=base_height, bg=main_bg_color)
        top_level.resizable(True, True)
        self.mainwindow = top_level

        # sidebar frame (left)
        self.sidebar_frame = Frame(top_level)
        self.sidebar_frame.configure(
            width=sidebar_width, height=base_height, bg=sidebar_color)
        self.sidebar_frame.pack(side="left")
        self.sidebar_frame.pack_propagate(False)

        # analysis frame (right)
        self.analysis_frame = Frame(top_level)
        self.analysis_frame.configure(
            width=base_width-sidebar_width, height=base_height, bg=main_bg_color)
        self.analysis_frame.pack(side="right")
        self.analysis_frame.pack_propagate(False)

        # load button images here due to garbage collection
        self.dnd_field_image = PhotoImage(
            file="ui/images/buttons/dnd_field_image.png", width=button_width, height=button_width)
        self.browse_files_button_image = PhotoImage(
            file="ui/images/buttons/browse_files_button_image.png", width=button_width, height=button_height)
        self.analysis_menu_button_image = PhotoImage(
            file="ui/images/buttons/analysis_menu_button_image.png", width=button_width, height=button_height)
        self.save_results_button_image = PhotoImage(
            file="ui/images/buttons/save_results_button_image.png", width=button_width, height=button_height)
        self.overview_button_image = PhotoImage(
            file="ui/images/buttons/overview_button_image.png", width=button_width, height=button_height)
        self.run_button_image = PhotoImage(
            file="ui/images/buttons/run_button_image.png", width=button_width, height=button_height)

        self.data_overview_label = None
        self.analysis_menu_window = None
        self.loading_label = None

        self.loading_image = PhotoImage(
            file="ui/images/panels/loading_panel.png")

        # create the UI elements
        self.create_file_name_label()
        self.create_file_dnd_field()
        self.create_buttons()
        self.place_images()
        self.create_visualization_canvas()

        self.mainwindow.mainloop()

    def create_file_name_label(self):
        # text indicating the currently loaded file
        self.file_name_label = Label(
            self.sidebar_frame)
        self.file_name_label.configure(
            text=self.fh.file_name, bg=default_label_color, fg=default_font_color)
        self.file_name_label.pack(pady=30, side="top")

    def create_file_dnd_field(self):
        # the field you can drop files into
        self.dnd_field = Label(self.sidebar_frame)
        self.dnd_field.configure(
            image=self.dnd_field_image, width=click_field_width, height=click_field_width, borderwidth=0)
        self.dnd_field.pack(pady=10, side="top")
        self.dnd_field.pack_propagate(False)
        self.dnd_field.drop_target_register(DND_FILES)
        self.dnd_field.dnd_bind("<<Drop>>", self.drop)

    def drop(self, event):
        # This function is called when stuff is dropped into the dnd_field
        if event.data:
            self.fh.open_file_by_dnd(event.data)

    def create_buttons(self):
        # button to browse the file explorer
        file_upload_button = Button(self.sidebar_frame)
        file_upload_button.configure(
            image=self.browse_files_button_image, width=click_field_width, height=click_field_height, borderwidth=0)
        file_upload_button.pack(pady=20, side="top")
        file_upload_button.configure(command=self.fh.open_file_by_explorer)

        # button to configure the analysis
        analysis_menu_button = Button(self.sidebar_frame)
        analysis_menu_button.configure(
            image=self.analysis_menu_button_image, width=click_field_width, height=click_field_height, borderwidth=0)
        analysis_menu_button.configure(command=self.analysis_menu)
        analysis_menu_button.pack(pady=20, side="top")

        # button to save the results
        save_result_button = Button(self.sidebar_frame)
        save_result_button.configure(
            image=self.save_results_button_image, width=click_field_width, height=click_field_height, borderwidth=0)
        save_result_button.pack(pady=20, side="top")

        # frame for the data overview and the run button
        run_frame = Frame(self.analysis_frame)
        run_frame.configure(width=200, height=200, bg=main_bg_color)
        run_frame.pack(side="top")

        # button to run the data overview
        self.overview_button = Button(run_frame)
        self.overview_button.configure(
            image=self.overview_button_image, width=click_field_width, height=click_field_height, borderwidth=0)
        self.overview_button.configure(command=self.data_overview_window)
        self.overview_button.pack(padx=10, pady=10, side="left")

        # button to run the full analysis
        self.run_button = Button(run_frame)
        self.run_button.configure(
            image=self.run_button_image, width=click_field_width, height=click_field_height, borderwidth=0)
        self.run_button.configure(command=self.run_analysis)
        self.run_button.pack(padx=10, pady=10, side="left")

        # analysis menu variables that indicate what options the user chose
        # UMAP
        self.umap_check = IntVar()
        self.umap_color_var = StringVar()
        self.umap_color_var.set("default")

    def place_images(self):
        bottom_bar = Frame(self.sidebar_frame)
        bottom_bar.configure(width=sidebar_width,
                             height=bottombar_height, bg=visualization_bg_color)
        bottom_bar.pack(side="bottom")
        bottom_bar.pack_propagate(False)

        image_claasenlab = PIL.Image.open(
            "ui/images/logos/claassen_lab_logo.png")
        image_claasenlab = image_claasenlab.resize(
            (claassenlab_image_height, claassenlab_image_height), PIL.Image.ANTIALIAS)
        image_claasenlab = ImageTk.PhotoImage(image_claasenlab)
        panel = Label(bottom_bar, image=image_claasenlab,
                      bg=Colors.ukt_greige)
        panel.image = image_claasenlab
        panel.pack(padx=5, side="left")

        image_uni_klinikum = PIL.Image.open(
            "ui/images/logos/UniklinikumTuebingen.png")
        # both should have the same height
        scale = image_uni_klinikum.height / claassenlab_image_height
        image_uni_klinikum = image_uni_klinikum.resize(
            (int(image_uni_klinikum.width / scale), int(image_uni_klinikum.height / scale)), PIL.Image.ANTIALIAS)
        image_uni_klinikum = ImageTk.PhotoImage(image_uni_klinikum)
        panel = Label(bottom_bar, image=image_uni_klinikum,
                      bg=Colors.ukt_greige)
        panel.image = image_uni_klinikum
        panel.pack(padx=5, side="right")

    def create_visualization_canvas(self):
        self.vis_canvas = Canvas(self.analysis_frame)
        self.vis_canvas.configure(width=500, height=200,
                                  bg=visualization_bg_color)
        self.vis_canvas.pack(side="bottom", fill=BOTH, expand=TRUE)

    def analysis_menu(self):
        """
        Creates a new window used to specify the analysis.
        Extend this whenever you add a new analysis functionality.
        """

        # prevents opening the window multiple times
        if self.analysis_menu_window:
            self.analysis_menu_window.destroy()

        self.analysis_menu_window = Toplevel(self.mainwindow)
        self.analysis_menu_window.title("Specify analysis")

        ws = self.analysis_menu_window.winfo_screenwidth()  # width of the screen
        hs = self.analysis_menu_window.winfo_screenheight()  # height of the screen

        # calculate x and y coordinates for the window
        x = (ws/2) - (analysis_menu_width/2) + analysis_menu_x_offset
        y = (hs/2) - (analysis_menu_height/2)

        self.analysis_menu_window.geometry(
            '%dx%d+%d+%d' % (analysis_menu_width, analysis_menu_height, x, y))
        self.analysis_menu_window.configure(
            width=analysis_menu_width, height=analysis_menu_height, bg=analysis_menu_bg_color)
        self.analysis_menu_window.resizable(False, False)

        self.umap_section()

        # TODO: extend

    def umap_section(self):
        """
        The UMAP section of the analysis menu.
        """

        umap_frame = Frame(self.analysis_menu_window)
        umap_frame.configure(
            width=analysis_menu_width - 20, height=analysis_menu_section_height, bg=default_label_color)
        umap_frame.pack_propagate(False)

        umap_checkbutton = Checkbutton(umap_frame)
        umap_checkbutton.configure(
            text="UMAP", font=analysis_menu_font, bg=analysis_menu_button_color,
            fg=analysis_menu_text_color, padx=10, variable=self.umap_check)

        self.umap_color_button = Menubutton(umap_frame)
        umap_menu = Menu(self.umap_color_button)
        self.umap_color_button.configure(menu=umap_menu)

        umap_menu.add_radiobutton(
            label="Default", variable=self.umap_color_var, value="default", command=self.update_umap_buttons)
        umap_menu.add_radiobutton(
            label="Louvain", variable=self.umap_color_var, value="louvain", command=self.update_umap_buttons)
        umap_menu.add_radiobutton(
            label="HES4", variable=self.umap_color_var, value="HES4", command=self.update_umap_buttons)
        umap_menu.add_radiobutton(
            label="TNFRSF4", variable=self.umap_color_var, value="TNFRSF4", command=self.update_umap_buttons)

        self.update_umap_buttons()

        umap_checkbutton.pack(side=LEFT)
        self.umap_color_button.pack(side=LEFT, padx=30)
        umap_frame.pack(side=TOP, pady=10)

    def update_umap_buttons(self):
        """
        Updates the menu buttons in the UMAP section of the analysis menu so that they show the selected option.
        """
        self.umap_color_button.configure(
            text="color = " + self.umap_color_var.get(), font=analysis_menu_font)

        self.mainwindow.update()

    def data_overview_window(self):
        """When the data overview button has been pressed."""

        # if there is no file loaded, notify the user and return
        if self.fh.no_file():
            self.show_message(no_file_message_type,
                              no_file_title, no_file_message)
            return

        # the frame that appears right after clicking the button
        self.do_pre_window = Frame(self.vis_canvas)
        self.do_pre_window.configure(
            bg=default_label_color, width=do_pre_window_w, height=do_pre_window_h)
        self.do_pre_window.pack_propagate(False)

        # the label showing the info text
        do_pre_window_label = Label(self.do_pre_window)
        do_pre_window_label.configure(
            bg=default_label_color, fg=default_font_color, text=do_pre_window_text, font=(standard_font, do_pre_window_font_size), pady=10)
        do_pre_window_label.pack()

        # the text input field
        text_input_field = Text(self.do_pre_window)
        text_input_field.configure(height=1, width=30, font=(
            standard_font, do_pre_window_font_size), bg=default_label_color, fg=default_font_color)
        text_input_field.pack()

        self.do_pre_window.place(relx=0.5, rely=0.5, anchor=CENTER)

        # focus the text input field right away
        text_input_field.focus_set()

        # enter starts the data overview with the current text input
        text_input_field.bind('<Return>', lambda event: self.data_overview(
            event, text_input_field.get("1.0", "end")))

    def data_overview(self, event, text_file_name: str):
        """
        Performs the data overview.
        """

        # destroy the potential old label first
        if self.data_overview_label:
            self.data_overview_label.destroy()

        try:
            self.do_pre_window.destroy()
            self.enable_loading_panel()

            do_string = self.analysis.data_overview(self.fh)

            self.disable_loading_panel()

            # add the data_overview_label
            self.data_overview_label = Label(self.vis_canvas)
            self.data_overview_label.configure(
                bg=default_label_color, text=do_string, font=(
                    standard_font, data_overview_font_size), fg=default_font_color, anchor=W, justify=LEFT)
            self.data_overview_label.pack(padx=10, pady=10, anchor=NW)

            # write the data overview to a text file with the given name
            # we have to delete the last two characters (\n) first
            last_index = len(text_file_name) - 1
            text_file_name = text_file_name[0:last_index]
            # only if input string was not empty
            if len(text_file_name) > 0:
                self.fh.write_data_overview_to_file(
                    text_file_name, do_string)
        except:
            self.show_message(data_overview_error_type,
                              data_overview_error_title, data_overview_error_message)
            self.disable_loading_panel()

    def run_analysis(self):
        """
        Manages the analysis when the user presses the "Run Analysis" button.
        """

        # if there is no file loaded, notify the user and return
        if self.fh.no_file():
            self.show_message(no_file_message_type,
                              no_file_title, no_file_message)
            return

        # manually check if there is at least one analysis option selected
        do_umap = self.umap_check.get()
        if not do_umap:  # and not everything else in the future
            self.show_message(no_analysis_selected_type,
                              no_analysis_selected_title, no_analysis_selected_message)
            return

        self.enable_loading_panel()

        # for now we only want a simple UMAP projection if the option was enabled
        try:
            if do_umap:
                self.analysis.umap(self.fh, self.umap_color_var.get())
        except:
            self.show_message(
                umap_error_type, umap_error_title, umap_error_message)
            self.disable_loading_panel()
            return

        self.disable_loading_panel()

    def enable_loading_panel(self):
        self.loading_label = Label(self.vis_canvas)
        self.loading_label.configure(image=self.loading_image)
        self.loading_label.place(relx=0.5, rely=0.5, anchor=CENTER)
        # important to ensure the update of the UI before the "heavy" methods
        self.mainwindow.update()

    def disable_loading_panel(self):
        if self.loading_label:
            self.loading_label.destroy()
            self.mainwindow.update()

    def show_message(self, type: str, title: str, message: str):
        """
        Creates a messagebox used for error messages.

        Args:
            type (str): "info", "warning" or "error"
            title (str): The title for the messagebox.
            message (str): The message to show.
        """

        if type == "warning":
            messagebox.showwarning(title, message)
        elif type == "error":
            messagebox.showerror(title, message)
        else:
            messagebox.showinfo(title, message)

    def updateUI(self):
        """
        This method is called after any event that should cause a UI change
        """

        # update the file name label and the tooltip
        self.file_name_label.config(text=self.fh.file_name)
        tooltip.CreateToolTip(self.file_name_label, text=self.fh.file_path)
