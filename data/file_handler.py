from tkinter import filedialog as fd
from io import TextIOWrapper
import os


explorer_title = "Select a file to open"
filetypes = [("TXT (only testing)", "*.txt"),
             ("BAM", "*.bam"), ("FASTQ", "*.fastq")]


class FileHandler:
    """
    Handles file opening, loading and processing

    Args:
        ui (UI): The user interface instance of the program
    """

    def __init__(self, ui):
        self.file_data: TextIOWrapper = None
        self.file_name: str = "No file opened!"
        self.ui = ui

    def open_file_by_explorer(self):
        """
        Opens the file explorer and lets the user select a file.
        Saves the corresponding TextIOWrapper for further processing.
        """

        # open the file explorer and let the user select a file
        file_path = fd.askopenfilename(
            title=explorer_title,
            filetypes=filetypes
        )

        # if the user cancelled
        if len(file_path) == 0:
            return

        # save the data and the last part of the file name
        self.file_data = open(file_path)
        self.file_name = os.path.basename(file_path)

        # update the UI to show the name of the loaded file
        self.ui.updateUI()

    def open_file_by_dnd(self, file_path):
        """
        This method is called when there is a drag-and-drop event.
        It checks whether the file type is valid and processes the file
        if that's the case.
        """

        # exclude brackets
        file_path = file_path[1:len(file_path)-1]

        # get the file type
        file_extension = os.path.splitext(file_path)[1]

        # check whether it matches a valid one
        for t in filetypes:
            # exclude the * from t[1]
            t = t[1][1:len(t[1])]

            if t == file_extension:
                self.file_name = os.path.basename(file_path)
                self.ui.updateUI()
                return
