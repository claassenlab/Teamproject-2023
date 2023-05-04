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
        filePath = fd.askopenfilename(
            title=explorer_title,
            filetypes=filetypes
        )

        # if the user cancelled
        if len(filePath) == 0:
            return

        # save the data and the last part of the file name
        self.file_data = open(filePath)
        self.file_name = os.path.basename(filePath)

        # update the UI to show the name of the loaded file
        self.ui.updateUI()
