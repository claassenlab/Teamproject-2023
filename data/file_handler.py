from tkinter import filedialog as fd
from io import TextIOWrapper
import os
from tkinter import Listbox


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
        self.file_path = None
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

        self.file_path = file_path

        # save the data and the last part of the file name
        self.file_data = open(file_path)
        self.file_name = os.path.basename(file_path)

        # update the UI to show the name of the loaded file
        self.ui.updateUI()

    def open_file_by_dnd(self, file_paths):
        """
        This method is called when there is a drag-and-drop event.
        It checks whether the file type is valid and processes the file
        if that's the case.
        """

        # convert the file paths to a list using a function from Listbox
        listbox = Listbox()
        files = listbox.tk.splitlist(file_paths)

        # pick the first file path, even if the user accidentally dropped multiple files
        file_path = files[0]

        # get the file type
        file_extension = os.path.splitext(file_path)[1]

        # if the file name is empty, we have to get the file extension in a dirty way
        # but in general, it is not that bad if one cannot load such files via dnd
        if len(file_extension) == 0:
            for i in reversed(range(len(file_path))):
                currentChar = file_path[i]
                # for Windows this works, maybe add alternatives for other OS
                if currentChar == "/":
                    break
                file_extension = file_path[i] + file_extension

        # check whether the extension matches a valid one
        for t in filetypes:
            # exclude the *
            t = t[1][1:len(t[1])]

            if t == file_extension:
                self.file_path = file_path
                self.file_name = os.path.basename(file_path)
                self.ui.updateUI()
                return
