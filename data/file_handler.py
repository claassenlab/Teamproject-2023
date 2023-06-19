from tkinter import filedialog as fd
from io import TextIOWrapper
import os
from tkinter import Listbox


explorer_title = "Select a file to open"

# TODO: remove .txt later, for now only loom files can be processed
# filetypes = [("Loom", "*.loom"), ("TXT (only testing)", "*.txt"),
# ("BAM", "*.bam"), ("FASTQ", "*.fastq")]
filetypes = [("Loom", "*.loom")]


class FileHandler:
    """
    Handles file opening, loading and processing

    Args:
        ui (UI): The user interface instance of the program
    """

    def __init__(self, ui):
        # saves the current file path and extension
        self.file_path = None
        self.file_extension = None
        self.ui = ui

        # saves the file name + extension of the currently loaded file
        # default message when no file loaded
        # TODO: remove this line and uncomment the one below
        self.file_name: str = "Only .loom files for now!"
        # self.file_name: str = "No file opened!"

    def no_file(self):
        """Whether there is a file loaded right now."""
        return self.file_path == None

    def open_file_by_explorer(self):
        """
        Opens the file explorer and lets the user select a file.
        """

        # open the file explorer and let the user select a file
        file_path = fd.askopenfilename(
            title=explorer_title,
            filetypes=filetypes
        )

        # if the user cancelled
        if len(file_path) == 0:
            return

        # save the file path and extension
        self.file_path = file_path
        self.file_extension = self.get_extension(file_path)

        # save the file name and update the UI label
        self.file_name = os.path.basename(file_path)
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

        # get the file extension
        file_extension = self.get_extension(file_path)

        # check whether the extension matches a valid one
        for t in filetypes:
            # exclude the *
            t = t[1][1:len(t[1])]

            if t == file_extension:
                # save the file path and extension
                self.file_path = file_path
                self.file_extension = file_extension

                # save the file name and update the UI label
                self.file_name = os.path.basename(file_path)
                self.ui.updateUI()
                return

    def get_extension(self, file_path: str):
        """
        Gets the file extension of the provided file path.

        Args:
            file_path (str): The path of the file.

        Returns:
            string: The file extension of the file, even for nameless files.
        """

        file_extension = os.path.splitext(file_path)[1]

        # if the file name is empty (e. g. '.bam'), we still want to get the file
        # extension so that those files can be processed -> do it in a dirty way
        if len(file_extension) == 0:
            for i in reversed(range(len(file_path))):
                currentChar = file_path[i]
                # for Windows this works, maybe add alternatives for other OS
                if currentChar == "/":
                    break
                file_extension = file_path[i] + file_extension

        return file_extension

    def get_loom(self):
        """
        Converts the currently loaded file to a loom file if it isn't already one.

        Returns:
            string: The file path of the corresponding loom file.
            None: If no file is currently loaded.
        """

        # check the extension of the currently loaded file
        if self.file_extension == ".bam":
            # convert it to .loom using velocyto
            # set self.file_path to this loom file path
            # set self.file_extension to ".loom"
            pass
        if self.file_extension == ".fastq":
            # convert it to .loom using velocyto
            # set self.file_path to this loom file path
            # set self.file_extension to ".loom"
            pass
        # if the loaded file is a .loom file, return the path
        if self.file_extension == ".loom":
            return self.file_path

        # return None if there is no currently loaded file
        return None

    def write_data_overview_to_file(self, file_name: str, do: str):
        """
        Writes the data overview string to a text file with the given name.
        Location: user/data_overview

        Args:
            file_name (str): The name of the file to be created (without .txt)
            do (str): The full, multi-line data overview string
        """

        f = open("user/data_overview/" + file_name + ".txt", "w")
        # write the full file path to the text file as well
        f.write("File: " + self.file_path + "\n\n" + do)
        f.close()
