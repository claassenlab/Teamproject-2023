import scanpy as sc
from data.file_handler import FileHandler


# Please note: Converting and reading needs a lot of time. Therefore, it is a good idea
# to use already converted/read data again and again, as long as it is still requested.
# This class follows this principle as it covers the whole analysis and only updates
# the AnnData object if the user changed the file.

class Analysis:
    """
    Does all the analysis when the user requests it.
    """

    def __init__(self):
        # no data yet
        self.file_path = None
        self.adata = None

    def update(self, fh: FileHandler):
        """
        When requesting analysis, this method checks if the stored data is still
        matching the current requested file path. If that is not the case,
        it needs to create a new AnnData object by reading the .loom file.

        Args:
            fh (FileHandler): The file handler storing the currently loaded file path.
        """

        # get the path of the currently loaded file
        current_path = fh.get_loom()

        # if the current requested path is not the current path of this Analysis object
        if current_path != self.file_path:
            # we need to read in a new AnnData object from the loom file
            self.adata = sc.read_loom(current_path)
            # and update the file path of this object
            self.file_path = current_path

    def no_data(self):
        """Returns whether there is data to process."""
        return self.file_path == None

    def data_overview(self, fh: FileHandler):
        """
        Computes the information to display when the user requests a data overview.
        Creates a string out of it which is then used by the UI.

        Args:
            fh (FileHandler): The file handler storing the currently loaded file.

        Returns:
            string: A multi-line string which the UI uses to display the information.
        """

        self.update(fh)

        # if there is no data
        if self.no_data():
            return None

        # each row corresponds to a cell
        n_cells = self.adata.n_obs

        # each column corresponds to a gene
        n_genes = self.adata.n_vars

        output = ""
        output += "Dataset overview:" + "\n"
        output += "----------------------------------------" + "\n"
        output += "Number of cells: " + str(n_cells) + "\n"
        output += "Number of genes: " + str(n_genes)

        return output

    def umap(self, fh: FileHandler, col: str):
        """
        Creates the UMAP projections.
        """

        self.update(fh)

        if self.no_data():
            return

        # create the UMAP projection from the loaded and converted file
        sc.pp.neighbors(self.adata)
        sc.tl.umap(self.adata)
        if col == "default":
            sc.pl.umap(self.adata)
        else:
            # Our example file does not work for colors other than default
            # so let us take a working dataset for demonstrational purposes
            adata = sc.datasets.pbmc68k_reduced()
            sc.pl.umap(adata, color=col)
