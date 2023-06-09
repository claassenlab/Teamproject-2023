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
        self.avg_non_0_signal_genes = None
        self.min_num_non_0_genes = None
        self.min_cell_non_0_genes = None
        self.max_num_non_0_genes = None
        self.max_cell_non_0_genes = None

    def update(self, fh: FileHandler):
        """
        When requesting analysis, this method checks if the stored data is still
        matching the current requested file path. If that is not the case,
        it needs to create a new AnnData object by reading the .loom file.

        Args:
            fh (FileHandler): The file handler storing the currently loaded file path.

        Returns:
            Whether one needs to update the stored data because the file is new.
        """

        # get the path of the currently loaded file
        current_path = fh.get_loom()

        # if the current requested path is not the current path of this Analysis object
        if current_path != self.file_path:
            # we need to read in a new AnnData object from the loom file
            self.adata = sc.read_loom(current_path)
            # and update the file path of this object
            self.file_path = current_path
            return True

        return False

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

        maybe_upate = self.update(fh)

        # if there is no data
        if self.no_data():
            return None

        # each row corresponds to a cell
        n_cells = self.adata.n_obs

        # each column corresponds to a gene
        n_genes = self.adata.n_vars

        # only do costly operations when there is an update required
        if (maybe_upate):

            total_non_0_genes = 0

            # get the cell with the minimum number of non-zero signal genes
            self.min_num_non_0_genes = n_genes + 1
            self.min_cell_non_0_genes = None

            # get the cell with the maximum number of non-zero signal genes
            self.max_num_non_0_genes = -1
            self.max_cell_non_0_genes = None

            for i in range(n_cells):
                # count the non-zero signal genes in the current cell
                non_0_genes_cell = 0
                for j in range(n_genes):
                    if (self.adata.X[i, j] != 0):
                        non_0_genes_cell += 1

                # check if the min/max cell needs to be updated
                if non_0_genes_cell < self.min_num_non_0_genes:
                    self.min_num_non_0_genes = non_0_genes_cell
                    self.min_cell_non_0_genes = self.adata.obs_names[i]
                if non_0_genes_cell > self.max_num_non_0_genes:
                    self.max_num_non_0_genes = non_0_genes_cell
                    self.max_cell_non_0_genes = self.adata.obs_names[i]

                # update total non-zero signal genes count
                total_non_0_genes += non_0_genes_cell

            # get the average amount of non-zero genes per cell by counting how many
            # genes have a non-zero signal and dividing by the amount of cells
            self.avg_non_0_signal_genes = total_non_0_genes/n_cells

        output = ""
        output += "Dataset overview:" + "\n"
        output += "----------------------------------------------------------------------------------------------------------" + "\n"
        output += "Number of cells: " + str(n_cells) + "\n"
        output += "Number of genes: " + str(n_genes) + "\n"
        output += "Average amount of non-zero signal genes: " + \
            str(self.avg_non_0_signal_genes) + "\n"
        output += "Cell with fewest non-zero signal genes: " + \
            str(self.min_cell_non_0_genes) + \
            " (" + str(self.min_num_non_0_genes) + " genes)" + "\n"
        output += "Cell with most non-zero signal genes: " + \
            str(self.max_cell_non_0_genes) + \
            " (" + str(self.max_num_non_0_genes) + " genes)"

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
            sc.pl.umap(self.adata, color=col)
