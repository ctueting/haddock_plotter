# Haddock Plotter

The Haddock Plotter is a Python script used to visualize and analyze HADDOCK result folders. HADDOCK (High Ambiguity Driven biomolecular DOCKing) is an information-driven flexible docking approach for the modeling of biomolecular complexes.

## Description

This script parses data from HADDOCK result folders and plots the information in an organized manner. The plots include information about the docking results such as Van der Waals energy, Electrostatic energy, and the Buried Surface Area. Furthermore, it is possible to include the HADDOCK score.

## Dependencies

- Python (version 3.6 or higher)
- Modules:
    - pandas
    - matplotlib
    - seaborn
    - numpy
    - os
    - sys
    - math
    - statistics

## Usage

### Import

Firstly, import the script:

\```python
import haddock_plotter
\```

### Running

To run the `haddock_plotter`, use the following command:

\```python
haddock_plotter.haddock_plotter(paths, plot_type="single", min_cluster_size=0.1, include_haddock_score=False,
                                plot_single=False, modus="top", palette="husl", figsize=None,
                                close=False, legend_labels=None, ncol=1, show_n=True,
                                save=False, filename=None, filetype="png", dpi=100)
\```

The parameters are as follows:

- `paths`: List of paths to the HADDOCK result folders.
- `plot_type`: Either "single" for individual plots per path or "multi" for all data in a single plot.
- `min_cluster_size`: Minimum cluster size to be considered, in proportion of total models.
- `include_haddock_score`: If True, includes HADDOCK score in the plot.
- `plot_single`: If True and `plot_type="single"`, each cluster gets its own subpanel in the plot.
- `modus`: Only applies for plot type "multi". "top" to only plot the largest clusters, "all" to plot all clusters.
- `palette`: Seaborn palette name or a list of colors with the same length as plotted clusters.
- `figsize`: Size of the figure in inches.
- `close`: If True, closes the figure after plotting (not recommended).
- `legend_labels`: List of labels for the legend. Must be the same length as plotted clusters.
- `ncol`: Number of columns in the legend.
- `show_n`: If True, shows the n value of each cluster.
- `save`: If True, saves the figure.
- `filename`: The name of the file to save the figure. If None, defaults to 'image.png'.
- `filetype`: The type of the file to save the figure.
- `dpi`: The resolution of the saved image, in dots per inch.

## Contribution

Feel free to contribute to this project by providing bug reports, feature requests, or code improvements via pull requests. Please ensure that your code passes existing unit tests and, if applicable, add new tests for your features.


This README was written by ChatGPT, based on the original code and a fruitful conversation.
