# import the required modules
import sys
import os
import requests
import tarfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import statistics
from matplotlib.colors import rgb2hex

# figure settings
plt.rcParams['pdf.fonttype'] = 42
font = {'size'   : 12}
plt.rc('font', **font)
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.bf'] = 'Arial:bold'
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['font.serif'] = 'Arial'
plt.rcParams['font.family'] = 'Arial'


# parse the input data
def parseDocking(paths, min_cluster_size = 0.1):
    
    if isinstance(paths, str):
        paths = [paths]

    # define the dict
    cluster_stats = {}
    
    for docking in paths:

        water = f"{docking}/structures/it1/water"
        
        if not os.path.isdir(water):
            print(f"The directory {water} does not exist.")
            continue
        
        # get the number of molecules
        n_models = len([model for model in next(os.walk(water))[2] if model.endswith(".pdb")])
        
        parsed_cluster = {}
        with open(f"{water}/cluster_rmsd.txt") as clusters:
            for cluster in clusters:
                if cluster[0] == "#":
                    continue

                parsed_cluster[cluster.split()[0]] = 0
                
        for cluster in parsed_cluster:
            with open(f"{water}/{cluster}") as f:
                parsed_cluster[cluster] = len(f.readlines())


        tmp = []
        for cluster, size in parsed_cluster.items():
            
            if size <= n_models * min_cluster_size:
                continue
            
            # all_stats
            all_stats = pd.read_csv(f"{water}/{cluster}.stat", sep = " ")

            # read desolvation
            ds = pd.read_csv(f"{water}/{cluster}_Edesolv", sep = " ")

            # append ds to the all_stats df
            all_stats["Edesolv"] = all_stats.apply(lambda x: ds.loc[ds["#struc"] == x["#Structure"], "Edesolv"].item(), axis=1)

            # drop unused columns
            remove = [ 'rmsd_all', 'rmsd_Emin', 'Einter', 'Enb', 'Evdw+0.1Eelec', 'Ecdih', 'Ecoup', 'Esani', 'Evean', 'Edani', '#NOEviol', '#Dihedviol', '#Coupviol', '#Veanviol', '#Daniviol']
            all_stats = all_stats.drop(remove, axis=1) 
            
            # add cluster id
            all_stats.insert(0, 'cluster', cluster)

            #         # add the haddock score
            #        * HADDOCKscore-water =  1.0 Evdw + 0.2 Eelec + 1.0 Edesol +  0.1 Eair
            all_stats["HADDOCK_score"] = all_stats.apply(lambda x: 1.0*x.Evdw + 0.2*x.Eelec + 1.0*x.Edesolv + 0.1*x.Eair, axis = 1)

            tmp.append(all_stats)
        
        # merge all dataframes
        if len(tmp) != 0:
            cluster_stats[docking.replace(";", ".")] = pd.concat(tmp)
        else:
            cluster_stats[docking.replace(";", ".")] = None
            
    return cluster_stats


def cols_and_ticks(include_haddock_score):
    if include_haddock_score:
        cols = ["Evdw", "Eelec", "BSA", "Edesolv", "HADDOCK_score"]
        ticks = ["VdW","ES","DS", "HS"]
    else:
        cols = ["Evdw", "Eelec", "BSA", "Edesolv"]
        ticks = ["VdW","ES","DS"]
        
    return cols, ticks

def beautify_plots(axs, ticks):
    ax0, ax1 = axs
    
    # beautify
    ax0.set_xlabel("")
    ax1.set_xlabel("")

    ax0.set_ylabel("Energetics $[a.u.]$")
    ax1.set_ylabel("Buried Surface Area [Å]")
    
    # remove the legends
    legend = ax0.get_legend_handles_labels()
    
    ax0.legend_.remove()
    ax1.legend_.remove()
    
    ax0.set_xticklabels(ticks)
    
    return legend

def plot_datapoints(data, axs, hue, boxplot_args, stripplot_args):
    
    if isinstance(hue, list):
        # create a new hue parameter
        data["hue"] = data[hue].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)
    else:
        data["hue"] = data[hue]
    
    
    ax0, ax1 = axs
    # plot vdw, es, ds
    sns.boxplot(data = data.loc[data.variable != "BSA"], x = "variable", y = "value", hue = "hue", ax = ax0, **boxplot_args)
    sns.stripplot(data = data.loc[data.variable != "BSA"], x = "variable", y = "value", hue = "hue", ax = ax0,  **stripplot_args)

    # plot bsa
    sns.boxplot(data = data.loc[data.variable == "BSA"], x = "variable", y = "value", hue = "hue", ax = ax1 ,  **boxplot_args)
    sns.stripplot(data = data.loc[data.variable == "BSA"], x = "variable", y = "value", hue = "hue", ax = ax1 , **stripplot_args)
    
    return axs

def make_plot(data, axs,  hue, ticks, legend_labels, show_n, boxplot_args, stripplot_args, legend_args):
    
    # plot the data
    axs = plot_datapoints(data, axs, hue = hue, boxplot_args = boxplot_args, stripplot_args = stripplot_args)

    # beautify
    legend = beautify_plots(axs, ticks)

    # make the legend
    make_legend(axs, data, legend, legend_labels, show_n, legend_args) 
    return
    
def make_legend(axs, data, legend, legend_labels, show_n, legend_args):
    # get the sample size(s)
    #n_all_clusters = t.loc[t.variable == "BSA", "cluster"].value_counts().to_dict()

    # optimize the legend
    handles = legend[0]
    
    handles0 = handles[:len(handles) // 2]
    handles1 = handles[len(handles) // 2:]
    
    l_colors0 = len(set([rgb2hex(np.around(np.array(h.get_fc()), decimals=0)) for h in handles0]))
    l_colors1 = len(set([rgb2hex(np.around(np.array(h.get_fc()), decimals=0)) for h in handles1]))
    
    if l_colors1 >= l_colors0:
        handles = handles1
    else:
        handles = handles0
    
    # handles = handles[len(handles) // 2:]
    
    labels = legend[1]
    labels = labels[len(labels) // 2:]
    
    new_labels = []
    n_list = []

    for label in labels:
        sublabel_values = label.split(";")
        
        if len(sublabel_values) == 1:
            # single cluster naming
            n = data.loc[(data.variable == "BSA") & (data.cluster == label)].shape[0]
            new_label = "Cluster " + label.replace("file.nam_clust", "")


            
        elif len(sublabel_values) == 2:
            # docking / cluster co-naming
            docking, cluster = sublabel_values
            n = data.loc[(data.variable == "BSA") & (data.cluster == cluster) & (data.docking == docking)].shape[0]
            new_label = docking.split("/")[-1] + ": Cluster " + cluster.replace("file.nam_clust", "")

        
        if show_n:
            new_label += f" $_{{(n = {n})}}$"
            n_list.append(f" $_{{(n = {n})}}$")
            
        
        new_labels.append(new_label)
            
    if legend_labels is not None:
        # verify, that legend_labels have the correct length
        if len(legend_labels) != len(new_labels):
            print(f"Provided legend_label wrong length. Expecting {len(new_labels)} entries!")
        else:
            user_defined_labels = []
             
            for n, label in enumerate(new_labels):
                user_defined_label = legend_labels[n]
                
                if show_n:
                    user_defined_label += n_list[n]
                
                print(f"Replacing {label} with {user_defined_label}")
                
                user_defined_labels.append(user_defined_label)
            
                      
            new_labels = user_defined_labels
    # add the legend 
    # set the axis
    y_min = axs[0].get_ylim()[0]
    y_max = axs[0].get_ylim()[1]
    axs[0].set_ylim(y_max,y_min*1.2)
    
    axs[0].legend(handles, new_labels, **legend_args)

    return 

def plot_all(data,  modus = "top", include_haddock_score = False, figsize = None, legend_labels = None, show_n = True, save = False, filename = None, filetype = "png", dpi = "100", boxplot_args = {}, stripplot_args = {}, legend_args = {}, rcParams = {}):
    
    if modus not in ["top", "all"]:
        print("Error. Unkown modi. Allowed ”all” or ”top”")
        return False
    
    # prepare the data
    dockings = list(data.keys())
        
    cols, ticks = cols_and_ticks(include_haddock_score)
    
    tmp = []
    for dock in dockings:
        d = data[dock]
        d["docking"] = dock
        
        # select the respective clusters
        if modus == "top":
            most_frequent_cluster = d['cluster'].value_counts().idxmax()
            tmp.append(d.loc[d.cluster == most_frequent_cluster])
        else:
            tmp.append(d)
    
    data = pd.melt(pd.concat(tmp), id_vars = ["docking", "cluster", "#Structure"], value_vars = cols)
    
    if figsize is None:
        figsize = (5 * data.groupby(['docking', 'cluster']).ngroups * 0.7 ,5)
    with plt.rc_context(rcParams):
        fig, axs = plt.subplots(1,2, figsize= figsize,  gridspec_kw={'width_ratios': [3, 1]})
        
        # plot the data
        hue = ['docking', 'cluster']
        make_plot(data, axs, hue, ticks, legend_labels, show_n, boxplot_args, stripplot_args, legend_args)  
        
        # finalize
        sns.despine()
        plt.tight_layout()
        
        
        
        if save:
            if filename is None:
                filename = "image"
            
            plt.savefig(filename + "." + filetype, dpi = dpi)
            
        plt.show()
    
    return 

def plot_single_data(data, include_haddock_score = False, plot_single = False, figsize = None, legend_labels = None, show_n = True, save = False, filename = None, filetype = "png", dpi = 100, boxplot_args = {}, stripplot_args = {}, legend_args = {}, rcParams = {}):
    
    cols, ticks = cols_and_ticks(include_haddock_score)
    
    data = pd.melt(data, id_vars = ["cluster", "#Structure"], value_vars = cols)

    if plot_single:
        if figsize is None:
            figsize = (5 * 0.7, 5 * len(data.cluster.unique()))
        with plt.rc_context(rcParams):
            fig, axsn = plt.subplots(len(data.cluster.unique()),2, figsize= figsize,  gridspec_kw={'width_ratios': [3, 1]})
                      
            for k, cluster in enumerate(data.cluster.unique()):
                
                axs = axsn[k]
                single_data = data.loc[data.cluster == cluster].copy()
                if legend_labels is not None:
                    single_legend_label = [legend_labels[k]]
                else:
                    single_legend_label = None

                # plot the data
                hue = "cluster"
                make_plot(single_data, axs, hue, ticks, single_legend_label, show_n, boxplot_args, stripplot_args, legend_args )
                
            sns.despine()
            plt.tight_layout()
            
            if save:
                if filename is None:
                    filename = "image"
                
                plt.savefig(filename + "." + filetype, dpi = dpi)
            
            plt.show()
            
    else:
        if figsize is None:
            figsize = (5 * len(data.cluster.unique()) * 0.7 ,5)
        with plt.rc_context(rcParams):
            fig, axs = plt.subplots(1,2, figsize= figsize,  gridspec_kw={'width_ratios': [3, 1]})
    
            # plot the data
            hue = "cluster"
            make_plot(data, axs, hue, ticks, legend_labels, show_n, boxplot_args, stripplot_args, legend_args)  
        
            sns.despine()
            plt.tight_layout()
            
            
            
            if save:
                if filename is None:
                    filename = "image"
                  
                plt.savefig(filename + "." + filetype, dpi = dpi)
            plt.show()
            
    return 

def plotter(paths, plot_type="single", min_cluster_size=0.1, include_haddock_score=False,
                    plot_single=False, modus="top", figsize=None,
                    legend_labels=None, show_n=True,
                    save=False, filename=None, filetype="png", dpi=100,
                    boxplot_args=None, stripplot_args=None, legend_args=None, rcParams={}):
    """
    Parses and plots data from HADDOCK result folders specified in 'paths'.
    
    Args:
        paths: List of paths to the HADDOCK result folders.
        plot_type (str, optional): Either "single" for individual plots per path or "multi" for all data in a single plot.
                                  Defaults to "single".
        min_cluster_size (float, optional): Minimum cluster size to be considered, in proportion of total models.
                                            Defaults to 0.1.
        include_haddock_score (bool, optional): If True, includes HADDOCK score in the plot. Defaults to False.
        plot_single (bool, optional): If True and plot_type="single", each cluster gets its own subpanel in the plot.
                                      Defaults to False.
        modus (str, optional): Only applies for plot type "multi". "top" to only plot the largest clusters, "all" to plot all clusters.
                               Defaults to "top".
        palette (str or list, optional): Seaborn palette name or a list of colors with the same length as plotted clusters.
                                         Defaults to "husl".
        figsize (tuple, optional): Size of the figure in inches. Defaults to None.
        legend_labels (list, optional): List of labels for the legend. Must be the same length as plotted clusters.
                                        Defaults to None.
        show_n (bool, optional): If True, shows the n value of each cluster. Defaults to True.
        save (bool, optional): If True, saves the figure. Defaults to False.
        filename (str, optional): The name of the file to save the figure. If None, defaults to 'image.png'.
        filetype (str, optional): The type of the file to save the figure. Defaults to "png".
        dpi (int, optional): The resolution of the saved image, in dots per inch. Defaults to 100.
        boxplot_args (dict, optional): Dictionary of arguments to be passed to the seaborn boxplot function. 
                                        If None, uses default values: {"showfliers" : False, "color" : "w", "notch" : True}.
        stripplot_args (dict, optional): Dictionary of arguments to be passed to the seaborn stripplot function. 
                                         If None, uses default values: {"s" : 4,  "alpha" : 0.25, "jitter" : True, "dodge" : True}.
        legend_args (dict, optional): Dictionary of arguments to be passed to the legend. 
                                      If None, uses default values: {"frameon" : False, "loc" : "upper left", "ncol" : 1}.
        rcParams (dict, optional): Dictionary of rc parameters to be updated before plotting. Defaults to None.

    Returns:
        data: Parsed data from the HADDOCK result folders.
    """

    # parse the data:
    data = parseDocking(paths, min_cluster_size=min_cluster_size)

    # parse the plotting parameter
    
    #  default plotting parameter
    default_boxplot_args = {"showfliers" : False, "color" : "w", "notch" : True}
    default_stripplot_args = {"palette" : "husl", "s" : 4,  "alpha" : 0.25, "jitter" : True, "dodge" : True}
    default_legend_args = {"frameon" : False, "loc" : "upper left", "ncol" : 1} 
    
    if plot_single:
        default_stripplot_args.update({'palette' : 'dark:k'})
    
    
    if boxplot_args is None:
        boxplot_args = default_boxplot_args
    else:
        default_boxplot_args.update(boxplot_args)
        boxplot_args = default_boxplot_args
    
    if stripplot_args is None:
        stripplot_args = default_stripplot_args
    else:
        default_stripplot_args.update(stripplot_args)
        stripplot_args = default_stripplot_args
    
    if legend_args is None:
        legend_args = default_legend_args
    else:
        default_legend_args.update(legend_args)
        legend_args = default_legend_args
    
    if plot_type == "single":
        
        for idx, (path, parsed_data) in enumerate(data.items()):
            if save:
                modified_filename = filename + "_" + f"{idx:02d}"
            else:
                modified_filename = filename
                        
            plot_single_data(parsed_data, include_haddock_score=include_haddock_score, plot_single=plot_single,  
                figsize=figsize, legend_labels=legend_labels,  
                show_n=show_n, save=save, filename=modified_filename, filetype=filetype, dpi=dpi,
                boxplot_args=boxplot_args, stripplot_args=stripplot_args, legend_args=legend_args, rcParams=rcParams)
            
    elif plot_type == "multi":
        plot_all(data,  modus=modus, include_haddock_score=include_haddock_score, 
                 figsize=figsize, legend_labels=legend_labels,  
                 show_n=show_n, save=save, filename=filename, filetype=filetype, dpi=dpi,
                 boxplot_args=boxplot_args, stripplot_args=stripplot_args, legend_args=legend_args, rcParams=rcParams)
    else:
        raise ValueError(f"Unknown plot_type: {plot_type}")
    return data


def download_results(url, download_location=".", delete=True):
    archive = url + ".tgz"

    # Download the archive
    response = requests.get(archive, stream=True)
    download_path = os.path.join(download_location, os.path.basename(archive))
    with open(download_path, 'wb') as f:
        f.write(response.content)

    # Unpack the archive
    try:
        with tarfile.open(download_path, 'r:gz') as tar:
            tar.extractall(path=download_location)
    except tarfile.ReadError:
        with tarfile.open(download_path, 'r:') as tar:
            tar.extractall(path=download_location)

    # If delete is True, remove the tgz file
    if delete:
        os.remove(download_path)