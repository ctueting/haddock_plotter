# import the required modules
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import statistics


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

def parse_palette(palette):

    if isinstance(palette, str):
        return palette
    
    # Define your own color palette
    return sns.color_palette(palette)


def plot_datapoints(data, axs, palette, hue):
    
    if isinstance(hue, list):
        # create a new hue parameter
        data["hue"] = data[hue].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)
    else:
        data["hue"] = data[hue]
    
    
    palette = parse_palette(palette)
    
    ax0, ax1 = axs
    # plot vdw, es, ds
    sns.boxplot(data = data.loc[data.variable != "BSA"], x = "variable", y = "value", hue = "hue", ax = ax0, showfliers=False, color = "w", notch=True)
    sns.stripplot(data = data.loc[data.variable != "BSA"], x = "variable", y = "value", hue = "hue", ax = ax0, s=4, palette=palette, alpha = 0.25, jitter=True, dodge = True)

    # plot bsa
    sns.boxplot(data = data.loc[data.variable == "BSA"], x = "variable", y = "value", hue = "hue", ax = ax1 , showfliers=False, color = "w", notch=True)
    sns.stripplot(data = data.loc[data.variable == "BSA"], x = "variable", y = "value", hue = "hue", ax = ax1 , s=4, palette=palette, alpha = 0.25, jitter=True, dodge = True)
    
    return axs

def make_plot(data, axs, palette, hue, ticks, legend_labels, show_n, ncol):
    # plot the data
    axs = plot_datapoints(data, axs, palette = palette, hue = hue)

    # beautify
    legend = beautify_plots(axs, ticks)

    # make the legend
    make_legend(axs, data, legend, legend_labels, show_n, ncol) 
    return
    
def make_legend(axs, data, legend, legend_labels, show_n, ncol):
    # get the sample size(s)
    #n_all_clusters = t.loc[t.variable == "BSA", "cluster"].value_counts().to_dict()

    # optimize the legend
    handles = legend[0]
    handles = handles[len(handles) // 2:]
    
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
    
    axs[0].legend(handles, new_labels, frameon = False, loc = "upper left", ncol = ncol)

    return 

def plot_all(data,  modus = "top", include_haddock_score = False, palette = "husl", figsize = None, close = False, legend_labels = None, ncol = 1, show_n = True, save = False, filename = None, filetype = "png", dpi = "100"):
    
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
    
    fig, axs = plt.subplots(1,2, figsize= figsize,  gridspec_kw={'width_ratios': [3, 1]})
    
    # plot the data
    hue = ['docking', 'cluster']
    make_plot(data, axs, palette, hue, ticks, legend_labels, show_n, ncol)  
    
    # finalize
    sns.despine()
    plt.tight_layout()
    
    if save:
        if filename is None:
            filename = "image"
        
        plt.savefig(filename + "." + filetype, dpi = dpi)
        
    if close:
        plt.close()
    
    return 

def plot_single_data(data, include_haddock_score = False, plot_single = False, palette = "husl", figsize = None, close = False, legend_labels = None, ncol = 1, show_n = True, save = False, filename = None, filetype = "png", dpi = "100"):
    
    cols, ticks = cols_and_ticks(include_haddock_score)
    
    data = pd.melt(data, id_vars = ["cluster", "#Structure"], value_vars = cols)

    if plot_single:
        if figsize is None:
            figsize = (5 * 0.7, 5 * len(data.cluster.unique()))
        
        fig, axsn = plt.subplots(len(data.cluster.unique()),2, figsize= figsize,  gridspec_kw={'width_ratios': [3, 1]})
        
        palette = ["black"]
        
        for k, cluster in enumerate(data.cluster.unique()):
            axs = axsn[k]
            single_data = data.loc[data.cluster == cluster].copy()
            if legend_labels is not None:
                single_legend_label = [legend_labels[k]]
            else:
                single_legend_label = None

            # plot the data
            hue = "cluster"
            make_plot(single_data, axs, palette, hue, ticks, single_legend_label, show_n, ncol)
    else:
        if figsize is None:
            figsize = (5 * len(data.cluster.unique()) * 0.7 ,5)

        fig, axs = plt.subplots(1,2, figsize= figsize,  gridspec_kw={'width_ratios': [3, 1]})

        # plot the data
        hue = "cluster"
        make_plot(data, axs, palette, hue, ticks, legend_labels, show_n, ncol)  
    
    sns.despine()
    plt.tight_layout()
    
    if save:
        if filename is None:
            filename = "image"
        
        plt.savefig(filename + "." + filetype, dpi = dpi)
        
    if close:
        plt.close()
    

    
    return 

def plotter(paths, plot_type="single", min_cluster_size=0.1, include_haddock_score=False,
                    plot_single=False, modus="top", palette="husl", figsize=None,
                    close=False, legend_labels=None, ncol=1, show_n=True,
                    save=False, filename=None, filetype="png", dpi=100):
    """
    Parses and plots data from HADDOCK result folders specified in 'paths'.
    
    Args:
        paths: List of paths to the HADDOCK result folders.
        plot_type (str): Either "single" for individual plots per path or "multi" for all data in a single plot.
        min_cluster_size (float): Minimum cluster size to be considered, in proportion of total models.
        include_haddock_score (bool): If True, includes HADDOCK score in the plot.
        plot_single (bool): If True and plot_type="single", each cluster gets its ownprint subpanel in the plot.
        modus (str): Only applies for plot type "multi". "top" to only plot the largest clusters, "all" to plot all clusters.
        palette (str or list): Seaborn palette name or a list of colors with the same length as plotted clusters.
        figsize (tuple): Size of the figure in inches.
        close (bool): If True, closes the figure after plotting (not recommended).
        legend_labels (list): List of labels for the legend. Must be the same length as plotted clusters.
        ncol (int): Number of columns in the legend.
        show_n (bool): If True, shows the n value of each cluster.
        save (bool): If True, saves the figure.
        filename (str): The name of the file to save the figure. If None, defaults to 'image.png'.
        filetype (str): The type of the file to save the figure.
        dpi (int): The resolution of the saved image, in dots per inch.

    Returns:
        data: Parsed data from the HADDOCK result folders.
    """

    # parse the data:
    data = parseDocking(paths, min_cluster_size=min_cluster_size)
    
    if plot_type == "single":
        
        for path, parsed_data in data.items():
            if save:
                modified_filename = path + "_" + filename
            else:
                modified_filename = filename
                
            plot_single_data(parsed_data, include_haddock_score=include_haddock_score, plot_single=plot_single, palette=palette, 
                figsize=figsize, close=close, legend_labels=legend_labels, ncol=ncol, 
                show_n=show_n, save=save, filename=modified_filename, filetype=filetype, dpi=dpi)
            
    elif plot_type == "multi":
        plot_all(data,  modus=modus, include_haddock_score=include_haddock_score, palette=palette, 
                 figsize=figsize, close=close, legend_labels=legend_labels, ncol=ncol, 
                 show_n=show_n, save=save, filename=filename, filetype=filetype, dpi=dpi)
    else:
        raise ValueError(f"Unknown plot_type: {plot_type}")
        
    return data