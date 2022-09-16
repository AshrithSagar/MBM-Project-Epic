'''
Created on 13 Jul 2017

@author: Amaurys Ávila Ibarra
'''
import matplotlib.pyplot as plt
import numpy as np


bars_colour="cadetblue"
medians_colour="midnightblue"
title_axis_size='x-large'
yticks_size='large'
labels_size='medium'
res_name_size='medium'
mut_name_size='medium'
lower_labels_colour = 'royalblue'
bars_width=0.2

def plot_model_bar(r_numbers, r_names, r_ddGs, sd, pdb_id, mol_dg, model_count=1, is_avg=False):

    '''
    Plotting ddG per residues for single model PDBs or
    average for multiple models PDBs list interchangeable with tuples.
    
    This function does NOT make a call to show(). 
    
    deltaGs are meant to be passed in a list of all ddG per residue
    
    r_numbers -- List of residue numbers.         
    r_names -- List of residues names.
    r_ddGs -- List of delta-deltaG for the whole molecule per residue  
    sd -- List of standard deviation per residue.
    pdb_id --  Basename of the input PDB.
    mol_dg -- Molecule deltaG
    model_count -- Number of models.
    is_avg -- Whether we are using average ddGs from multiple models or not.
    '''

    y_top = 35
    y_bottom = -5
    res_count = len(r_ddGs)
    x_values = np.arange(1, res_count + 1)
    fig, ax = plt.subplots(figsize=(100, 60))

    ax.patch.set_alpha(0)
    fig.patch.set_alpha(0)
    fig.canvas.set_window_title('%s - Residues ddG' % (pdb_id))
    
#     ax.autoscale_view(tight=False,scalex=False,scaley=True)
    plt.subplots_adjust(left=0.065, right=0.95, top=0.9, bottom=0.25)

    if is_avg:
        res_bars = ax.bar(x_values, r_ddGs, bars_width, color=bars_colour, align='center', yerr=sd)
        ax.set_title(r"Residues' Average $\Delta\Delta$G", fontsize=title_axis_size, weight='bold')
        ax.set_xlabel("%s     %s models     $\Delta$Gwt = %.2f kJ/mol" % (pdb_id, model_count, mol_dg), fontsize=title_axis_size, weight='bold')
    else:
        res_bars = ax.bar(x_values, r_ddGs, bars_width, color=bars_colour, align='center')
        ax.set_title(r"Residues' $\Delta\Delta$G", fontsize=title_axis_size, weight='bold')
        ax.set_xlabel("%s     $\Delta$Gwt = %.2f kJ/mol" % (pdb_id, mol_dg), fontsize=title_axis_size, weight='bold')

    ax.set_ylim(y_bottom, y_top)
    ax.set_xlim(-1, res_count + 2)

    ax.set_ylabel(r'$\Delta\Delta$G  kJ/mol', fontsize=title_axis_size, weight='bold')

    ax.tick_params(axis='y', labelsize=yticks_size)
    
    ax.set_xticks(x_values)
    xtickNames = plt.setp(ax, xticklabels=r_names)
    plt.setp(xtickNames, rotation=45, size=res_name_size)

    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax.set_axisbelow(True)

    upperLabels = [str(np.round(s, 2)) for s in r_ddGs]
    set_labels(res_bars, ax, (x_values - bars_width / 6), r_numbers, upperLabels)
    
    if is_avg:
        plt.figtext(0.90, 0.06, '|', color='black', 
                    weight='roman', size='medium')
        plt.figtext(0.915, 0.058, ' Standard Deviation', color='black',
                     weight='roman', size='x-small')

    return


def set_labels(res_bars, ax, x_values, lowerLabels, upperLabels):
    '''
    Set the lower and upper labels in the graph.
    res_bars -- Bar for the residues.
    ax -- Axes object
    x_values -- X coordinates for lowerLabels and upperLabels.
    lowerLabels -- Labels to show at the bottom of the graph, residue number (B46)
    upperLabels -- Labels to show at the top of the graph, ddGs for the residues.
    '''

    ax_ylims = ax.get_ylim()
    top_lim = ax_ylims[1]
    bot_lim = ax_ylims[0]

    for i in range(len(res_bars)):
        ax.text(x_values[i], bot_lim - (.5 * bot_lim), lowerLabels[i], size='medium', weight='semibold',
                color=lower_labels_colour, rotation=45)
        ax.text(x_values[i], .9 * top_lim, upperLabels[i], horizontalalignment='center', size='medium',
                weight='semibold', color=bars_colour, rotation=45)


def plot_mutations_bar(mut_names, mut_ddGs, pdb_id, mol_dg, model_count):

    '''
    Plotting ddG per mutation for single model PDBs or
    average for multiple models PDBs list interchangeable with tuples.
    
    This function does NOT make a call to show(). 
    
    deltaGs are meant to be passed in a list of all ddG per residue
       
    mut_names -- List of mutations names.
    mut_ddGs -- List of delta-deltaG for mutations  
    pdb_id --  Basename of the input PDB.
    mol_dg -- Molecule deltaG
    model_count -- Number of models.
    '''
    mut_count = len(mut_ddGs)
    x_values = np.arange(1, mut_count + 1)
    fig, ax = plt.subplots(figsize=(100, 60))
    
    ax.patch.set_alpha(0)
    fig.patch.set_alpha(0)

    fig.canvas.set_window_title('%s - Residues ddG' % (pdb_id))
    plt.subplots_adjust(left=0.065, right=0.95, top=0.9, bottom=0.25)

    mut_bars = ax.bar(x_values, mut_ddGs, bars_width, color=bars_colour, align='center')
    ax.set_title(r"Mutations' $\Delta\Delta$G", fontsize=title_axis_size, weight='bold')
    ax.set_xlabel("%s     %s models     $\Delta$Gwt = %.2f kJ/mol" % (pdb_id, model_count, mol_dg), fontsize=title_axis_size, weight='bold')

    y_top = int(max(mut_ddGs)) + 10
    y_bottom = -5
    ax.set_ylim(y_bottom, y_top)
    ax.set_xlim(-1, mut_count + 2)

    ax.set_ylabel(r'$\Delta\Delta$G kJ/mol', fontsize=title_axis_size, weight='bold')

    ax.tick_params(axis='y', labelsize=yticks_size)
    ax.set_xticks(x_values)
    xtickNames = plt.setp(ax, xticklabels=mut_names)
    plt.setp(xtickNames, rotation=45, fontsize=mut_name_size)

    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax.set_axisbelow(True)

    upperLabels = [str(np.round(s, 2)) for s in mut_ddGs]

    lx_values = x_values - bars_width / 6

    for i in range(len(mut_bars)):
#         ax.text(lx_values[i],  y_bottom - (.5*y_bottom), upperLabels[i], size='x-small', weight='semibold',
#             color='royalblue',rotation=45, fontsize=8)
        ax.text(lx_values[i], .95 * y_top, upperLabels[i], horizontalalignment='center', size='medium',
            weight='bold', color=bars_colour, rotation=45)
        

    return

