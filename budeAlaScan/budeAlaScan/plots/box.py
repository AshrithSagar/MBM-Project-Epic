'''
Created on 13 Jul 2017

@author: Amaurys Ávila Ibarra
'''

from matplotlib.patches import Polygon
from qtconsole.mainwindow import background

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

def plot_models_box(r_numbers, r_names, r_ddGs, pdb_id, mol_dg, model_count=1):

    '''
    Plotting ddG per residues for multiple models PDBs
    list interchangeable with tuples.
    
    This function does NOT make a call to show(). 
    
    deltaGs are meant to be passed in a list of lists where has the format
    ((R1_M1, R1_M2, R1_Mn), (Rn_M1, Rn_M2, Rn_Mn))

    r_numbers -- List of residue numbers.         
    r_names -- List of residues names.
    r_ddGs -- List of delta-deltaG for the whole molecule per residue  
    pdb_id --  Basename of the input PDB.
    mol_dg -- Molecule deltaG
    model_count -- Number of models.
    '''

#     limit_numbers = 80
#     r_numbers = r_numbers[0:limit_numbers]
#     r_names  = r_names[0:limit_numbers]
#     r_ddGs = r_ddGs[0:limit_numbers]



    y_top = 35
    y_bottom = -5

    res_count = len(r_ddGs)

    fig, ax1 = plt.subplots(figsize=(100, 60))
    ax1.patch.set_alpha(0)
    fig.patch.set_alpha(0)
    
    fig.canvas.set_window_title('%s - Residues ddG' % (pdb_id))
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = plt.boxplot(r_ddGs, notch=0, sym='+', vert=1, whis=1.5, widths=bars_width)

    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='.')

    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax1.set_axisbelow(True)

    ax1.set_title(r"Residues' $\Delta\Delta$G per Model", fontsize=title_axis_size, weight='bold')
    ax1.set_xlabel("%s     %s models     $\Delta$Gwt = %.2f kJ/mol" % (pdb_id, model_count, mol_dg), fontsize=title_axis_size, weight='bold')
    ax1.set_ylabel(r'$\Delta\Delta$G kJ/mol', fontsize=title_axis_size, weight='bold')

    medians = list(range(res_count))

    for i in range(res_count):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])

        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=bars_colour)
        ax1.add_patch(boxPolygon)

        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]

        # Plot averages
        plt.plot([np.average(med.get_xdata())], [np.average(r_ddGs[i])],
                 color='w', marker='*', markeredgecolor='k')

    # Set the axes ranges and axes labels
    ax1.set_xlim(-1, res_count + 2)
    ax1.set_ylim(y_bottom, y_top)

    ax1.tick_params(axis='y', labelsize=yticks_size)
    xtickNames = plt.setp(ax1, xticklabels=r_names)
    plt.setp(xtickNames, rotation=45, size=res_name_size)

    # Adding medians values as upperLabels, two decimal places,
    # and residue number as lowerLabels
    pos = np.arange(res_count) + 1
    upperLabels = [str(np.round(s, 2)) for s in medians]
    lowerLabels = r_numbers

    ax_ylims = ax1.get_ylim()
    top_lim = ax_ylims[1]
    bot_lim = ax_ylims[0]

    for tick, label in zip(range(res_count), ax1.get_xticklabels()):
        ax1.text(pos[tick], .9 * top_lim, upperLabels[tick],
                 horizontalalignment='center', size=labels_size, weight='semibold',
                 color=medians_colour, rotation=45)

        ax1.text(pos[tick], bot_lim - (.5 * bot_lim), lowerLabels[tick],
                 horizontalalignment='center', size=labels_size, weight='semibold',
                 color=lower_labels_colour, rotation=45)

    plt.figtext(0.90, 0.06, '*', color='grey', 
                weight='roman', size='medium')
    plt.figtext(0.915, 0.06, ' Average Value', color='black', weight='roman',
                size='x-small')
    
    plt.figtext(0.897, 0.04, '0.0', color=medians_colour,
                weight='semibold', size='small')
    plt.figtext(0.915, 0.04, ' Median Value', color='black', weight='roman',
                size='x-small')

    return


def plot_mutations_box(mut_names, mut_ddGs, pdb_id, mol_dg, model_count):

    '''
    Plotting ddG per mutation for multiple models PDBs
    list interchangeable with tuples.
    
    This function does NOT make a call to show(). 
    
    deltaGs are meant to be passed in a list of lists where has the format
    ((R1_M1, R1_M2, R1_Mn), (Rn_M1, Rn_M2, Rn_Mn))
     
    mut_names -- List of mutations names.
    mut_ddGs -- List of delta-deltaG for all mutations  
    pdb_id --  Basename of the input PDB.
    mol_dg -- Molecule deltaG
    model_count -- Number of models.
    '''
    a_max = 0

    for a_mut in mut_ddGs:
        mut_max = max(a_mut)
        if a_max < mut_max:
            a_max = mut_max

    y_top = int(a_max) + 10
    y_bottom = -10

    mut_count = len(mut_ddGs)

    fig, ax = plt.subplots(figsize=(100, 60))
    ax.patch.set_alpha(0)
    fig.patch.set_alpha(0)
    fig.canvas.set_window_title('%s - Residues ddG' % (pdb_id))
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = plt.boxplot(mut_ddGs, notch=0, sym='+', vert=1, whis=1.5, widths=bars_width)

    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='.')

    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax.set_axisbelow(True)

    ax.set_title(r"Mutations' $\Delta\Delta$G per Model", fontsize=title_axis_size, weight='bold')
    ax.set_xlabel("%s     %s models     $\Delta$Gwt = %.2f kJ/mol" % (pdb_id, model_count, mol_dg), fontsize=title_axis_size, weight='bold')
    ax.set_ylabel(r'$\Delta\Delta$G kJ/mol', fontsize=title_axis_size, weight='bold')

    medians = list(range(mut_count))

    for i in range(mut_count):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])

        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=bars_colour)
        ax.add_patch(boxPolygon)

        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]

        # Plot averages
        plt.plot([np.average(med.get_xdata())], [np.average(mut_ddGs[i])],
                 color='w', marker='*', markeredgecolor='k')

    # Set the axes ranges and axes labels
    ax.set_xlim(-1, mut_count + 2)
    ax.set_ylim(y_bottom, y_top)

    ax.tick_params(axis='y', labelsize=yticks_size)
    xtickNames = plt.setp(ax, xticklabels=mut_names)
    plt.setp(xtickNames, rotation=45, size=mut_name_size)

    # Adding medians values as upperLabels, two decimal places,
    # and residue number as lowerLabels
    pos = np.arange(mut_count) + 1
    upperLabels = [str(np.round(s, 2)) for s in medians]
#     lowerLabels = r_numbers

    ax_ylims = ax.get_ylim()
    top_lim = ax_ylims[1]
#     bot_lim = ax_ylims[0]

    for tick, label in zip(range(mut_count), ax.get_xticklabels()):
        ax.text(pos[tick], .98 * top_lim, upperLabels[tick],
                 horizontalalignment='center', size=labels_size, weight='semibold',
                 color=medians_colour, rotation=45)

    plt.figtext(0.90, 0.06, '*', color='grey', 
                weight='roman', size='medium')
    plt.figtext(0.915, 0.06, ' Average Value', color='black', weight='roman',
                size='x-small')
    
    plt.figtext(0.897, 0.04, '0.0', color=medians_colour,
                weight='semibold', size='small')
    plt.figtext(0.915, 0.04, ' Median Value', color='black', weight='roman',
                size='x-small')

    return
