import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

def create_figure(data):
    dataframe_data = get_rates(data)
    fig = plt.figure()
    gs = GridSpec(5,14)
    
    ax_left = fig.add_subplot(gs[0:6,0:6])
    ax_right_1, ax_right_2 = fig.add_subplot(gs[0,6:10]), fig.add_subplot(gs[0,10:14])
    ax_right_3, ax_right_4 = fig.add_subplot(gs[1,6:10]), fig.add_subplot(gs[1,10:14])
    ax_right_5, ax_right_6 = fig.add_subplot(gs[2,6:10]), fig.add_subplot(gs[2,10:14])
    ax_right_7, ax_right_8 = fig.add_subplot(gs[3,6:10]), fig.add_subplot(gs[3,10:14])
    ax_right_9, ax_right_10 = fig.add_subplot(gs[4,6:10]), fig.add_subplot(gs[4,10:14])

    make_figure_left(ax_left, dataframe_data)
    for i, axis in enumerate([ax_right_1, ax_right_2, ax_right_3, ax_right_4, ax_right_5, ax_right_6, ax_right_7, ax_right_8, ax_right_9]):
        make_confusion_matrices(axis, i, data)
    
    #fig.subplots_adjust(wspace=0.2, hspace=0.2)
    fig.set_size_inches((12, 5.5), forward=False)
    #plt.savefig("performance_plot_ultimate.png", dpi=500, bbox_inches="tight")
    
    plt.show()

def make_confusion_matrices(axis, i, data):
    font_bold_title = {"fontsize": 12, "weight": "bold"}
    font_numbers = {"fontsize": 13}
    font_labels = {"family": "monospace", "fontsize": 12}
    axis.set_axis_off()
    vertical_lines_coordinates = [0.35, 0.56, 0.77, 0.79, 0.99]
    horizontal_lines_coordinates = [0.35, 0.67, 0.99]

    for vertical in vertical_lines_coordinates:
        axis.axvline(x=vertical, ymin=vertical_lines_coordinates[0], ymax=vertical_lines_coordinates[-1], color="black")
    for horizontal in horizontal_lines_coordinates:
        axis.axhline(y=horizontal, xmin=vertical_lines_coordinates[0], xmax=vertical_lines_coordinates[2], color="black")
        axis.axhline(y=horizontal, xmin=vertical_lines_coordinates[3], xmax=vertical_lines_coordinates[4], color="black")
    
    axis.text(x=np.mean(vertical_lines_coordinates[:2]), y=1.08, s="C+", horizontalalignment="center", fontdict=font_labels) # + symbol
    axis.text(x=np.mean(vertical_lines_coordinates[1:3]), y=1.08, s="C―", horizontalalignment="center", fontdict=font_labels) # - symbol
    axis.text(x=np.mean(vertical_lines_coordinates[3:]), y=1.08, s="MCC", horizontalalignment="center", fontdict=font_labels) # MCC symbol
    axis.text(x=vertical_lines_coordinates[0]-0.10, y=np.mean(horizontal_lines_coordinates[1:]), s="P+", horizontalalignment="center", verticalalignment="center", fontdict=font_labels) # + Pred symbol
    axis.text(x=vertical_lines_coordinates[0]-0.10, y=np.mean(horizontal_lines_coordinates[:2]), s="P―", horizontalalignment="center", verticalalignment="center", fontdict=font_labels) # - Pred symbol
    
    ### Plot name ###
    plot_source = list(data.keys())[i]
    axis.text(vertical_lines_coordinates[0]-0.24, y=np.mean([horizontal_lines_coordinates[0], horizontal_lines_coordinates[-1]]), s=plot_source, horizontalalignment="center", verticalalignment="center", rotation="vertical", fontdict=font_bold_title)

    ### Plot values ###
    plot_values = data[plot_source]
    axis.text(x=np.mean(vertical_lines_coordinates[:2]), y=np.mean(horizontal_lines_coordinates[1:]), s=plot_values[0], horizontalalignment="center", verticalalignment="center", fontdict=font_numbers) # TP
    axis.text(x=np.mean(vertical_lines_coordinates[1:3]), y=np.mean(horizontal_lines_coordinates[1:]), s=plot_values[1], horizontalalignment="center", verticalalignment="center", fontdict=font_numbers) # FP
    axis.text(x=np.mean(vertical_lines_coordinates[3:5]), y=np.mean(horizontal_lines_coordinates[1:]), s=plot_values[2], horizontalalignment="center", verticalalignment="center", fontdict=font_numbers) # MCC
    axis.text(x=np.mean(vertical_lines_coordinates[:2]), y=np.mean(horizontal_lines_coordinates[:2]), s=plot_values[3], horizontalalignment="center", verticalalignment="center", fontdict=font_numbers) # FN
    axis.text(x=np.mean(vertical_lines_coordinates[1:3]), y=np.mean(horizontal_lines_coordinates[:2]), s=plot_values[4], horizontalalignment="center", verticalalignment="center", fontdict=font_numbers) # TN
    axis.text(x=np.mean(vertical_lines_coordinates[3:5]), y=np.mean(horizontal_lines_coordinates[:2]), s=plot_values[5], horizontalalignment="center", verticalalignment="center", fontdict=font_numbers) # MCC
    
def make_figure_left(ax_left, data):
    #Plot the line
    ax_left.plot(data, linestyle="-", linewidth=7, marker="o", markersize=12, markeredgecolor="black", markeredgewidth=1.5)
    #Shrink bottom axis height
    box = ax_left.get_position()
    ax_left.set_position([box.x0, box.y0 + box.height * 0.28, box.width, box.height * 0.79]) #0.26, 76
    #Formatting ticks
    plt.sca(ax_left)
    plt.xticks([0,1,2,3,4,5,6,7], data.index, rotation=90, color="black", weight="bold")
    ax_left.legend(ax_left.get_lines(), ["TPR", "TNR", "PPV", "NPV"], loc='upper center', bbox_to_anchor=(0.5, -0.32), ncol=4, fancybox=True, shadow=False, prop=dict(weight='bold'))
    plt.ylabel("Performance", color="black", fontsize=11, weight="bold")
    ax_left.yaxis.set_label_coords(-0.1,0.5)

def get_rates(data):
    rate_data = {"TPR": [], "TNR": [], "PPV": [], "NPV": []}
    for method, d in data.items():
        tp, fp, mcc, fn, tn, time = d[0], d[1], d[2], d[3], d[4], d[5]
        tpr, tnr, ppv, npv = tp/(fn+tp), tn/(fp+tn), tp/(tp+fp), tn/(fn+tn)
        rate_data["TPR"].append(tpr)
        rate_data["TNR"].append(tnr)
        rate_data["PPV"].append(ppv)
        rate_data["NPV"].append(npv)
    data_df = pd.DataFrame(rate_data, index = list(data.keys()))
    return data_df

def main():
    data = {'Strict\nconsensus': [97, 117, 0.42, 34, 312, "2h"], #changes here
            'Light\nConsensus': [212, 343, 0.31, 81, 607, "2h"],
            'UEP': [209, 416, 0.23, 85, 536, "15s"], 
            'pyDock': [204, 403, 0.23, 89, 547, "2h"], 
            'PRODIGY': [157, 394, 0.11, 131, 550, "2h"],
            'FoldX': [156, 276, 0.22, 137, 674, "2h"], 
            'BeAtMuSiC': [82, 139, 0.15, 211, 810, "Web"], 
            'mCSM\nuntrained': [16, 28, 0.05, 118, 284, "Web"],
            'mCSM\ntrained': [84, 20, 0.59, 75, 619, "Web"]}
    create_figure(data)

main()
