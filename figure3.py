import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.gridspec import GridSpec

def make_figure(ax, data, string, use_legend):
    data_df = pd.DataFrame(data, index = ["All", "Gain", "Neutral", "Loss"])
    ax.plot(data_df, linestyle="-", linewidth=7, marker="o", markersize=12, markeredgecolor="black", markeredgewidth=1.5)
    plt.xticks([0, 1, 2, 3], data_df.index, color="black")#, weight="bold")
    ax.set_xlabel(string, color="black", fontsize=11, weight="bold")
    ax.set_ylim([-0.05, 0.35])
    if use_legend:
        ax.legend(list(data.keys()), fancybox=True, shadow=False, bbox_to_anchor=(0.55, -0.07), ncol=5, prop=dict(weight='bold'))
    else:
        ax.set_ylabel("MCC performance", color="black", fontsize=11, weight="bold")

def main():
    data_mcc_volume = {"UEP": [0.23, 0.14, 0.23, 0.09], "pyDock": [0.23, 0.16, 0.26, 0.15], "FoldX": [0.22, 0.19, 0.16, 0.18], "PRODIGY": [0.11, 0.13, 0.06, -0.03], "BeAtMuSiC": [0.15, 0.08, 0.19, 0.10]}
    data_mcc_hydrop = {"UEP": [0.23, 0.20, 0.25, 0.27], "pyDock": [0.23, 0.19, 0.33, 0.24], "FoldX": [0.22, 0.23, 0.14, 0.20], "PRODIGY": [0.11, 0.13, 0.21, 0.06], "BeAtMuSiC": [0.15, 0.15, 0.12, 0.13]}
    
    volume_string = "Amino acid size change"
    hydrop_string = "Hydrophobicity change"

    fig = plt.figure()
    gs = GridSpec(1,2)
    ax_left = fig.add_subplot(gs[0,0])
    ax_right = fig.add_subplot(gs[0,1])
    
    make_figure(ax_left, data_mcc_volume, volume_string, False)
    make_figure(ax_right, data_mcc_hydrop, hydrop_string, True)
    fig.set_size_inches((11, 9.6), forward=False)
    fig.savefig("figure3.png", dpi=500, bbox_inches="tight")

main()
