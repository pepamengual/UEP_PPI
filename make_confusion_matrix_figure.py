import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def create_figure(data):
    fig, ax = plt.subplots(4, 2)
    vertical_lines_coordinates = [0.25, 0.48, 0.71, 0.73, 0.99]
    horizontal_lines_coordinates = [0.25, 0.62, 0.99]
    
    mono = {'family': 'monospace', 'fontsize': 13}
    for axis in ax.flat:
        for vertical in vertical_lines_coordinates:
            axis.axvline(x=vertical, ymin=vertical_lines_coordinates[0], ymax=vertical_lines_coordinates[-1], color="black")
        for horizontal in horizontal_lines_coordinates:
            axis.axhline(y=horizontal, xmin=vertical_lines_coordinates[0], xmax=vertical_lines_coordinates[2], color="black")
            axis.axhline(y=horizontal, xmin=vertical_lines_coordinates[3], xmax=vertical_lines_coordinates[4], color="black")
        
        axis.text(x=np.mean(vertical_lines_coordinates[:2]), y=1.05, s="+", horizontalalignment='center', fontdict=mono) # + symbol
        axis.text(x=np.mean(vertical_lines_coordinates[1:3]), y=1.05, s="-", horizontalalignment='center', fontdict=mono) # - symbol
        axis.text(x=np.mean(vertical_lines_coordinates[3:]), y=1.05, s="MCC", horizontalalignment='center', fontdict=mono) # MCC symbol
        axis.text(x=0.15, y=np.mean(horizontal_lines_coordinates[1:]), s="+ Pred", horizontalalignment='center', verticalalignment='center', fontdict=mono) # + Pred symbol
        axis.text(x=0.15, y=np.mean(horizontal_lines_coordinates[:2]), s="- Pred", horizontalalignment='center', verticalalignment='center', fontdict=mono) # - Pred symbol
        axis.text(x=0, y=np.mean([horizontal_lines_coordinates[0], horizontal_lines_coordinates[-1]]), s="UEP", horizontalalignment='center', verticalalignment='center', fontdict=mono, rotation="vertical")
    [axis.set_axis_off() for axis in ax.ravel()] #set all axis off
    plt.show()

def main():
    data = [10, 9, 8, 7]
    create_figure(data)

main()
