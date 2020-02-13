import matplotlib.pyplot as plt
import matplotlib.patches as patches

def create_figure(data):
    fig, ax = plt.subplots(4, 2)
    
    ax[0,0].axvline(x=0.01, ymax=0.98, color="black")
    ax[0,0].axvline(x=0.49, ymax=0.98, color="black")
    ax[0,0].axvline(x=0.99, ymax=0.98, color="black")
    ax[0,0].axhline(y=0.01, xmax=0.98, color="black")
    ax[0,0].axhline(y=0.49, xmax=0.98, color="black")
    ax[0,0].axhline(y=0.99, xmax=0.98, color="black")

    #ax[0,0].add_patch(patches.Rectangle((0.01, 0.01), 0.98, 0.98, fill=False)) #fill false removes background


    [axis.set_axis_off() for axis in ax.ravel()] #set all axis off
    plt.show()

def main():
    data = [10, 9, 8, 7]
    create_figure(data)

main()
