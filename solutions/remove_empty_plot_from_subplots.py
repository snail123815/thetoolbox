import matplotlib.pyplot as plt


columns = ["a", "b", "c", "d"]
fig, axes = plt.subplots(2,2, sharex=True)

# To reliably get which plot is filled, better use a switch dictionary
plotted = {}
for c, ax in zip(columns, axes.ravel()):
    plotted[ax] = 0
    if c == "d":
        print("I didn't actually need 'd'")
        continue
    ax.plot([1,2,3,4,5,6,5,4,6,7])
    ax.set_title(c)
    plotted[ax] = 1

if axes.ndim == 2:
    for a, axs in enumerate(reversed(axes)):
        for b, ax in enumerate(reversed(axs)):
            if plotted[ax] == 0:
                # if a switch dict is not possible
                # ax.get_lines(), ax.get_images(), ax.findobj()
                ax.set_axis_off()
                # now find the plot above
                axes[-2-a][-1-b].xaxis.set_tick_params(which='both', labelbottom=True)
            else:
                break # usually only the last few plots are empty, but delete this line if not the case
else:
    for i, ax in enumerate(reversed(axes)):
        if plotted[ax] == 0:
            ax.set_axis_off()
            axes[-2-i].xaxis.set_tick_params(which='both', labelbottom=True)
            # should also work with horizontal subplots
            # all modifications to the tick params should happen after this
        else:
            break

plt.show()
