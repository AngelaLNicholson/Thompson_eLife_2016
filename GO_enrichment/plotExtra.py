#Function to add ticks
#from here:http://stackoverflow.com/questions/10615960/matplotlib-add-strings-as-custom-x-ticks-but-also-keep-existing-numeric-tick
import matplotlib.pyplot as plt

def addticks(ax,newLocs,newLabels,pos='x'):
    # Draw to get ticks
    plt.draw()

    # Get existing ticks
    if pos=='x':
        locs = ax.get_xticks().tolist()
        labels=[x.get_text() for x in ax.get_xticklabels()]
    elif pos =='y':
        locs = ax.get_yticks().tolist()
        labels=[x.get_text() for x in ax.get_yticklabels()]
    else:
        print("WRONG pos. Use 'x' or 'y'")
        return

    # Build dictionary of ticks
    Dticks=dict(zip(locs,labels))
    print 'Dticks', Dticks
    # Add/Replace new ticks
    for Loc,Lab in zip(newLocs,newLabels):
        Dticks[Loc]=Lab

    # Get back tick lists
    locs=list(Dticks.keys())
    labels=list(Dticks.values())

    # Generate new ticks
    if pos=='x':
        ax.set_xticks(locs)
        ax.set_xticklabels(labels)
    elif pos =='y':
        ax.set_yticks(locs)
        ax.set_yticklabels(labels)

'''
#Get numpy arrays
x=np.linspace(0,2)
y=np.sin(4*x)

#Start figure
fig = plt.figure()
ax=fig.add_subplot(111)

#Plot Arrays
ax.plot(x,y)
#Add a twin axes
axr=ax.twinx()

#Add more ticks
addticks(ax,[1/3,0.75,1.0],['1/3','3/4','Replaced'])
addticks(axr,[0.5],['Miguel'],'y')

#Save figure
plt.savefig('MWE.pdf')
'''