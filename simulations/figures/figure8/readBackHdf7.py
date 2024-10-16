import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

def makeHeatmap( fig, ax, data, dt, label, idx, maxval ):
    cbarlabels = {'A': '[Ca](uM)', 'C':'[Ca4.CaM](nM)', 'E':'[MAPK_P](uM)'}
    maxval = max( np.max( data ), maxval )
    im = ax.imshow(data, cmap='hot', interpolation='nearest', aspect='auto', clim = (0, maxval) )
    if idx == 1:
        cbar = fig.colorbar( im, ax = ax )
        cbar.ax.tick_params(labelsize=14)
        cbar.set_label(cbarlabels[label], fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    #im.tick_params(labelsize=14)
    # Add labels, title and ticks
    #ax.set_xticks(np.arange(len(cadf[0])))
    #ax.set_yticks(np.arange(len(cadf[0])))
    #ax.set_xticklabels(df.columns)
    # Convert row indices to seconds for y-axis labels
    num_labels = 5  # Adjust the number of labels as needed
    y_tick_indices = np.linspace(0, len(data) - 1, num_labels, dtype=int)
    y_labels = [i * dt for i in y_tick_indices]
    ax.set_yticks(y_tick_indices)
    ax.set_yticklabels(y_labels, fontsize = 14)
    #ax.set_title('Heatmap of DataFrame')
    ax.set_ylabel('Time (seconds)', fontsize = 16)
    if label == "A":
        ax.set_xlabel('Position (Synapse #)', fontsize = 16)
    else:
        ax.set_xlabel('Position (microns)', fontsize = 16)
    newlabel = chr( ord( label ) + idx )
    ax.text( -0.10, 1.05, newlabel, fontsize = 18, weight = "bold", transform=ax.transAxes )


def main():
    parser = argparse.ArgumentParser( description = "Display run from allNonlin?.py stored in specified hd5f file" )
    parser.add_argument( "file1", type = str, help = "Name of HDF5 file1", )
    parser.add_argument( "file2", type = str, help = "Name of HDF5 file2", )
    args = parser.parse_args()

    fig = plt.figure( figsize = (16, 20))
    gs = fig.add_gridspec( 4, 2, width_ratios = [ 1, 1.2] )
    maxes = [-100, -100, -100, -100]
    maxes = loadDisplay( args.file1, fig, gs, 0, maxes )
    maxes = loadDisplay( args.file2, fig, gs, 1, maxes )

    plt.tight_layout()
    plt.show()

def loadDisplay( filename, fig, gs, idx, maxes ):
    # Read from HDF5
    with pd.HDFStore(filename, mode='r') as store:
        # Read the DataFrame
        spinedf = store['spineCa'].values * 1000
        camdf = store['dendCaM'].values * 1e6
        mapkdf = store['dendMAPK'].values * 1000
        Vm = store['somaVm'].values * 1000
        branchVm = store['branchVm'].values * 1000
        metadata = store.get_storer('dendCaM').attrs.metadata
    
    # Print the results
    #print("DataFrame:\n", camdf)
    #print("\nNumPy array:\n", Vm)
    #print("\nMetadata:\n", metadata)

    makeHeatmap( fig, fig.add_subplot( gs[0,idx] ), spinedf,metadata['chemDt'], "A", idx, maxes[0] )
    makeHeatmap( fig, fig.add_subplot( gs[1,idx] ), camdf, metadata['chemDt'], "C", idx, maxes[1] )
    makeHeatmap( fig, fig.add_subplot( gs[2,idx] ), mapkdf, metadata['chemDt'], "E", idx, maxes[2] )

    ax3 = fig.add_subplot( gs[3,idx] )
    t = np.linspace( 0, metadata['runtime'], len( Vm ) )
    ax3.plot( t, Vm, label = "soma" )
    ax3.plot( t, branchVm + 20, label = "dendrite" )
    ax3.tick_params(axis='both', which='major', labelsize=14)
    ax3.set_xlabel('Time (seconds)', fontsize = 16)
    ax3.set_ylabel('Memb Potl (mV)', fontsize = 16)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.legend( fontsize = 14, frameon = False )
    ax3.text( -0.08, 1.05, chr( ord("G") + idx), fontsize = 18, weight = "bold", transform=ax3.transAxes )
    ax3.set_ylim( -75, max( np.max( branchVm ), maxes[3] ) + 20 )
    #print( max( np.max( branchVm ), maxes[3] ), np.max( branchVm ), maxes[3] )
    return [ np.max( spinedf ), np.max( camdf ), np.max( mapkdf ), np.max( branchVm) ]


if __name__ == '__main__':
    main()
