from __future__ import print_function
import moose
import pylab
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import rdesigneur as rd
from scipy import stats
import time
import argparse
import itertools
import pandas as pd

displayMoogli = False
displayScatter = False
displayRuns = False
dumpSpatialPlotsFlag = False
scatterParamToUseInStats = 0
realStartTime = time.time()
settleTime = 2

# Stim amplitude is unitless, defined in units of A mol conc.
# Stim Width is unitless, defined as multiple of diffusion length.
# Stim Vel is unitless, defined in terms of diffusion length by  time units of diff constt.
# diffConst is defined in terms of square of diffusion length by unit time.
# diffLength here is in SI units: m^2/s
#

params = {
    'diffusionLength':1.0e-6,  # Diffusion characteristic length, used as voxel length too.
    'dendDiameter': 0.6e-6,  # Diameter of section of dendrite in model
    'dendLength': 100e-6,   # Length of section of dendrite in model
    'spineSizeScale': 1.0,  # Length scaling for spines. Vol wil be x^3.
    'diffConstCa':100.0e-12,      # Diffusion constant of Ca
    'diffConstMAPK': 5e-12, # Diffusion constant for MAPK
    'diffConstPP': 4e-12, # Diff constant for MAPK-activated phosphatase
    'diffConstCaM': 20e-12, # Diff constant for CaM
    'CaActivateRafKf': 6e6,	# 1/sec.mM: rate for activation of Raf by Ca
    'blankVoxelsAtEnd':10,  # of voxels to leave blank at end of cylinder
    'preStimTime':1.0,      # Time to run before turning on stimulus.
    'stimBurstTime':2.0,    # Time for a stimulus burst
    'postStimTime':10.0,    # Time to run after stimulus. ~3x decay time
    'runtime':20.0,         # Simulation run time
    'checkPoint':1.0,       # How often to do a checkpoint.
    'chemDt':0.0005,        # compute timestep for chemical signaling.
    'chemPlotDt':0.005,     # Plotting timestep for chemical signaling.
    'elecPlotDt':0.1e-3,    # Plotting timestep for electrical signaling.
    'spineSpacing':2.0e-6,  # Spacing between spines.
    'stimSpacing':4,        # Stimulus spacing, in terms of # of spines.
    'stimFuncDt':0.005,     # 5 ms resolution for stimulus function.
    'meanSpikeRate':0.000001,    # Basal mean rate for all synapses.
    'activeSpikeRate':20.0,  # Active input rate on specified synapses.
    'baseGabaSpikeRate':1.0,    # 1 Hz.
    'thetaGabaSpikeRate':0.5, # This is the peak, not the average theta.
    'thetaFreq':8.0, 		# 
    'amparSynapseWeight': 20.0,    #    Overriden by args
    'nmdarSynapseWeight': 20.0,    #    Overriden by args
    'gabarSynapseWeight': 20.0,    #    Not overridden
    'LCaDensity': 8.0,     # Dend Channel density for LCa channels.
    'adaptorScale':60.0e3, # Adaptor scale factor from conc to Na density in Seimens/m^2.
    'CaPsdScale': 0.16,    # Adaptor scale from psd elec Ca to chem conc.
    'Em': -70.0e-3,     #   Resting potential of neuron
    'refractoryPeriod':0.010,  # 10 ms refractory time.
    'cellModel': './cells/VHC-neuron.CNG.swc',  # Cell morphology file
    'chemModel': './NN_mapk16.g',  # Chemical model file.
    'fnumber': 0,  # Output file index
    'seed': 1234,   # Seeds random numbers
    'seqDx': 4.0e-6,    # Sequence spatial interval
    'seqDt': 3.0,       # Sequence time interval
}

# Here we define a string for each of the 5 stimulus timings in sequence:
# Comment in the seq order you want.
seq = [0, 1, 2, 3, 4]   # ordered
#seq = [4, 2, 0, 3, 1]   #scrambled
stimStrList = ["{0} + (t>{1}) * (t<{2}) * {3}".format( 0, params['preStimTime'] + params['seqDt']*seq[idx], params['preStimTime'] + params['stimBurstTime']+ params['seqDt']*seq[idx], params['activeSpikeRate'] ) for idx in range(5) ]

gabaRateExpression = "{} + 2*{} * cos(3.14159*t*{})^2".format( params['baseGabaSpikeRate'], params['thetaGabaSpikeRate'], params['thetaFreq'] )
#print( "GRE = ", gabaRateExpression )
gabaRateExpression = "0.001"

plotlist = [
            ['soma', '1', '.', 'Vm', 'Soma Vm'],
            ['head12', '1', 'Ca_conc', 'Ca', 'head2 eCa'],
            ['dend36', '1', 'dend/DEND/P_MAPK', 'conc', 'dend36_P_MAPK'],
        ]

moogList = [
    #['#', '1', '.', 'Vm', 'Memb potential'],
    ['#', '1', 'DEND/P_MAPK', 'conc', '[P_MAPK] (uM)',0, 0.15],
    ['#', '1', 'DEND/Ca', 'conc', '[Ca] (uM)',0, 0.5],
    ['#', '1', 'SPINE/Ca', 'conc', '[Ca] (uM)', 0, 0.5, True, 2],
]

moogList = [[]]

def setGluVec( idx, freq ):
    #gluVec = moose.wildcardFind( "/model/elec/head#/glu/sh/synapse/synInput_rs" )
    #nmdaVec = moose.wildcardFind( "/model/elec/head#/NMDA/sh/synapse/synInput_rs" )

    if idx >= 0 and idx < 80:
        #gluVec[idx].rate = freq
        #nmdaVec[idx].rate = freq
        moose.element( "/model/elec/head{}/glu/sh/synapse/synInput_rs".format(idx) ).rate = freq
        moose.element( "/model/elec/head{}/NMDA/sh/synapse/synInput_rs".format(idx) ).rate = freq

def zeroGluVec():
    gluVec = moose.wildcardFind( "/model/elec/head#/glu/sh/synapse/synInput_rs" )
    nmdaVec = moose.wildcardFind( "/model/elec/head#/NMDA/sh/synapse/synInput_rs" )
    for gg in gluVec:
        gg.rate = 0
    for nn in nmdaVec:
        nn.rate = 0


def stimFunc( freq, seqDt, seqOnT, seqDx, seqN, doGroup, xoffset, perm ):
    eps = 1e-5
    t = moose.element( '/clock' ).currentTime
    if abs( t - int( t ) ) < eps:
        print( ".", end = "", flush = True )
    #print( "{:.3f}\t{}".format( t, [ee.rate for ee in moose.wildcardFind("/model/elec/head#/glu/sh/synapse/synInput_rs" )] ) )
    if t < settleTime-eps:
        return
    idx = int( (eps + t - settleTime) / seqDt )
    if doGroup :
        if idx == 0:
            for ii in range(seqN ):
                setGluVec( ii * seqDx + xoffset, freq )
        elif idx == 1:
            zeroGluVec()
        else:
            return
    else:
        if idx == (seqN + 1):
            zeroGluVec()
            return
        if idx > seqN:
            return

        idx2 = int( (eps + t - settleTime - seqOnT) / seqDt )
        if perm == 0:
            xidx = seqDx * idx
        else:
            pp = itertools.permutations( list( range( seqN ) ) )
            pn = next( itertools.islice( pp, perm, None ) )
            xidx = seqDx * pn[ idx%len(pn) ]
        #print( t, idx, idx2, xidx )
        if idx2 < idx:
            setGluVec( xidx + xoffset, freq )
        else:
            setGluVec( xidx + xoffset, 0 )
        setGluVec( xidx +xoffset - seqDx, 0 )
    #print( "{:.3f}\t{}".format( t, [ee.rate for ee in moose.wildcardFind("/model/elec/head#/glu/sh/synapse/synInput_rs" )] ) )

def dumpPlotsToHDF( params, args ):
    somaPlot = moose.element( '/model/graphs/plot0' )
    branchPlot = moose.element( '/model/graphs/plot4' )
    CaMplot = {}
    mapkplot = {}
    spineplot = {}
    camDt = 0.0
    for idx, pp in enumerate( moose.wildcardFind( '/model/graphs/plot1[]')):
        CaMplot["d{}".format(idx)] = pp.vector
        camDt = pp.dt

    for idx, pp in enumerate( moose.wildcardFind( '/model/graphs/plot2[]')):
        mapkplot["d{}".format(idx)] = pp.vector

    for idx, pp in enumerate( moose.wildcardFind( '/model/graphs/plot3[]')):
        spineplot["d{}".format(idx)] = pp.vector
    #CaMplot['time'] = np.array( np.linspace( 0, args.runtime, len( caplot[0] ) ) )
    print( idx, camDt, len( CaMplot["d0"] ) )
    df = pd.DataFrame(CaMplot)
    df2 = pd.DataFrame( mapkplot )
    df3 = pd.DataFrame( spineplot )
    metadata = { **params, **vars(args) }
    metadata[ 'chemDt' ] = camDt
    metadata[ 'vDt' ] = somaPlot.dt
    with pd.HDFStore(args.outputFile, mode='w') as store:
        store.put('dendCaM', df, format='table', data_columns=True)
        store.put('dendMAPK', df2, format='table', data_columns=True)
        store.put('spineCa', df3, format='table', data_columns=True)
        store.put('somaVm', pd.Series(somaPlot.vector))  # Store as a Series
        store.put('branchVm', pd.Series(branchPlot.vector))  # Store as a Series
        store.get_storer('dendCaM').attrs.metadata = metadata

    


'''
# Sample DataFrame
data = {'col1': [1, 2], 'col2': [3, 4]}
df = pd.DataFrame(data)

# 1-D NumPy array
long_array = np.random.rand(1000)  # Adjust the size as needed

# Metadata dictionary
metadata = {'version': 1.0, 'description': 'Sample data with long array'}

# Save to HDF5
with pd.HDFStore('my_data.h5', mode='w') as store:
    store.put('my_table', df, format='table', data_columns=True)
    store.put('long_array', pd.Series(long_array))  # Store as a Series
    store.get_storer('my_table').attrs.metadata = metadata
'''


def runModel( params, args ):
    buildModel( params, args )
    t0 = time.time()
    moose.reinit()
    zeroGluVec()
    moose.start( args.runtime )
    dumpPlotsToHDF( params, args ) # This could be a utility func for nsdf
    moose.delete( '/model' )
    moose.delete( '/library' )
    print( "wallclock time = {:.3f}".format( time.time() -t0 ) )

def buildModel( params, args ):
    diffusionLength = params['diffusionLength']
    dendLength = params['dendLength']
    diffusionLength = params['diffusionLength']
    library = moose.Neutral( '/library' )

    #chanpath = os.path.dirname( os.path.realpath(__file__)) + '/chans/proto22.'
    chanpath = ''
    moose.seed( params['seed'] )
    rdes = rd.rdesigneur(
        useGssa = False,
        turnOffElec = False,
        chemDt = params['chemDt'],
        chemPlotDt = params['chemPlotDt'],
        diffusionLength = diffusionLength,
        spineProto = [['makeExcSpine()', 'spine']],
        chanProto = [
            [ chanpath + 'make_K_AHP()', 'K_AHP' ],
            [ chanpath + 'make_K_A()', 'K_A' ],
            [ chanpath + 'make_K_C()', 'K_C' ],
            [ chanpath + 'make_K_DR()', 'K_DR' ],
            [ chanpath + 'make_Na()', 'Na' ],
            [ chanpath + 'make_Ca_conc()', 'Ca_conc' ],
            [ chanpath + 'make_Ca()', 'Ca' ],
            [ chanpath + 'make_NMDA()', 'NMDA' ],
            [ chanpath + 'make_glu()', 'glu' ],
            [ chanpath + 'make_GABA()', 'GABA' ],
        ],
        chemProto = [[params['chemModel'], 'chem']],
        # ibran, name, somaDia, somaLen, dendDia, dendLen, dendNumSeg,  branchDia, branchLen, branchNumSeg
        cellProto = [['branchedCell', 'soma', 10e-6, 20e-6, 1.2e-6, 200e-6, 5, 0.6e-6, 160e-6, 16]],
        chanDistrib = [
            ["Ca_conc", "#", "tau", "0.0133" ],
            ["Ca", "#dend#,#basal#,#apical#,#branch#", "Gbar", str( params["LCaDensity"] ) ],
            ["Ca", "#soma#", "Gbar", "40" ],
            ["Na", "#dend#,#basal#", "Gbar", "60" ],
            ["Na", "#soma#", "Gbar", "600" ],
            ["Na", "#apical#,#branch#", "Gbar", "40+40*exp(-p/200e-6)" ],
            ["K_DR", "#dend#,#basal#", "Gbar", "(p < 400e-6)*200" ],
            ["K_DR", "#soma#", "Gbar", "360" ],
            ["K_DR", "#apical#,#branch#", "Gbar", "60+40*(p < 125e-6)" ],
            ["K_AHP", "#", "Gbar", "8" ],
            ["K_C", "#basal#,#dend#,#apical#,#branch#", "Gbar", "50+150*exp(-p/200e-6)" ],
            ["K_C", "#soma#", "Gbar", "100" ],
            ["K_A", "#soma#", "Gbar", "50" ],
            ["K_A", "#dend#,#apical#,#branch#", "Gbar", "50*(1 + 2.0e-6/(dia + 0.1e-6))" ],
            ["GABA", "#branch#", "Gbar", "10 + 30*(p < 125e-6)" ],
        ],
        spineDistrib = [['spine','branch1_#', str(params['spineSpacing']),'-1e-7', str( params['spineSizeScale'] ), '0.0', '0', '0' ]],
        chemDistrib = [
            ['DEND', 'branch1_#', 'dend', '1', diffusionLength ],
            ['SPINE', '#', 'spine', '1', 'DEND' ],
            ['PSD', '#', 'psd', '1', 'DEND' ]
        ],
        # Ideally should be synced. There should be a way to do this.
        stimList = [
                [ '#', str( params['gabarSynapseWeight'] ), 'GABA', 'randsyn', str(gabaRateExpression) ],
                [ 'head#',  str( args.weight ), 'glu', 'periodicsyn', "0.001"],
                [ 'head#',  str( args.weight ), 'NMDA', 'periodicsyn', "0.001"],
            ],
        adaptorList = [
            [ 'Ca_conc', 'Ca', 'PSD/Ca_input', 'concInit', 2e-6, params['CaPsdScale'] ],
            ['Ca_conc','Ca','DEND/Ca_input','concInit',2.0e-6, 0.0001],
            [ 'DEND/channel_p', 'conc', 'Na', 'modulation', 1.0, params['adaptorScale']],
        ],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Vm'],
            ['#', '1', 'DEND/CaM/CaM_Ca4', 'conc', 'Dend Ca4CaM conc'],
            ['#', '1', 'DEND/P_MAPK', 'conc', 'P_MAPK conc'],
            ['#', '1', 'SPINE/Ca', 'conc', 'Spine Ca conc'],
            ['branch1_9', '1', '.', 'Vm', 'branch Vm'],
        ],
        moogList = moogList
    )

    ############## Set Ca diffusion const ##########################
    for ca in moose.wildcardFind( '/library/##/Ca[ISA=PoolBase]' ):
        ca.diffConst = params['diffConstCa']

    ############## Set MAPK diffusion const ##########################
    temp = params['diffConstMAPK']
    moose.element( '/library/chem/kinetics/DEND/P_MAPK' ).diffConst = temp
    moose.element( '/library/chem/kinetics/DEND/MAPK' ).diffConst = temp
    ############## Set PP diffusion const ##########################
    temp = params['diffConstPP']
    moose.element( '/library/chem/kinetics/DEND/reg_phosphatase' ).diffConst = temp
    moose.element( '/library/chem/kinetics/DEND/inact_phosphatase' ).diffConst = temp
    ############## Set CaM diffusion const ##########################
    temp = params['diffConstCaM']
    moose.element( '/library/chem/kinetics/DEND/CaM/CaM_Ca4' ).diffConst = temp
    moose.element( '/library/chem/kinetics/DEND/CaM/CaM_Ca3' ).diffConst = temp
    moose.element( '/library/chem/kinetics/DEND/CaM/CaM_Ca2' ).diffConst = temp
    moose.element( '/library/chem/kinetics/DEND/CaM/CaM_Ca' ).diffConst = temp
    moose.element( '/library/chem/kinetics/DEND/CaM/CaM' ).diffConst = temp
    ############## Set resting potential ##########################
    for i in moose.wildcardFind( "/library/##[][ISA=CompartmentBase]" ):
        i.Em = params[ 'Em' ]
        i.initVm = params[ 'Em' ]
    ############## Set sensitivity to Ca ##########################
    moose.element( '/library/chem/kinetics/DEND/Ca_activate_Raf' ).Kf = 6e6

    #################### Build the model ##########################
    rdes.buildModel()
    moose.reinit()
    moose.seed( 1 )
    '''
    moose.le( '/model/elec' )
    moose.le( '/model/elec/head23' )
    moose.le( '/model/elec/head23/glu/sh/synapse' )
    moose.le( '/model/elec/head23/NMDA' )
    moose.showfield( '/model/elec/head23/glu/sh/synapse/synInput_rs' )
    moose.showmsg( '/model/elec/head23/glu/sh/synapse/synInput_rs' )
    '''
    moose.delete ( '/model/stims' )
    for ee in moose.wildcardFind( "/model/elec/#/GABA/sh/synapse/synInput_rs" ):
        ee.doPeriodic = False
        ee.rate = params['baseGabaSpikeRate']

    gluVec = moose.wildcardFind( "/model/elec/head#/glu/sh/synapse/synInput_rs" )
    nmdaVec = moose.wildcardFind( "/model/elec/head#/NMDA/sh/synapse/synInput_rs" )
    #rdes.displayMoogli( 0.01, params['runtime'], 0.0, colormap = 'plasma', mergeDisplays = True, bg = 'default' )
    moose.Neutral( '/model/stims' )
    pr = moose.PyRun( '/model/stims/stimulus' )
    pr.tick = 14
    moose.setClock( 14, params[ 'stimFuncDt' ] )
    pr.runString = 'stimFunc({}, {}, {}, {}, {}, {}, {}, {})'.format( 
            args.frequency, args.seqDt, args.seqOnT, args.seqDx, 
            args.seqN, args.groupStimulus, args.xOffset, args.permutation )


def main():
    parser = argparse.ArgumentParser( description = "Deliver patterned synaptic stims to neuron with chemical and electrical nonlinearities in its dendrites" )
    parser.add_argument( "-n", "--numProcesses", type = int, help = "Number of processes to launch, default = 1", default = 1 )
    parser.add_argument( '-spk', '--spiking', action="store_true", help ='Flag: when set, use high Na/K channel densities in soma to get spiking.' )
    parser.add_argument( '-v', '--voltage_clamp', action="store_true", help ='Flag: when set, do voltage clamp for glu and GABA currents respectively.')
    parser.add_argument( '-d', '--deterministic', action="store_true", help ='Flag: when set, use deterministic ODE solver. Normally uses GSSA stochastic solver.')
    parser.add_argument( "-c", "--chemModelName", type = str, help = "Optional: specify name of chemical model file, assumed to be in ./Models dir.", default = "NN_mapk16.g" )
    parser.add_argument( "-m", "--morphologyFile", type = str, help = "Optional: specify name of morphology file. Default uses ball-and-stick model." )
    parser.add_argument( "-f", "--frequency", type = float, help = "Optional: Frequency of spiking input to stimulated synapses. Default = 20", default = 20 )
    parser.add_argument( "--seqDt", type = float, help = "Optional: Timing between successive input patterns. Default = 3", default = 3 )
    parser.add_argument( "--seqOnT", type = float, help = "Optional: On time of input patterns. Default = 2", default = 2 )
    parser.add_argument( "--seqDx", type = int, help = "Optional: Increment in position of stimulus, in units of # synapses. Default = 1", default = 1 )
    parser.add_argument( "--seqN", type = int, help = "Optional: Number of stimuli to deliver, in units of # synapses. Default = 5", default = 5 )
    parser.add_argument( '-g', '--groupStimulus', action="store_true", help ='Flag: when set, deliver all seqN stimuli as a group, spaced at seqDx, in a single interval of seqDt. Else deliver a sequence. Default: deliver sequence.')
    parser.add_argument( "--xOffset", type = int, help = "Optional: x position of stimulus start point, in units of # synapses. Default = 5", default = 5 )
    parser.add_argument( "-p", "--permutation", type = int, help = "Optional: Index of permutation to use from set of all permutations of length seqN. Default = 0, ie, no permutation.", default = 0 )
    parser.add_argument( "-w", "--weight", type = float, help = "Optional: Weight of AMPAR and NMDNAR synapses. Default = 20", default = 20 )
    parser.add_argument( "-o", "--outputFile", type = str, help = "Optional: name hdf file in which to dump simulation results. Default= out7_<freq>_<seqDt>_<seqDx>_<seqN>_[grp/seq]_<perm>.h5" )
    parser.add_argument( '-r', "--runtime", type = float, help = "Optional: Runtime. Default = 25", default = 25 )
    args = parser.parse_args()
    if not args.outputFile:
        args.outputFile = 'out7_frq{}_Dt{}_Dx{}_N{}_{}_p{}_w{}.h5'.format( 
                args.frequency, args.seqDt, args.seqDx, 
                args.seqN, "grp" if args.groupStimulus else "seq",
                args.permutation, args.weight )
    runModel( params, args)

if __name__ == '__main__':
    main()
