MIN_M = 2
MAX_M = 9

T_POP_POST = 10

hippoChem = {
    'name' : 'hippo-chem',
    'populationSize' : 400000,
    'connectionProb' : 0.05,
    'ensembleSize' : 1000,
    'inputSpacing' : 2e-6,
    'dendriteLength' : 10e-3,
    'backgroundRate' : 1e-2,
    'inputDuration' : 2.0,
    'flankLength' : 0.0
}
    
hippoCICR = {
    'name' : 'hippo-CICR',
    'populationSize' : 400000,
    'connectionProb' : 0.05,
    'ensembleSize' : 1000,
    'inputSpacing' : 2e-6,
    'dendriteLength' : 10e-3,
    'backgroundRate' : 1e-2,
    'inputDuration' : 2e-1, 
    'flankLength' : 0.0
}

hippoElec = {
    'name' : 'hippo-elec',
    'populationSize' : 400000,
    'connectionProb' : 0.05,
    'ensembleSize' : 1000,
    'inputSpacing' : 10e-6,
    'dendriteLength' : 10e-3,
    'backgroundRate' : 1e-1,
    'inputDuration' : 4e-3,
    'flankLength' : 100e-6
}

cortexChem = {
    'name' : 'cortex-chem',
    'populationSize' : 100000,
    'connectionProb' : 0.2,
    'ensembleSize' : 1000,
    'inputSpacing' : 2e-6,
    'dendriteLength' : 10e-3,
    'backgroundRate' : 1e-1,
    'inputDuration' : 2.0,
    'flankLength' : 0.0
}

cortexCICR = {
    'name' : 'cortex-CICR',
    'populationSize' : 100000,
    'connectionProb' : 0.2,
    'ensembleSize' : 1000,
    'inputSpacing' : 2e-6,
    'dendriteLength' : 10e-3,
    'backgroundRate' : 1e-1,
    'inputDuration' : 2e-1, 
    'flankLength' : 0.0
}

cortexElec = {
    'name' : 'cortex-elec',
    'populationSize' : 100000,
    'connectionProb' : 0.2,
    'ensembleSize' : 1000,
    'inputSpacing' : 10e-6,
    'dendriteLength' : 10e-3,
    'backgroundRate' : 1.0,
    'inputDuration' : 4e-3, 
    'flankLength' : 100e-6
}

