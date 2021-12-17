BuToKee_v3_2021Nov23.txt:
- Uses /BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-Custom_RK_BParking_for_RK_102X_upgrade2018_realistic_v15-v2/MINIAODSIM
- Uses BToKee.filterBySelection = cms.bool(False), so no pre-sel in BToKee builder
- Uses bainbrid/BParkingNANO:run3_trigger_2021Nov23 (?)

BuToKee_v3_2021Dec01.txt:
- Uses /BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-Custom_RK_BParking_for_RK_102X_upgrade2018_realistic_v15-v2/MINIAODSIM
- NANOAOD pre-selection opened up as much as possible !!!
- Including drForCleaning and dzForCleaning set to -1, so isPFoverlap is always false?
- All electrons (PFPF + PFLP) 
- Otherwise unknown, now obsolete
- Uses bainbrid/BParkingNANO:run3_trigger_open

BuToKee_v3_2021Dec16.txt:
- Uses /BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-Custom_RK_BParking_for_RK_102X_upgrade2018_realistic_v15-v2/MINIAODSIM
- Using 176 out of 197 files
- NANOAOD pre-selection opened up as much as possible !!!
- Including drForCleaning and dzForCleaning set to -1, so isPFoverlap is always false?
- Only PF electrons, so only PFPF category
- Current latest-and-greatest
- Uses bainbrid/BParkingNANO:run3_trigger_open_onlyPFPF

TO DO:

- Repeat above with default NANOAOD pre-selection?
