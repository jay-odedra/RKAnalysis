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
- Using 176 out of 197 files, 4971934 events out of 5690878
- NANOAOD pre-selection opened up as much as possible !!!
- Including drForCleaning and dzForCleaning set to -1, so isPFoverlap is always false?
- Only PF electrons, so only PFPF category
- Uses bainbrid/BParkingNANO:run3_trigger_open_onlyPFPF

BuToKee_v3_2021Dec20.txt:
- Uses /BuToKee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-Custom_RK_BParking_for_RK_102X_upgrade2018_realistic_v15-v2/MINIAODSIM
- Using 154 out of 197 files, 3181200 events out of 5690878
- NANOAOD pre-selection applied as default
- Only PF electrons, so only PFPF category
- Uses bainbrid/BParkingNANO:run3_trigger_onlyPFPF
- Current latest-and-greatest

2022Sep05 production:
- Used 12_4_X version of BParkingNANO ...
- Branch: DiElectronX/BParkingNANO:main (will tag it)
- B->Kee MC samples for rare, J/psi, Psi(2S) modes
- Run2022C and Run2022D(v2) data sets
- Note Run2022D is infact v2, v1 was not processed
- Used no JSON file, but appropriate GTs
