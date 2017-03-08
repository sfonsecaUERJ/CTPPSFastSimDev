## Instructions:

```
cmsrel CMSSW_9_0_0_pre4
git clone https://github.com/sfonsecaUERJ/CTPPSFastSimDev.git CMSSW_9_0_0_pre3/src
cd CMSSW_9_0_0_pre3/src
cmsenv
scram b -j8
```
#### Some examples using CMSDriver commands:

## MinBias samples including CTPPS FastSim

```
cmsDriver.py MinBias_13TeV_pythia8_TuneCUETP8M1_cfi -n 1000 --fast --conditions auto:run2_mc --magField 38T_PostLS1 --beamspot Realistic25ns13TeV2016Collision -s GEN,SIM,RECOBEFMIX --eventcontent FASTPU --datatier GEN-SIM-RECO --era Run2_25ns --fileout minbias.root --customise SimTransport/CTPPSHectorProducer/FastSimWithCTPPS_cff.customise

```

## NO PU + CTPPS FastSim

```
cmsDriver.py GluGluTo2Jets_M_100_7TeV_exhume_cff.py --conditions auto:run2_mc --fast -n 100 --eventcontent AODSIM -s GEN,SIM,RECOBEFMIX,DIGI:pdigi_valid,RECO --datatier GEN-SIM-DIGI-RECO --beamspot Realistic25ns13TeV2016Collision --pileup NoPileUp --era Run2_25ns --customise SimTransport/CTPPSHectorProducer/FastSimWithCTPPS_cff.customise --no_exec
 
```
## PU + CTPPS FastSim

```
cmsDriver.py GluGluTo2Jets_M_100_7TeV_exhume_cff.py --conditions auto:run2_mc --fast -n 100 --eventcontent AODSIM -s GEN,SIM,RECOBEFMIX,DIGI:pdigi_valid,RECO --datatier AODSIM --beamspot Realistic25ns13TeV2016Collision  --era Run2_25ns --pileup_input file:minbias.root --pileup AVE_35_BX_25ns --customise SimTransport/CTPPSHectorProducer/FastSimWithCTPPS_cff.customise --no_exec ```


#### You should modify the ExHume fragment (GluGluTo2Jets_M_100_7TeV_exhume_cff.py) - ECM (13TeV) and Mass range (300-2000 GeV). You can find parameters below: 
```
MassRangeHigh = cms.double(2000.0)
MassRangeLow = cms.double(300.0)
comEnergy = cms.double(13000.0)
```
