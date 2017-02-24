import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/d/dilson/public/pps/SimPPS/CMSSW_8_1_0_pre5/src/simevent_CTPPS_100k.root'
        #'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
    )
)

process.demo = cms.EDAnalyzer('SimTransportValidation',
#process.demo = cms.EDAnalyzer('CTPPSValidation',
    #lable of MC event
    MCEvent = cms.untracked.InputTag("source"),
    #label of charged MC particles
    ChgGenPartCollectionName = cms.untracked.InputTag("chargeParticles")
)

process.TFileService = cms.Service("TFileService",
                fileName = cms.string('ctpps_validation.root')
)

process.p = cms.Path(process.demo)
