import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_1.root',
    ),
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.validation = cms.EDAnalyzer('CTPPSFastValidation',
    MCEvent = cms.untracked.InputTag("LHCTransport"),
    ChgGenPartCollectionName = cms.untracked.InputTag("genParticles"),
    psimHitTag = cms.InputTag('CTPPSSimHits','CTPPSHits'),
    recHitTag = cms.InputTag("CTPPSFastRecHits","CTPPSFastRecHits"),
    tracksPPSTag = cms.InputTag("CTPPSFastTracks","CTPPSFastTrack"),
    jetsTag = cms.InputTag('ak4PFJets'),
    fPhysChannelTag = cms.string('ExHuMe (ggTo2Jets)')  
)

process.TFileService = cms.Service("TFileService",
                #fileName = cms.string('/eos/user/d/dilson/pps/outValidation_LHCT.root')
                fileName = cms.string('outValidation_LHCT.root')
)

process.options = cms.untracked.PSet(
            SkipEvent = cms.untracked.vstring('ProductNotFound')
            )

process.p = cms.Path(process.validation)
