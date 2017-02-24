import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_1.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_10.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_11.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_12.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_13.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_14.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_16.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_17.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_18.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_19.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_2.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_20.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_21.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_3.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_4.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_5.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_6.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_7.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_8.root',
            'file:/eos/user/d/dilson/pps/root/GluGluTo2Jets_LHCT_NoPU_9.root'
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
                fileName = cms.string('outValidation_LHCT_ggTo2Jets.root')
)

process.options = cms.untracked.PSet(
            SkipEvent = cms.untracked.vstring('ProductNotFound')
            )

process.p = cms.Path(process.validation)
