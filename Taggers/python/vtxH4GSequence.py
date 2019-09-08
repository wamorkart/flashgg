import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggTkVtxMap_cfi import flashggVertexMapUnique

vtxH4GSequence = cms.Sequence(flashggVertexMapUnique)
