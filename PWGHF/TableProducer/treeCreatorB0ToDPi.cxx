// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file treeCreatorB0ToDPi.cxx
/// \brief Writer of the 2 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
/// \note Extended from treeCreatorBplusToD0Pi.cxx
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_b0;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
// DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);
DECLARE_SOA_COLUMN(ImpactParameter0, impactParameter0, float);
DECLARE_SOA_COLUMN(ImpactParameter1, impactParameter1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(Chi2PCA, chi2PCA, float);
DECLARE_SOA_COLUMN(NSigmaTOFBachPi, nSigmaTOFBachPi, float);
DECLARE_SOA_COLUMN(NSigmaTOFBachKa, nSigmaTOFBachKa, float);
DECLARE_SOA_COLUMN(NSigmaTPCBachPi, nSigmaTPCBachPi, float);
DECLARE_SOA_COLUMN(NSigmaTPCBachKa, nSigmaTPCBachKa, float);
DECLARE_SOA_COLUMN(MCflag, mcflag, int8_t);
// D (Prong0) selection variable
DECLARE_SOA_COLUMN(DM, dM, float);
DECLARE_SOA_COLUMN(DCt, dCt, float);
DECLARE_SOA_COLUMN(DPtProng0, dptProng0, float);
DECLARE_SOA_COLUMN(DPtProng1, dptProng1, float);
DECLARE_SOA_COLUMN(DPtProng2, dptProng2, float);
DECLARE_SOA_COLUMN(DY, dY, float);
DECLARE_SOA_COLUMN(DEta, dEta, float);
DECLARE_SOA_COLUMN(DCPA, dCPA, float);
DECLARE_SOA_COLUMN(DCPAXY, dCPAXY, float);
DECLARE_SOA_COLUMN(DChi2PCA, dChi2PCA, float);
DECLARE_SOA_COLUMN(DDecayLength, dDecayLength, float);
DECLARE_SOA_COLUMN(DDecayLengthXY, dDecayLengthXY, float);
DECLARE_SOA_COLUMN(DDecayLengthNormalised, dDecayLengthNormalised, float);
DECLARE_SOA_COLUMN(DDecayLengthXYNormalised, ddecayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(DImpactParameterProduct, dimpactParameterProduct, float);
DECLARE_SOA_COLUMN(DImpactParameter0, dimpactParameter0, float);
DECLARE_SOA_COLUMN(DImpactParameter1, dimpactParameter1, float);
DECLARE_SOA_COLUMN(DImpactParameter2, dimpactParameter2, float);
DECLARE_SOA_COLUMN(DImpactParameterNormalised0, dimpactParameterNormalised0, float);
DECLARE_SOA_COLUMN(DImpactParameterNormalised1, dimpactParameterNormalised1, float);
DECLARE_SOA_COLUMN(DImpactParameterNormalised2, dimpactParameterNormalised2, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk0Ka, nSigmaTOFTrk0Ka, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk0Pi, nSigmaTOFTrk0Pi, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk0Ka, nSigmaTPCTrk0Ka, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk0Pi, nSigmaTPCTrk0Pi, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk1Ka, nSigmaTOFTrk1Ka, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk1Pi, nSigmaTOFTrk1Pi, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk1Ka, nSigmaTPCTrk1Ka, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk1Pi, nSigmaTPCTrk1Pi, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk2Ka, nSigmaTOFTrk2Ka, float);
DECLARE_SOA_COLUMN(NSigmaTOFTrk2Pi, nSigmaTOFTrk2Pi, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk2Ka, nSigmaTPCTrk2Ka, float);
DECLARE_SOA_COLUMN(NSigmaTPCTrk2Pi, nSigmaTPCTrk2Pi, float);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

// put the arguments into the table
DECLARE_SOA_TABLE(HfCandB0Full, "AOD", "HFCANDBPFull",
                  full::RSecondaryVertex,
                  full::PtProng0,
                  full::PProng0,
                  full::PtProng1,
                  full::PProng1,
                  // full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::P,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::CPA,
                  full::CPAXY,
                  full::ImpactParameterProduct,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  full::ImpactParameterNormalised0,
                  full::ImpactParameterNormalised1,
                  hf_cand::PxProng0,
                  hf_cand::PyProng0,
                  hf_cand::PzProng0,
                  hf_cand::PxProng1,
                  hf_cand::PyProng1,
                  hf_cand::PzProng1,
                  hf_cand::Chi2PCA,
                  full::NSigmaTOFBachPi,
                  full::NSigmaTOFBachKa,
                  full::NSigmaTPCBachPi,
                  full::NSigmaTPCBachKa,
                  full::MCflag,
                  full::DM,
                  full::DPtProng0,
                  full::DPtProng1,
                  full::DPtProng2,
                  full::DY,
                  full::DEta,
                  full::DCPA,
                  full::DCPAXY,
                  full::DChi2PCA,
                  full::DDecayLength,
                  full::DDecayLengthXY,
                  full::DDecayLengthNormalised,
                  full::DDecayLengthXYNormalised,
                  full::DImpactParameterProduct,
                  full::DImpactParameter0,
                  full::DImpactParameter1,
                  full::DImpactParameter2,
                  full::DImpactParameterNormalised0,
                  full::DImpactParameterNormalised1,
                  full::DImpactParameterNormalised2,
                  full::NSigmaTOFTrk0Pi, // FIXME : change naming convention
                  full::NSigmaTOFTrk0Ka,
                  full::NSigmaTPCTrk0Pi,
                  full::NSigmaTPCTrk0Ka,
                  full::NSigmaTOFTrk1Pi,
                  full::NSigmaTOFTrk1Ka,
                  full::NSigmaTPCTrk1Pi,
                  full::NSigmaTPCTrk1Ka,
                  full::NSigmaTOFTrk2Pi,
                  full::NSigmaTOFTrk2Ka,
                  full::NSigmaTPCTrk2Pi,
                  full::NSigmaTPCTrk2Ka);

DECLARE_SOA_TABLE(HfCandB0FullEvents, "AOD", "HFCANDBPFullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandB0FullParticles, "AOD", "HFCANDBPFullP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorB0ToDPi {
  Produces<o2::aod::HfCandB0Full> rowCandidateFull;
  Produces<o2::aod::HfCandB0FullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCandB0FullParticles> rowCandidateFullParticles;

  Configurable<int> isSignal{"isSignal", 1, "save only MC matched candidates"};

  void init(InitContext const&)
  {
  }

  void process(aod::Collisions const& collisions,
               aod::McCollisions const& mccollisions,
               soa::Join<aod::HfCandB0, aod::HfCandB0McRec, aod::HfSelB0ToDPi> const& candidates,
               soa::Join<aod::McParticles_000, aod::HfCandB0McGen> const& particles,
               aod::BigTracksPID const& tracks,
               aod::HfCand2Prong const&)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto& collision : collisions) {
      rowCandidateFullEvents(
        collision.bcId(),
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        1);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto& candidate : candidates) {
      auto fillTable = [&](int CandFlag,
                           // int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY) {
        auto dCand = candidate.prong0();
        auto piCand = candidate.prong1_as<aod::BigTracksPID>();
        // adding D daughters to the table
        auto dDaughter0 = dCand.prong0_as<aod::BigTracksPID>();
        auto dDaughter1 = dCand.prong1_as<aod::BigTracksPID>();
        auto dDaughter2 = dCand.prong2_as<aod::BigTracksPID>();

        auto invMassD = 0.;
        //FIXME : adapt this to Dplus
        if (piCand.sign() > 0) {
          invMassD0 = o2::aod::hf_cand_2prong::invMassD0barToKPi(d0Cand);
        } else if (piCand.sign() < 0) {
          invMassD0 = o2::aod::hf_cand_2prong::invMassD0ToPiK(d0Cand);
        }

        // if (FunctionSelection >= 1) {
        if (std::abs(candidate.flagMcMatchRec()) >= isSignal) {

          rowCandidateFull(
            candidate.rSecondaryVertex(),
            candidate.ptProng0(),
            RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
            candidate.ptProng1(),
            RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
            // 1 << CandFlag,
            FunctionInvMass,
            candidate.pt(),
            candidate.p(),
            FunctionCt,
            candidate.eta(),
            candidate.phi(),
            FunctionY,
            candidate.decayLength(),
            candidate.decayLengthXY(),
            candidate.decayLengthNormalised(),
            candidate.decayLengthXYNormalised(),
            candidate.cpa(),
            candidate.cpaXY(),
            candidate.impactParameterProduct(),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.impactParameterNormalised0(),
            candidate.impactParameterNormalised1(),
            candidate.pxProng0(),
            candidate.pyProng0(),
            candidate.pzProng0(),
            candidate.pxProng1(),
            candidate.pyProng1(),
            candidate.pzProng1(),
            candidate.chi2PCA(),
            piCand.tofNSigmaPi(),
            piCand.tofNSigmaKa(),
            piCand.tpcNSigmaPi(),
            piCand.tpcNSigmaKa(),
            candidate.flagMcMatchRec(),
            invMassD,
            dCand.ptProng0(),
            dCand.ptProng1(),
            dCand.ptProng2(),
            o2::aod::hf_cand_2prong::yD0(d0Cand), // FIXME : to adapt to Dplus
            dCand.eta(),
            dCand.cpa(),
            dCand.cpaXY(),
            dCand.chi2PCA(),
            dCand.decayLength(),
            dCand.decayLengthXY(),
            dCand.decayLengthNormalised(),
            dCand.decayLengthXYNormalised(),
            dCand.impactParameterProduct(),
            dCand.impactParameter0(),
            dCand.impactParameter1(),
            dCand.impactParameter2(),
            dCand.impactParameterNormalised0(),
            dCand.impactParameterNormalised1(),
            dCand.impactParameterNormalised2(),
            dDaughter0.tofNSigmaPi(),
            dDaughter0.tofNSigmaKa(),
            dDaughter0.tpcNSigmaPi(),
            dDaughter0.tpcNSigmaKa(),
            dDaughter1.tofNSigmaPi(),
            dDaughter1.tofNSigmaKa(),
            dDaughter1.tpcNSigmaPi(),
            dDaughter1.tpcNSigmaKa(),
            dDaughter2.tofNSigmaPi(),
            dDaughter2.tofNSigmaKa(),
            dDaughter2.tpcNSigmaPi(),
            dDaughter2.tpcNSigmaKa());
        }
      };

      // fillTable(0, candidate.isSelB0ToDPi(), invMassB0ToDPi(candidate), ctB0(candidate), yB0(candidate));
      fillTable(0, invMassB0ToDPi(candidate), ctB0(candidate), yB0(candidate));
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (auto& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::B0ToDPi) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode())),
          particle.flagMcMatchGen());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorB0ToDPi>(cfgc));
  return workflow;
}
