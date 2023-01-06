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

/// \file candidateCreatorB0.cxx
/// \brief Reconstruction of B0 candidates
/// \note Adapted from candidateCreatorXicc.cxx
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/DerivedTablesForSkimming.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_2prong;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::aod::hf_cand_b0; // from CandidateReconstructionTables.h
using namespace o2::framework::expressions;

/// Reconstruction of B0 candidates
struct HfCandidateCreatorB0 {
  Produces<aod::HfCandB0Base> rowCandidateBase; // table defined in CandidateReconstructionTables.h

  // vertexing
  Configurable<double> bz{"bz", 5., "magnetic field"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any B0 is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
  // selection
  Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};
  Configurable<int> selectionFlagD{"selectionFlagD", 1, "Selection Flag for D"};  // 1 or 3 or 7 ?
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};

  double massPi = RecoDecay::getMassPDG(kPiPlus);
  double massD = RecoDecay::getMassPDG(pdg::Code::kDMinus);
  double massDPi = 0.;

  int iterator=0;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagD);

  OutputObj<TH1F> hMassDToPiKPi{TH1F("hMassB0ToPiKPi", "D^{#minus} candidates;inv. mass (p^{#minus} K^{#plus} #pi^{#minus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  OutputObj<TH1F> hPtD{TH1F("hPtD", "D^{#minus} candidates;D^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hPtPion{TH1F("hPtPion", "#pi^{#plus} candidates;#pi^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPAD{TH1F("hCPAD", "D^{#minus} candidates;D^{#minus} cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassB0ToDPi{TH1F("hMassB0ToDPi", "2-prong candidates;inv. mass (B^{0} #rightarrow D^{#minus}#pi^{#plus} #rightarrow #pi^{#minus}K^{#plus}#pi^{#minus}#pi^{#plus}) (GeV/#it{c}^{2});entries", 500, 3., 8.)};
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx. position (cm^{2});entries", 100, 0., 1.e-4)};
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx. position (cm^{2});entries", 100, 0., 0.2)};

  // process function using preselected D Pi candidates stored in AO2D tables
  void process(aod::HfPvRefit const& pvRefits,
                aod::HfTrack0 const& track0s,
                aod::HfTrack1 const& track1s,
                aod::HfTrack2 const& track2s,
                aod::HfTrack3 const& trackPions) //, aod::HfSelD const& candDs, aod::HfSelPi const& candPis)
  {
    LOG(info) << "Process function of B0 candidate creator";
    for (const auto& pvRefit : pvRefits) {
      o2::dataformats::VertexBase primaryVertex = hf_pv_refit::getPrimaryVertex(pvRefit);
    }
    //for (const auto& [pvRefit, track0, track1, track2, trackPion] : combinations(o2::soa::CombinationsFullIndexPolicy(pvRefits, track0s, track1s, track2s, trackPions))) {
    for (const auto& track0 : track0s) {
      auto ptProng0 = sqrt(track0.px()*track0.px() + track0.py()*track0.py());
      hPtPion->Fill(ptProng0);

      o2::track::TrackParametrizationWithError<float> trackParCov0 = getTrackParCov(track0);
      LOG(info) << "trackParCov0" << trackParCov0.getX() ;
      //iterator++;
      //LOG(info) << iterator;
      //LOG(info) << track0.prong0() << track1.prong1() << track2.prong2() << trackPion.prong3();
    }
    //o2::track::TrackParametrizationWithError<TrackPrecision> getTrackParCov(const T& track)
  }


  void processFull(aod::Collision const&,
               soa::Filtered<soa::Join<
                 aod::HfCand3Prong,
                 aod::HfSelDplusToPiKPi>> const& candDs,
               aod::BigTracks const& tracks)
  {
    // Initialise fitter for B vertex (2-prong vertex filter)
    o2::vertexing::DCAFitterN<2> df2;
    df2.setBz(bz);
    df2.setPropagateToPCA(propagateToPCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(true);

    // Initial fitter to redo D-vertex to get extrapolated daughter tracks (3-prong vertex filter)
    o2::vertexing::DCAFitterN<3> df3;
    df3.setBz(bz);
    df3.setPropagateToPCA(propagateToPCA);
    df3.setMaxR(maxR);
    df3.setMaxDZIni(maxDZIni);
    df3.setMinParamChange(minParamChange);
    df3.setMinRelChi2Change(minRelChi2Change);
    df3.setUseAbsDCA(true);

    // loop over D candidates
    for (const auto& candD : candDs) {
      if (!TESTBIT(candD.hfflag(), hf_cand_3prong::DecayType::DplusToPiKPi)) { // FIXME: useless, this condition gives a flag=0 in DPlusToPiKPi selector
        continue;
      }
      if (candD.isSelDplusToPiKPi() >= selectionFlagD) {  // FIXME: useless as candidates are already filtered
        hMassDToPiKPi->Fill(invMassDplusToPiKPi(candD), candD.pt());
      }
      hPtD->Fill(candD.pt());
      hCPAD->Fill(candD.cpa());

      // track0 <-> pi, track1 <-> K, track2 <-> pi
      auto track0 = candD.prong0_as<aod::BigTracks>();
      auto track1 = candD.prong1_as<aod::BigTracks>();
      auto track2 = candD.prong2_as<aod::BigTracks>();
      auto trackParCov0 = getTrackParCov(track0);
      auto trackParCov1 = getTrackParCov(track1);
      auto trackParCov2 = getTrackParCov(track2);
      auto collision = track0.collision();

      // reconstruct 3-prong secondary vertex (D±)
      if (df3.process(trackParCov0, trackParCov1, trackParCov2) == 0) {
        continue;
      }

      const auto& secondaryVertex = df3.getPCACandidate();
      trackParCov0.propagateTo(secondaryVertex[0], bz);
      trackParCov1.propagateTo(secondaryVertex[0], bz);
      trackParCov2.propagateTo(secondaryVertex[0], bz);

      // D∓ → π∓ K± π∓
      array<float, 3> pVecpiK = {track0.px() + track1.px(), track0.py() + track1.py(), track0.pz() + track1.pz()};
      array<float, 3> pVecD = {pVecpiK[0] + track2.px(), pVecpiK[1] + track2.py(), pVecpiK[2] + track2.pz()};
      auto trackParCovPiK = o2::dataformats::V0(df3.getPCACandidatePos(), pVecpiK, df3.calcPCACovMatrixFlat(),
                                                trackParCov0, trackParCov1, {0, 0}, {0, 0});
      auto trackParCovD = o2::dataformats::V0(df3.getPCACandidatePos(), pVecD, df3.calcPCACovMatrixFlat(),
                                              trackParCovPiK, trackParCov2, {0, 0}, {0, 0});

      int index0D = track0.globalIndex();
      int index1D = track1.globalIndex();
      int index2D = track2.globalIndex();
      // int charge = track0.sign() + track1.sign() + track2.sign();

      // loop over pions
      for (const auto& trackPion : tracks) {
        // minimum pT selection
        if (trackPion.pt() < ptPionMin) {
          continue;
        }
        // reject pions that are D daughters
        if (trackPion.globalIndex() == index0D || trackPion.globalIndex() == index1D || trackPion.globalIndex() == index2D) {
          continue;
        }
        // reject pi and D with same sign
        if (trackPion.sign()*track0.sign() > 0) {
          continue;
        }

        hPtPion->Fill(trackPion.pt());
        array<float, 3> pVecPion= {trackPion.px(), trackPion.py(), trackPion.pz()};;
        auto trackParCovPi = getTrackParCov(trackPion);

        // ---------------------------------
        // reconstruct the 2-prong B0 vertex
        if (df2.process(trackParCovD, trackParCovPi) == 0) {
          continue;
        }

        // calculate relevant properties
        const auto& secondaryVertexB0 = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();

        // must be called before getTrack query
        df2.propagateTracksToVertex();
        // track.getPxPyPzGlo(pVec) modifies pVec of track
        df2.getTrack(0).getPxPyPzGlo(pVecD);
        df2.getTrack(1).getPxPyPzGlo(pVecPion);

        auto primaryVertex = getPrimaryVertex(collision);
        auto covMatrixPV = primaryVertex.getCov();
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackParCovD.propagateToDCA(primaryVertex, bz, &impactParameter0);
        trackParCovPi.propagateToDCA(primaryVertex, bz, &impactParameter1);

        hCovSVXX->Fill(covMatrixPCA[0]);
        hCovPVXX->Fill(covMatrixPV[0]);

        // get uncertainty of the decay length
        double phi, theta;
        // getPointDirection modifies phi and theta
        getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexB0, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        int hfFlag = BIT(hf_cand_b0::DecayType::B0ToDPi);

        // fill the candidate table for the B0 here:
        rowCandidateBase(collision.globalIndex(),
                          collision.posX(), collision.posY(), collision.posZ(),
                          secondaryVertexB0[0], secondaryVertexB0[1], secondaryVertexB0[2],
                          errorDecayLength, errorDecayLengthXY,
                          chi2PCA,
                          pVecD[0], pVecD[1], pVecD[2],
                          pVecPion[0], pVecPion[1], pVecPion[2],
                          impactParameter0.getY(), impactParameter1.getY(),
                          std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                          candD.globalIndex(), trackPion.globalIndex(),
                          hfFlag);

        // calculate invariant mass
        auto arrayMomenta = array{pVecD, pVecPion};
        massDPi = RecoDecay::m(std::move(arrayMomenta), array{massD, massPi});
        if (candD.isSelDplusToPiKPi() > 0) {
          hMassB0ToDPi->Fill(massDPi);
        }
      } // pi loop
    } // D loop
  } // process
  PROCESS_SWITCH(HfCandidateCreatorB0, processFull, "Process full", false);
}; // struct

/// Extends the base table with expression columns.
struct HfCandidateCreatorB0Expressions {
  Spawns<aod::HfCandB0Ext> rowCandidateB0;

  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HfCandidateCreatorB0Mc {
  Produces<aod::HfCandB0McRec> rowMcMatchRec; // table defined in CandidateReconstructionTables.h
  Produces<aod::HfCandB0McGen> rowMcMatchGen; // table defined in CandidateReconstructionTables.h

  void processMc(aod::HfCandB0 const& candidates,
                 aod::HfCand3Prong const&,
                 aod::BigTracksMC const& tracks,
                 aod::McParticles const& particlesMC)
  {
    int indexRec = -1;
    int8_t signB0 = 0;
    int8_t signD = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;

    // Match reconstructed candidates.
    for (const auto& candidate : candidates) {
      // Printf("New rec. candidate");
      flag = 0;
      origin = 0;
      debug = 0;
      auto candD = candidate.prong0();
      auto arrayDaughters = array{candD.prong0_as<aod::BigTracksMC>(),
                                  candD.prong1_as<aod::BigTracksMC>(),
                                  candD.prong2_as<aod::BigTracksMC>(),
                                  candidate.prong1_as<aod::BigTracksMC>()};
      auto arrayDaughtersD = array{candD.prong0_as<aod::BigTracksMC>(),
                                   candD.prong1_as<aod::BigTracksMC>(),
                                   candD.prong2_as<aod::BigTracksMC>()};
      // B0 → D- π+ → (π- K+ π-) π+
      // Printf("Checking B0 → D- π+");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kB0, array{-kPiPlus, +kKPlus, -kPiPlus, +kPiPlus}, true, &signB0, 2);
      if (indexRec > -1) {
        // D- → π- K+ π-
        // Printf("Checking D- → π- K+ π-");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersD, pdg::Code::kDMinus, array{-kPiPlus, +kKPlus, -kPiPlus}, true, &signD, 1);
        if (indexRec > -1) {
          flag = signB0 * BIT(hf_cand_b0::DecayType::B0ToDPi);
        } else {
          debug = 1;
          LOGF(info, "WARNING: B0 in decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }
      }
      rowMcMatchRec(flag, origin, debug);
    }

    // Match generated particles.
    for (const auto& particle : particlesMC) {
      // Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      // B0 → D- π+
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kB0, array{-int(pdg::Code::kDPlus), +kPiPlus}, true, &signB0)) {
        // Match D- -> π- K+ π-
        auto candDMC = particlesMC.rawIteratorAt(particle.daughtersIds().front());
        // Printf("Checking D- -> π- K+ π-");
        if (RecoDecay::isMatchedMCGen(particlesMC, candDMC, -int(pdg::Code::kDPlus), array{-kPiPlus, +kKPlus, -kPiPlus}, true, &signD)) {
          flag = signB0 * BIT(hf_cand_b0::DecayType::B0ToDPi);
        }
      }
      rowMcMatchGen(flag, origin);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorB0Mc, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfCandidateCreatorB0>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorB0Expressions>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorB0Mc>(cfgc)};
  return workflow;
}