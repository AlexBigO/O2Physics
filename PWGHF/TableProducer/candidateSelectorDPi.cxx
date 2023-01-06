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

/// \file candidateSelectorDPi.cxx
/// \brief Preselection of (PiKPi)Pi pairs for reconstruction of B0 candidates
/// \note
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
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::aod::hf_pv_refit;
using namespace o2::aod::hf_track_par_cov;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_dplus_to_pi_k_pi;
using namespace o2::framework::expressions;

struct HfCandidateSelectorDPi {
    // Produces AOD tables to store track information
    Produces<aod::HfPvRefit> hfPvRefit;
    Produces<aod::HfTrack0> hfTrack0;
    Produces<aod::HfTrack1> hfTrack1;
    Produces<aod::HfTrack2> hfTrack2;
    Produces<aod::HfTrack3> hfTrackPion;

    Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};
    Configurable<int> selectionFlagD{"selectionFlagD", 1, "Selection Flag for D"}; // 1 or 3 or 7?
    Filter filterDCandidates = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagD);

    /// Pion selection (D Pi <-- B0)
    /// \param trackPion is a track with the pion hypothesis
    /// \param track0 is prong0 of selected D candidate
    /// \param track1 is prong1 of selected D candidate
    /// \param track2 is prong2 of selected D candidate
    /// \param indexD is an array of globalIndex of each D prong
    /// \return true if trackPion passes all cuts
    template <typename T1, typename T2, typename T3>
    bool isPionSelected(const T1& trackPion,
        const T2& track0, const T2& track1, const T2& track2,
        const T3& indexD)
    {
        // minimum pT selection
        if (trackPion.pt() < ptPionMin) {
            return false;
        }
        // reject pions that are D daughters
        if (trackPion.globalIndex() == indexD[0] || trackPion.globalIndex() == indexD[1] || trackPion.globalIndex() == indexD[2]) {
            return false;
        }
        // reject pi D with same sign as D
        if (trackPion.sign()*track0.sign() > 0) {
            return false;
        }
        return true;
    }

    void process(aod::Collision const&,
                    soa::Filtered<soa::Join<
                    aod::HfCand3Prong,
                    aod::HfSelDplusToPiKPi>> const& candDs,
                aod::BigTracks const& tracks)
    {
        // loop over D candidates
        for (const auto& candD : candDs) {
            // track0 <-> pi, track1 <-> K, track2 <-> pi
            auto track0 = candD.prong0_as<aod::BigTracks>();
            auto track1 = candD.prong1_as<aod::BigTracks>();
            auto track2 = candD.prong2_as<aod::BigTracks>();
            array<int64_t,3> indexD = {track0.globalIndex(), track1.globalIndex(), track2.globalIndex()};

            // compute D momentum vector
            //array<float, 3> pVecpiK = {track0.px() + track1.px(), track0.py() + track1.py(), track0.pz() + track1.pz()};
            //array<float, 3> pVecD = {pVecpiK[0] + track2.px(), pVecpiK[1] + track2.py(), pVecpiK[2] + track2.pz()};

            for (const auto& trackPion : tracks) {
                // Pion selection
                if (!isPionSelected(trackPion, track0, track1, track2, indexD)) {
                    continue;
                }

                // fill AOD table with collision information (taken from track0)
                auto collision = track0.collision();
                hfPvRefit(collision.globalIndex(),
                            collision.posX(), collision.posY(), collision.posZ(),
                            collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());

                // fill the AOD tables with Pi,K,Pi and Pion track information
                auto trackParCov0 = getTrackParCovAttributes(track0);
                hfTrack0(trackParCov0[kX], trackParCov0[kAlpha],
                            trackParCov0[kY], trackParCov0[kZ], trackParCov0[kSnp],
                            trackParCov0[kTgl], trackParCov0[kQ2Pt],
                            trackParCov0[kSigY2], trackParCov0[kSigZY], trackParCov0[kSigZ2],
                            trackParCov0[kSigSnpY], trackParCov0[kSigSnpZ],
                            trackParCov0[kSigSnp2], trackParCov0[kSigTglY], trackParCov0[kSigTglZ],
                            trackParCov0[kSigTglSnp], trackParCov0[kSigTgl2],
                            trackParCov0[kSigQ2PtY], trackParCov0[kSigQ2PtZ], trackParCov0[kSigQ2PtSnp],
                            trackParCov0[kSigQ2PtTgl], trackParCov0[kSigQ2Pt2],
                            track0.px(), track0.py(), track0.pz(),

                            0);

                auto trackParCov1 = getTrackParCovAttributes(track1);
                hfTrack1(trackParCov1[kX], trackParCov1[kAlpha],
                            trackParCov1[kY], trackParCov1[kZ], trackParCov1[kSnp],
                            trackParCov1[kTgl], trackParCov1[kQ2Pt],
                            trackParCov1[kSigY2], trackParCov1[kSigZY], trackParCov1[kSigZ2],
                            trackParCov1[kSigSnpY], trackParCov1[kSigSnpZ],
                            trackParCov1[kSigSnp2], trackParCov1[kSigTglY], trackParCov1[kSigTglZ],
                            trackParCov1[kSigTglSnp], trackParCov1[kSigTgl2],
                            trackParCov1[kSigQ2PtY], trackParCov1[kSigQ2PtZ], trackParCov1[kSigQ2PtSnp],
                            trackParCov1[kSigQ2PtTgl], trackParCov1[kSigQ2Pt2],
                            track1.px(), track1.py(), track1.pz(),
                            1);
                
                auto trackParCov2 = getTrackParCovAttributes(track2);
                hfTrack2(trackParCov2[kX], trackParCov2[kAlpha],
                            trackParCov2[kY], trackParCov2[kZ], trackParCov2[kSnp],
                            trackParCov2[kTgl], trackParCov2[kQ2Pt],
                            trackParCov2[kSigY2], trackParCov2[kSigZY], trackParCov2[kSigZ2],
                            trackParCov2[kSigSnpY], trackParCov2[kSigSnpZ],
                            trackParCov2[kSigSnp2], trackParCov2[kSigTglY], trackParCov2[kSigTglZ],
                            trackParCov2[kSigTglSnp], trackParCov2[kSigTgl2],
                            trackParCov2[kSigQ2PtY], trackParCov2[kSigQ2PtZ], trackParCov2[kSigQ2PtSnp],
                            trackParCov2[kSigQ2PtTgl], trackParCov2[kSigQ2Pt2],
                            track2.px(), track2.py(), track2.pz(),
                            2);
                
                auto trackParCovPion = getTrackParCovAttributes(trackPion);
                hfTrackPion(trackParCovPion[kX], trackParCovPion[kAlpha],
                            trackParCovPion[kY], trackParCovPion[kZ], trackParCovPion[kSnp],
                            trackParCovPion[kTgl], trackParCovPion[kQ2Pt],
                            trackParCovPion[kSigY2], trackParCovPion[kSigZY], trackParCovPion[kSigZ2],
                            trackParCovPion[kSigSnpY], trackParCovPion[kSigSnpZ],
                            trackParCovPion[kSigSnp2], trackParCovPion[kSigTglY], trackParCovPion[kSigTglZ],
                            trackParCovPion[kSigTglSnp], trackParCovPion[kSigTgl2],
                            trackParCovPion[kSigQ2PtY], trackParCovPion[kSigQ2PtZ], trackParCovPion[kSigQ2PtSnp],
                            trackParCovPion[kSigQ2PtTgl], trackParCovPion[kSigQ2Pt2],
                            trackPion.px(), trackPion.py(), trackPion.pz(),
                            3);
            } // pion loop
        } // D loop
    } // process
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorDPi>(cfgc)};
}
