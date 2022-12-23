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
/// \brief Preselection of DPi pairs for reconstruction of B0 candidates
/// \note 
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_dplus_to_pi_k_pi;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace hf_sel_d_pi
{
DECLARE_SOA_COLUMN(PxD, pxD, float); //!
DECLARE_SOA_COLUMN(PyD, pyD, float); //!
DECLARE_SOA_COLUMN(PzD, pzD, float); //!
DECLARE_SOA_COLUMN(PxPi, pxPi, float); //!
DECLARE_SOA_COLUMN(PyPi, pyPi, float); //!
DECLARE_SOA_COLUMN(PzPi, pzPi, float); //!
} // namespace hf_sel_d_pi
DECLARE_SOA_TABLE(HfSelD, "AOD", "HFSELD", //!
                hf_sel_d_pi::PxD,
                hf_sel_d_pi::PyD,
                hf_sel_d_pi::PzD);

DECLARE_SOA_TABLE(HfSelPi, "AOD", "HFSELPI", //!
                hf_sel_d_pi::PxPi,
                hf_sel_d_pi::PyPi,
                hf_sel_d_pi::PzPi);
} // namespace o2::aod

struct HfCandidateSelectorDPi {
    // Produces AOD tables to store D tracks and Pi tracks
    Produces<aod::HfSelD> hfSelD;
    Produces<aod::HfSelPi> hfSelPi;
    
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

    void process(soa::Filtered<soa::Join<
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
            array<float, 3> pVecpiK = {track0.px() + track1.px(), track0.py() + track1.py(), track0.pz() + track1.pz()};
            array<float, 3> pVecD = {pVecpiK[0] + track2.px(), pVecpiK[1] + track2.py(), pVecpiK[2] + track2.pz()};
            
            for (const auto& trackPion : tracks) {
                // Pion selection
                if (!isPionSelected(trackPion, track0, track1, track2, indexD)) {
                    continue;
                }
                // fille the AOD tables with D and Pi track information
                hfSelD(pVecD[0], pVecD[1], pVecD[2]);
                hfSelPi(trackPion.px(), trackPion.py(),trackPion.pz());
            } // pion loop
        } // D loop
    } // process
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateSelectorDPi>(cfgc)};
}
