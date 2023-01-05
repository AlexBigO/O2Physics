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

/// \file DerivedTablesForSkimming.h
/// \brief Header file with definition of variables, methods, and tables used in the candidateSelectorDPi.cxx and candidateCreatorB0.cxx tasks
/// \note
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#ifndef PWGHF_DATAMODEL_DERIVED_TABLES_FOR_SKIMMING_H_
#define PWGHF_DATAMODEL_DERIVED_TABLES_FOR_SKIMMING_H_

namespace o2::aod
{

/// Build TrackParametrizationWithError instance from track attributes
/// \param x is the X of track evaluation
/// \param alpha is the track frame angle
/// \param arraypar is an array of 5 parameters: Y,Z,sin(phi),tg(lambda),q/pT
/// \param covpar is an array of covariance matrix elements
/// \return TrackParametrizationWithError object
/*
template <typename TrackPrecision = float, typename T1 = float, typename T2 = std::array<TrackPrecision, 5>, typename T3 = std::array<TrackPrecision, 15>>
o2::track::TrackParametrizationWithError<TrackPrecision> getTrackParCov(const T1& x, const T1& alpha, const T2& arraypar, const T3& covpar)
{
    return o2::track::TrackParametrizationWithError<TrackPrecision>(x, alpha, std::move(arraypar), std::move(covpar));
}
*/







// TODO : see if this is needed (there is already hf_pv_refit_...) --> then move the Id in hf_track_par_cov
namespace hf_pv_refit
{

enum HfPVRefit : int {
    kPosX = 0,
    kPosY,
    kPosZ,
    kCovXX,
    kCovXY,
    kCovYY,
    kCovXZ,
    kCovYZ,
    kCovZZ
};

// Helper functions for collision handling

/// Extracts primary vertex position and covariance matrix from a collision.
template <typename T>
std::array<float, 9> getPrimaryVertexAttributes(const T& collision)
{
  return {collision.posX(), collision.posY(), collision.posZ(),
            collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()};
}

///
template <typename T>
o2::dataformats::VertexBase getPrimaryVertex(const T& arraycollision)
{
  o2::math_utils::Point3D<float> vtxXYZ(arraycollision[kPosX], arraycollision[kPosY], arraycollision[kPosZ]);
  std::array<float, 6> vtxCov{arraycollision[kCovXX], arraycollision[kCovXY], arraycollision[kCovYY], arraycollision[kCovXZ], arraycollision[kCovYZ], arraycollision[kCovZZ]};
  return o2::dataformats::VertexBase{std::move(vtxXYZ), std::move(vtxCov)};
}

DECLARE_SOA_COLUMN(Id, globalIndex, int64_t); //!
DECLARE_SOA_COLUMN(PosX, posX, float); //!
DECLARE_SOA_COLUMN(PosY, posY, float); //!
DECLARE_SOA_COLUMN(PosZ, posZ, float); //!
DECLARE_SOA_COLUMN(CovXX, covXX, float); //!
DECLARE_SOA_COLUMN(CovXY, covXY, float); //!
DECLARE_SOA_COLUMN(CovYY, covYY, float); //!
DECLARE_SOA_COLUMN(CovXZ, covXZ, float); //!
DECLARE_SOA_COLUMN(CovYZ, covYZ, float); //!
DECLARE_SOA_COLUMN(CovZZ, covZZ, float); //!
} // namespace hf_pv_refit

DECLARE_SOA_TABLE(HfPVRefit, "AOD", "HFPVREFIT", //!
                hf_pv_refit::Id,
                hf_pv_refit::PosX,
                hf_pv_refit::PosY,
                hf_pv_refit::PosZ,
                hf_pv_refit::CovXX,
                hf_pv_refit::CovXY,
                hf_pv_refit::CovYY,
                hf_pv_refit::CovXZ,
                hf_pv_refit::CovYZ,
                hf_pv_refit::CovZZ);

// store track information in AO2D tables
namespace hf_track_par_cov
{

enum HfTrackParCov : int {
    kX = 0,
    kAlpha,
    kY,
    kZ,
    kSnp,
    kTgl,
    kQ2Pt,
    kSigY2,
    kSigZY,
    kSigZ2,
    kSigSnpY,
    kSigSnpZ,
    kSigSnp2,
    kSigTglY,
    kSigTglZ,
    kSigTglSnp,
    kSigTgl2,
    kSigQ2PtY,
    kSigQ2PtZ,
    kSigQ2PtSnp,
    kSigQ2PtTgl,
    kSigQ2Pt2
};


// Helper functions for track handling

/// Extracts track parameters and covariance matrix from a track
/// \param track is a track
/// \return array with track information
template <typename TrackPrecision = float, typename T>
std::array<TrackPrecision, 22> getTrackParCovAttributes(const T& track)
{
    return {track.x(), track.alpha(),
            track.y(), track.z(), track.snp(),
            track.tgl(), track.signed1Pt(),
            track.cYY(), track.cZY(), track.cZZ(),
            track.cSnpY(), track.cSnpZ(),
            track.cSnpSnp(), track.cTglY(), track.cTglZ(),
            track.cTglSnp(), track.cTglTgl(),
            track.c1PtY(), track.c1PtZ(), track.c1PtSnp(),
            track.c1PtTgl(), track.c1Pt21Pt2()};
}

template <typename TrackPrecision = float, typename T>
o2::track::TrackParametrizationWithError<TrackPrecision> getTrackParCov(const T& trackParCovAttributes)
{
    float x = trackParCovAttributes[kX];
    float alpha = trackParCovAttributes[kAlpha];
    std::array<float, 5> arraypar;
    std::array<float, 15> covpar;

    for (int i = kY; i <= kQ2Pt; ++i) {
        arraypar[i-kY] = trackParCovAttributes[i];
    }

    for (int i = kSigY2; i <= kSigQ2Pt2; ++i) {
        covpar[i-kSigY2] = trackParCovAttributes[i];
    }

    return o2::track::TrackParametrizationWithError<TrackPrecision>(x, alpha, std::move(arraypar), std::move(covpar));
}

DECLARE_SOA_COLUMN(X, x, float); //!
DECLARE_SOA_COLUMN(Alpha, alpha, float); //!
DECLARE_SOA_COLUMN(Y, y, float); //!
DECLARE_SOA_COLUMN(Z, z, float); //!
DECLARE_SOA_COLUMN(Snp, snp, float); //!
DECLARE_SOA_COLUMN(Tgl, tgl, float); //!
DECLARE_SOA_COLUMN(Q2Pt, signed1Pt, float); //!

DECLARE_SOA_COLUMN(SigY2, cYY, float); //!
DECLARE_SOA_COLUMN(SigZY, cZY, float); //!
DECLARE_SOA_COLUMN(SigZ2, cZZ, float); //!
DECLARE_SOA_COLUMN(SigSnpY, cSnpY, float); //!
DECLARE_SOA_COLUMN(SigSnpZ, cSnpZ, float); //!
DECLARE_SOA_COLUMN(SigSnp2, cSnpSnp, float); //!
DECLARE_SOA_COLUMN(SigTglY, cTglY, float); //!
DECLARE_SOA_COLUMN(SigTglZ, cTglZ, float); //!
DECLARE_SOA_COLUMN(SigTglSnp, cTglSnp, float); //!
DECLARE_SOA_COLUMN(SigTgl2, cTglTgl, float); //!
DECLARE_SOA_COLUMN(SigQ2PtY, c1PtY, float); //!
DECLARE_SOA_COLUMN(SigQ2PtZ, c1PtZ, float); //!
DECLARE_SOA_COLUMN(SigQ2PtSnp, c1PtSnp, float); //!
DECLARE_SOA_COLUMN(SigQ2PtTgl, c1PtTgl, float); //!
DECLARE_SOA_COLUMN(SigQ2Pt2, c1Pt21Pt2, float); //!

DECLARE_SOA_COLUMN(Px, px, float); //!
DECLARE_SOA_COLUMN(Py, py, float); //!
DECLARE_SOA_COLUMN(Pz, pz, float); //!
DECLARE_SOA_COLUMN(IsProng0, prong0, float); //!
DECLARE_SOA_COLUMN(IsProng1, prong1, float); //!
DECLARE_SOA_COLUMN(IsProng2, prong2, float); //!
DECLARE_SOA_COLUMN(IsProng3, prong3, float); //! // for the pion
} // namespace hf_track_par_cov

// general columns
#define HFTRACKPARCOV_COLUMNS  \
    hf_track_par_cov::X, \
    hf_track_par_cov::Alpha, \
    hf_track_par_cov::Y, \
    hf_track_par_cov::Z, \
    hf_track_par_cov::Snp, \
    hf_track_par_cov::Tgl, \
    hf_track_par_cov::Q2Pt, \
    hf_track_par_cov::SigY2, \
    hf_track_par_cov::SigZY, \
    hf_track_par_cov::SigZ2, \
    hf_track_par_cov::SigSnpY, \
    hf_track_par_cov::SigSnpZ, \
    hf_track_par_cov::SigSnp2, \
    hf_track_par_cov::SigTglY, \
    hf_track_par_cov::SigTglZ, \
    hf_track_par_cov::SigTglSnp, \
    hf_track_par_cov::SigTgl2, \
    hf_track_par_cov::SigQ2PtY, \
    hf_track_par_cov::SigQ2PtZ, \
    hf_track_par_cov::SigQ2PtSnp, \
    hf_track_par_cov::SigQ2PtTgl, \
    hf_track_par_cov::SigQ2Pt2, \
    hf_track_par_cov::Px, \
    hf_track_par_cov::Py, \
    hf_track_par_cov::Pz

DECLARE_SOA_TABLE(HfTrack0, "AOD", "HFTRACK0", //!
                HFTRACKPARCOV_COLUMNS,
                hf_track_par_cov::IsProng0);

DECLARE_SOA_TABLE(HfTrack1, "AOD", "HFTRACK1", //!
                HFTRACKPARCOV_COLUMNS,
                hf_track_par_cov::IsProng1);

DECLARE_SOA_TABLE(HfTrack2, "AOD", "HFTRACK2", //!
                HFTRACKPARCOV_COLUMNS,
                hf_track_par_cov::IsProng2);

DECLARE_SOA_TABLE(HfTrack3, "AOD", "HFTRACK3", //!
                HFTRACKPARCOV_COLUMNS,
                hf_track_par_cov::IsProng3);
} // namespace o2::aod

#endif // PWGHF_DATAMODEL_DERIVED_TABLES_FOR_SKIMMING_H_