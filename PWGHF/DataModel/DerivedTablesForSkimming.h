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
/// \brief Header file with definition of variables, methods, and tables
//  used to fold (unfold) track and primary vertex information by writing (reading) AO2Ds
/// \note
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#ifndef PWGHF_DATAMODEL_DERIVED_TABLES_FOR_SKIMMING_H_
#define PWGHF_DATAMODEL_DERIVED_TABLES_FOR_SKIMMING_H_

namespace o2::aod
{

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

/// Extracts primary vertex position and covariance matrix from a collision
/// \param collision
/// \return array with primary vertex information
template <typename T>
std::array<float, 9> getPrimaryVertexAttributes(const T& collision)
{
  return {collision.posX(), collision.posY(), collision.posZ(),
            collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()};
}

/// Constructs primary vertex from its attributes
/// \param arraycollision is a std::array<float, 9> containing primary vertex attributes
/// \return VertexBase instance
template <typename T>
o2::dataformats::VertexBase getPrimaryVertex(const T& pvRefit)
{
  o2::math_utils::Point3D<float> vtxXYZ(pvRefit.posX(), pvRefit.posY(), pvRefit.posZ());
  std::array<float, 6> vtxCov{pvRefit.covXX(), pvRefit.covXY(), pvRefit.covYY(), pvRefit.covXZ(), pvRefit.covYZ(), pvRefit.covZZ()};
  return o2::dataformats::VertexBase{std::move(vtxXYZ), std::move(vtxCov)};
}

// CAREFUL: the getters names shall be the same as the ones of the getPrimaryVertex method in Common/Core/trackUtilities.h
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

DECLARE_SOA_TABLE(HfPvRefit, "AOD", "HFPVREFIT", //!
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

/// Extracts track parameters and covariance matrix from a track
/// \param track
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

/// Constructs TrackParametrizationWithError instance from trackTuple
/// \param trackuple is a row of joined tables HFTrack<0,1,2,3>
/// \return TrackParametrizationWithError object for prong0
template <typename TrackPrecision = float, typename T>
o2::track::TrackParametrizationWithError<TrackPrecision> getTrackParCov0(const T& trackTuple)
{
    std::array<TrackPrecision, 5> arrayparProng0 = {trackTuple.yProng0(), trackTuple.zProng0(), trackTuple.snpProng0(),
                                            trackTuple.tglProng0(), trackTuple.signed1PtProng0()};
    std::array<TrackPrecision, 15> covparProng0 = {trackTuple.cYYProng0(), trackTuple.cZYProng0(), trackTuple.cZZProng0(),
                                           trackTuple.cSnpYProng0(), trackTuple.cSnpZProng0(),
                                           trackTuple.cSnpSnpProng0(), trackTuple.cTglYProng0(), trackTuple.cTglZProng0(),
                                           trackTuple.cTglSnpProng0(), trackTuple.cTglTglProng0(),
                                           trackTuple.c1PtYProng0(), trackTuple.c1PtZProng0(), trackTuple.c1PtSnpProng0(),
                                           trackTuple.c1PtTglProng0(), trackTuple.c1Pt21Pt2Prong0()};

    return o2::track::TrackParametrizationWithError<TrackPrecision>(trackTuple.xProng0(), trackTuple.alphaProng0(), std::move(arrayparProng0), std::move(covparProng0));
}

/// Constructs TrackParametrizationWithError instance from trackTuple
/// \param trackuple is a row of joined tables HFTrack<0,1,2,3>
/// \return TrackParametrizationWithError object for prong1
template <typename TrackPrecision = float, typename T>
o2::track::TrackParametrizationWithError<TrackPrecision> getTrackParCov1(const T& trackTuple)
{
    std::array<TrackPrecision, 5> arrayparProng1 = {trackTuple.yProng1(), trackTuple.zProng1(), trackTuple.snpProng1(),
                                            trackTuple.tglProng1(), trackTuple.signed1PtProng1()};
    std::array<TrackPrecision, 15> covparProng1 = {trackTuple.cYYProng1(), trackTuple.cZYProng1(), trackTuple.cZZProng1(),
                                           trackTuple.cSnpYProng1(), trackTuple.cSnpZProng1(),
                                           trackTuple.cSnpSnpProng1(), trackTuple.cTglYProng1(), trackTuple.cTglZProng1(),
                                           trackTuple.cTglSnpProng1(), trackTuple.cTglTglProng1(),
                                           trackTuple.c1PtYProng1(), trackTuple.c1PtZProng1(), trackTuple.c1PtSnpProng1(),
                                           trackTuple.c1PtTglProng1(), trackTuple.c1Pt21Pt2Prong1()};

    return o2::track::TrackParametrizationWithError<TrackPrecision>(trackTuple.xProng1(), trackTuple.alphaProng1(), std::move(arrayparProng1), std::move(covparProng1));
}

/// Constructs TrackParametrizationWithError instance from trackTuple
/// \param trackuple is a row of joined tables HFTrack<0,1,2,3>
/// \return TrackParametrizationWithError object for prong2
template <typename TrackPrecision = float, typename T>
o2::track::TrackParametrizationWithError<TrackPrecision> getTrackParCov2(const T& trackTuple)
{
    std::array<TrackPrecision, 5> arrayparProng2 = {trackTuple.yProng2(), trackTuple.zProng2(), trackTuple.snpProng2(),
                                            trackTuple.tglProng2(), trackTuple.signed1PtProng2()};
    std::array<TrackPrecision, 15> covparProng2 = {trackTuple.cYYProng2(), trackTuple.cZYProng2(), trackTuple.cZZProng2(),
                                           trackTuple.cSnpYProng2(), trackTuple.cSnpZProng2(),
                                           trackTuple.cSnpSnpProng2(), trackTuple.cTglYProng2(), trackTuple.cTglZProng2(),
                                           trackTuple.cTglSnpProng2(), trackTuple.cTglTglProng2(),
                                           trackTuple.c1PtYProng2(), trackTuple.c1PtZProng2(), trackTuple.c1PtSnpProng2(),
                                           trackTuple.c1PtTglProng2(), trackTuple.c1Pt21Pt2Prong2()};

    return o2::track::TrackParametrizationWithError<TrackPrecision>(trackTuple.xProng2(), trackTuple.alphaProng2(), std::move(arrayparProng2), std::move(covparProng2));
}

/// Constructs TrackParametrizationWithError instance from trackTuple
/// \param trackuple is a row of joined tables HFTrack<0,1,2,3>
/// \return TrackParametrizationWithError object for prong3
template <typename TrackPrecision = float, typename T>
o2::track::TrackParametrizationWithError<TrackPrecision> getTrackParCov3(const T& trackTuple)
{
    std::array<TrackPrecision, 5> arrayparProng3 = {trackTuple.yProng3(), trackTuple.zProng3(), trackTuple.snpProng3(),
                                            trackTuple.tglProng3(), trackTuple.signed1PtProng3()};
    std::array<TrackPrecision, 15> covparProng3 = {trackTuple.cYYProng3(), trackTuple.cZYProng3(), trackTuple.cZZProng3(),
                                           trackTuple.cSnpYProng3(), trackTuple.cSnpZProng3(),
                                           trackTuple.cSnpSnpProng3(), trackTuple.cTglYProng3(), trackTuple.cTglZProng3(),
                                           trackTuple.cTglSnpProng3(), trackTuple.cTglTglProng3(),
                                           trackTuple.c1PtYProng3(), trackTuple.c1PtZProng3(), trackTuple.c1PtSnpProng3(),
                                           trackTuple.c1PtTglProng3(), trackTuple.c1Pt21Pt2Prong3()};

    return o2::track::TrackParametrizationWithError<TrackPrecision>(trackTuple.xProng3(), trackTuple.alphaProng3(), std::move(arrayparProng3), std::move(covparProng3));
}

/*
/// Constructs TrackParametrizationWithError instance from track attributes
/// \param trackParCovAttributes is a std::array<float, 22> containing track parameters and covariance matrix
/// \return TrackParametrizationWithError instance
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
*/

/*
// CAREFUL: the getters names shall be the same as the ones of the getTrackParCov method in Common/Core/trackUtilities.h
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
DECLARE_SOA_COLUMN(IsProng3, prong3, float); //!
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
*/
/*Prong0*/
DECLARE_SOA_COLUMN(XProng0, xProng0, float); //!
DECLARE_SOA_COLUMN(AlphaProng0, alphaProng0, float); //!
DECLARE_SOA_COLUMN(YProng0, yProng0, float); //!
DECLARE_SOA_COLUMN(ZProng0, zProng0, float); //!
DECLARE_SOA_COLUMN(SnpProng0, snpProng0, float); //!
DECLARE_SOA_COLUMN(TglProng0, tglProng0, float); //!
DECLARE_SOA_COLUMN(Q2PtProng0, signed1PtProng0, float); //!
DECLARE_SOA_COLUMN(SigY2Prong0, cYYProng0, float); //!
DECLARE_SOA_COLUMN(SigZYProng0, cZYProng0, float); //!
DECLARE_SOA_COLUMN(SigZ2Prong0, cZZProng0, float); //!
DECLARE_SOA_COLUMN(SigSnpYProng0, cSnpYProng0, float); //!
DECLARE_SOA_COLUMN(SigSnpZProng0, cSnpZProng0, float); //!
DECLARE_SOA_COLUMN(SigSnp2Prong0, cSnpSnpProng0, float); //!
DECLARE_SOA_COLUMN(SigTglYProng0, cTglYProng0, float); //!
DECLARE_SOA_COLUMN(SigTglZProng0, cTglZProng0, float); //!
DECLARE_SOA_COLUMN(SigTglSnpProng0, cTglSnpProng0, float); //!
DECLARE_SOA_COLUMN(SigTgl2Prong0, cTglTglProng0, float); //!
DECLARE_SOA_COLUMN(SigQ2PtYProng0, c1PtYProng0, float); //!
DECLARE_SOA_COLUMN(SigQ2PtZProng0, c1PtZProng0, float); //!
DECLARE_SOA_COLUMN(SigQ2PtSnpProng0, c1PtSnpProng0, float); //!
DECLARE_SOA_COLUMN(SigQ2PtTglProng0, c1PtTglProng0, float); //!
DECLARE_SOA_COLUMN(SigQ2Pt2Prong0, c1Pt21Pt2Prong0, float); //!
DECLARE_SOA_COLUMN(PxProng0, pxProng0, float); //!
DECLARE_SOA_COLUMN(PyProng0, pyProng0, float); //!
DECLARE_SOA_COLUMN(PzProng0, pzProng0, float); //!

/*Prong1*/
DECLARE_SOA_COLUMN(XProng1, xProng1, float); //!
DECLARE_SOA_COLUMN(AlphaProng1, alphaProng1, float); //!
DECLARE_SOA_COLUMN(YProng1, yProng1, float); //!
DECLARE_SOA_COLUMN(ZProng1, zProng1, float); //!
DECLARE_SOA_COLUMN(SnpProng1, snpProng1, float); //!
DECLARE_SOA_COLUMN(TglProng1, tglProng1, float); //!
DECLARE_SOA_COLUMN(Q2PtProng1, signed1PtProng1, float); //!
DECLARE_SOA_COLUMN(SigY2Prong1, cYYProng1, float); //!
DECLARE_SOA_COLUMN(SigZYProng1, cZYProng1, float); //!
DECLARE_SOA_COLUMN(SigZ2Prong1, cZZProng1, float); //!
DECLARE_SOA_COLUMN(SigSnpYProng1, cSnpYProng1, float); //!
DECLARE_SOA_COLUMN(SigSnpZProng1, cSnpZProng1, float); //!
DECLARE_SOA_COLUMN(SigSnp2Prong1, cSnpSnpProng1, float); //!
DECLARE_SOA_COLUMN(SigTglYProng1, cTglYProng1, float); //!
DECLARE_SOA_COLUMN(SigTglZProng1, cTglZProng1, float); //!
DECLARE_SOA_COLUMN(SigTglSnpProng1, cTglSnpProng1, float); //!
DECLARE_SOA_COLUMN(SigTgl2Prong1, cTglTglProng1, float); //!
DECLARE_SOA_COLUMN(SigQ2PtYProng1, c1PtYProng1, float); //!
DECLARE_SOA_COLUMN(SigQ2PtZProng1, c1PtZProng1, float); //!
DECLARE_SOA_COLUMN(SigQ2PtSnpProng1, c1PtSnpProng1, float); //!
DECLARE_SOA_COLUMN(SigQ2PtTglProng1, c1PtTglProng1, float); //!
DECLARE_SOA_COLUMN(SigQ2Pt2Prong1, c1Pt21Pt2Prong1, float); //!
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float); //!
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float); //!
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float); //!

/*Prong2*/
DECLARE_SOA_COLUMN(XProng2, xProng2, float); //!
DECLARE_SOA_COLUMN(AlphaProng2, alphaProng2, float); //!
DECLARE_SOA_COLUMN(YProng2, yProng2, float); //!
DECLARE_SOA_COLUMN(ZProng2, zProng2, float); //!
DECLARE_SOA_COLUMN(SnpProng2, snpProng2, float); //!
DECLARE_SOA_COLUMN(TglProng2, tglProng2, float); //!
DECLARE_SOA_COLUMN(Q2PtProng2, signed1PtProng2, float); //!
DECLARE_SOA_COLUMN(SigY2Prong2, cYYProng2, float); //!
DECLARE_SOA_COLUMN(SigZYProng2, cZYProng2, float); //!
DECLARE_SOA_COLUMN(SigZ2Prong2, cZZProng2, float); //!
DECLARE_SOA_COLUMN(SigSnpYProng2, cSnpYProng2, float); //!
DECLARE_SOA_COLUMN(SigSnpZProng2, cSnpZProng2, float); //!
DECLARE_SOA_COLUMN(SigSnp2Prong2, cSnpSnpProng2, float); //!
DECLARE_SOA_COLUMN(SigTglYProng2, cTglYProng2, float); //!
DECLARE_SOA_COLUMN(SigTglZProng2, cTglZProng2, float); //!
DECLARE_SOA_COLUMN(SigTglSnpProng2, cTglSnpProng2, float); //!
DECLARE_SOA_COLUMN(SigTgl2Prong2, cTglTglProng2, float); //!
DECLARE_SOA_COLUMN(SigQ2PtYProng2, c1PtYProng2, float); //!
DECLARE_SOA_COLUMN(SigQ2PtZProng2, c1PtZProng2, float); //!
DECLARE_SOA_COLUMN(SigQ2PtSnpProng2, c1PtSnpProng2, float); //!
DECLARE_SOA_COLUMN(SigQ2PtTglProng2, c1PtTglProng2, float); //!
DECLARE_SOA_COLUMN(SigQ2Pt2Prong2, c1Pt21Pt2Prong2, float); //!
DECLARE_SOA_COLUMN(PxProng2, pxProng2, float); //!
DECLARE_SOA_COLUMN(PyProng2, pyProng2, float); //!
DECLARE_SOA_COLUMN(PzProng2, pzProng2, float); //!

/*Prong3*/
DECLARE_SOA_COLUMN(XProng3, xProng3, float); //!
DECLARE_SOA_COLUMN(AlphaProng3, alphaProng3, float); //!
DECLARE_SOA_COLUMN(YProng3, yProng3, float); //!
DECLARE_SOA_COLUMN(ZProng3, zProng3, float); //!
DECLARE_SOA_COLUMN(SnpProng3, snpProng3, float); //!
DECLARE_SOA_COLUMN(TglProng3, tglProng3, float); //!
DECLARE_SOA_COLUMN(Q2PtProng3, signed1PtProng3, float); //!
DECLARE_SOA_COLUMN(SigY2Prong3, cYYProng3, float); //!
DECLARE_SOA_COLUMN(SigZYProng3, cZYProng3, float); //!
DECLARE_SOA_COLUMN(SigZ2Prong3, cZZProng3, float); //!
DECLARE_SOA_COLUMN(SigSnpYProng3, cSnpYProng3, float); //!
DECLARE_SOA_COLUMN(SigSnpZProng3, cSnpZProng3, float); //!
DECLARE_SOA_COLUMN(SigSnp2Prong3, cSnpSnpProng3, float); //!
DECLARE_SOA_COLUMN(SigTglYProng3, cTglYProng3, float); //!
DECLARE_SOA_COLUMN(SigTglZProng3, cTglZProng3, float); //!
DECLARE_SOA_COLUMN(SigTglSnpProng3, cTglSnpProng3, float); //!
DECLARE_SOA_COLUMN(SigTgl2Prong3, cTglTglProng3, float); //!
DECLARE_SOA_COLUMN(SigQ2PtYProng3, c1PtYProng3, float); //!
DECLARE_SOA_COLUMN(SigQ2PtZProng3, c1PtZProng3, float); //!
DECLARE_SOA_COLUMN(SigQ2PtSnpProng3, c1PtSnpProng3, float); //!
DECLARE_SOA_COLUMN(SigQ2PtTglProng3, c1PtTglProng3, float); //!
DECLARE_SOA_COLUMN(SigQ2Pt2Prong3, c1Pt21Pt2Prong3, float); //!
DECLARE_SOA_COLUMN(PxProng3, pxProng3, float); //!
DECLARE_SOA_COLUMN(PyProng3, pyProng3, float); //!
DECLARE_SOA_COLUMN(PzProng3, pzProng3, float); //!

/*?*/
DECLARE_SOA_COLUMN(IsProng0, prong0, float); //!
DECLARE_SOA_COLUMN(IsProng1, prong1, float); //!
DECLARE_SOA_COLUMN(IsProng2, prong2, float); //!
DECLARE_SOA_COLUMN(IsProng3, prong3, float); //!
} // namespace hf_track_par_cov

DECLARE_SOA_TABLE(HfTrack0, "AOD", "HFTRACK0", //!
                hf_track_par_cov::XProng0,
                hf_track_par_cov::AlphaProng0,
                hf_track_par_cov::YProng0,
                hf_track_par_cov::ZProng0,
                hf_track_par_cov::SnpProng0,
                hf_track_par_cov::TglProng0,
                hf_track_par_cov::Q2PtProng0,
                hf_track_par_cov::SigY2Prong0,
                hf_track_par_cov::SigZYProng0,
                hf_track_par_cov::SigZ2Prong0,
                hf_track_par_cov::SigSnpYProng0,
                hf_track_par_cov::SigSnpZProng0,
                hf_track_par_cov::SigSnp2Prong0,
                hf_track_par_cov::SigTglYProng0,
                hf_track_par_cov::SigTglZProng0,
                hf_track_par_cov::SigTglSnpProng0,
                hf_track_par_cov::SigTgl2Prong0,
                hf_track_par_cov::SigQ2PtYProng0,
                hf_track_par_cov::SigQ2PtZProng0,
                hf_track_par_cov::SigQ2PtSnpProng0,
                hf_track_par_cov::SigQ2PtTglProng0,
                hf_track_par_cov::SigQ2Pt2Prong0,
                hf_track_par_cov::PxProng0,
                hf_track_par_cov::PyProng0,
                hf_track_par_cov::PzProng0);

DECLARE_SOA_TABLE(HfTrack1, "AOD", "HFTRACK1", //!
                hf_track_par_cov::XProng1,
                hf_track_par_cov::AlphaProng1,
                hf_track_par_cov::YProng1,
                hf_track_par_cov::ZProng1,
                hf_track_par_cov::SnpProng1,
                hf_track_par_cov::TglProng1,
                hf_track_par_cov::Q2PtProng1,
                hf_track_par_cov::SigY2Prong1,
                hf_track_par_cov::SigZYProng1,
                hf_track_par_cov::SigZ2Prong1,
                hf_track_par_cov::SigSnpYProng1,
                hf_track_par_cov::SigSnpZProng1,
                hf_track_par_cov::SigSnp2Prong1,
                hf_track_par_cov::SigTglYProng1,
                hf_track_par_cov::SigTglZProng1,
                hf_track_par_cov::SigTglSnpProng1,
                hf_track_par_cov::SigTgl2Prong1,
                hf_track_par_cov::SigQ2PtYProng1,
                hf_track_par_cov::SigQ2PtZProng1,
                hf_track_par_cov::SigQ2PtSnpProng1,
                hf_track_par_cov::SigQ2PtTglProng1,
                hf_track_par_cov::SigQ2Pt2Prong1,
                hf_track_par_cov::PxProng1,
                hf_track_par_cov::PyProng1,
                hf_track_par_cov::PzProng1);

DECLARE_SOA_TABLE(HfTrack2, "AOD", "HFTRACK2", //!
                hf_track_par_cov::XProng2,
                hf_track_par_cov::AlphaProng2,
                hf_track_par_cov::YProng2,
                hf_track_par_cov::ZProng2,
                hf_track_par_cov::SnpProng2,
                hf_track_par_cov::TglProng2,
                hf_track_par_cov::Q2PtProng2,
                hf_track_par_cov::SigY2Prong2,
                hf_track_par_cov::SigZYProng2,
                hf_track_par_cov::SigZ2Prong2,
                hf_track_par_cov::SigSnpYProng2,
                hf_track_par_cov::SigSnpZProng2,
                hf_track_par_cov::SigSnp2Prong2,
                hf_track_par_cov::SigTglYProng2,
                hf_track_par_cov::SigTglZProng2,
                hf_track_par_cov::SigTglSnpProng2,
                hf_track_par_cov::SigTgl2Prong2,
                hf_track_par_cov::SigQ2PtYProng2,
                hf_track_par_cov::SigQ2PtZProng2,
                hf_track_par_cov::SigQ2PtSnpProng2,
                hf_track_par_cov::SigQ2PtTglProng2,
                hf_track_par_cov::SigQ2Pt2Prong2,
                hf_track_par_cov::PxProng2,
                hf_track_par_cov::PyProng2,
                hf_track_par_cov::PzProng2);

DECLARE_SOA_TABLE(HfTrack3, "AOD", "HFTRACK3", //!
                hf_track_par_cov::XProng3,
                hf_track_par_cov::AlphaProng3,
                hf_track_par_cov::YProng3,
                hf_track_par_cov::ZProng3,
                hf_track_par_cov::SnpProng3,
                hf_track_par_cov::TglProng3,
                hf_track_par_cov::Q2PtProng3,
                hf_track_par_cov::SigY2Prong3,
                hf_track_par_cov::SigZYProng3,
                hf_track_par_cov::SigZ2Prong3,
                hf_track_par_cov::SigSnpYProng3,
                hf_track_par_cov::SigSnpZProng3,
                hf_track_par_cov::SigSnp2Prong3,
                hf_track_par_cov::SigTglYProng3,
                hf_track_par_cov::SigTglZProng3,
                hf_track_par_cov::SigTglSnpProng3,
                hf_track_par_cov::SigTgl2Prong3,
                hf_track_par_cov::SigQ2PtYProng3,
                hf_track_par_cov::SigQ2PtZProng3,
                hf_track_par_cov::SigQ2PtSnpProng3,
                hf_track_par_cov::SigQ2PtTglProng3,
                hf_track_par_cov::SigQ2Pt2Prong3,
                hf_track_par_cov::PxProng3,
                hf_track_par_cov::PyProng3,
                hf_track_par_cov::PzProng3);


} // namespace o2::aod

#endif // PWGHF_DATAMODEL_DERIVED_TABLES_FOR_SKIMMING_H_