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

/// \file HfMlResponseDplusToPiKPi.h
/// \brief Class to compute the ML response for D± → π± K∓ π± analysis selections
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg

#ifndef PWGHF_CORE_HFMLRESPONSEDPLUSTOPIKPI_H_
#define PWGHF_CORE_HFMLRESPONSEDPLUSTOPIKPI_H_

#include <string>
#include <vector>

#include "PWGHF/Core/HfMlResponse.h"

namespace o2::analysis
{

enum class InputFeaturesDplusToPiKPi : int8_t {
  ptProng0 = 0,
  ptProng1,
  ptProng2,
  impactParameter0,
  impactParameter1,
  impactParameter2,
  decayLength,
  decayLengthXYNormalised,
  cpa,
  cpaXY,
  maxNormalisedDeltaIP,
  tpcTofNSigmaPi0,
  tpcTofNSigmaPi1,
  tpcTofNSigmaPi2,
  tpcTofNSigmaKa0,
  tpcTofNSigmaKa1,
  tpcTofNSigmaKa2,
  nFeatures
};

template <typename T = float>
class HfMlResponseDplusToPiKPi : public o2::analysis::HfMlResponse<T>
{
 public:
  /// Default constructor
  HfMlResponseDplusToPiKPi() = default;
  /// Default destructor
  virtual ~HfMlResponseDplusToPiKPi() = default;  

  /// Method to translate input features strings into int
  /// \param cfgInputFeatures array of input features names
  void cacheInputFeaturesIndices(std::vector<std::string> const& cfgInputFeatures)
  {
    uint8_t counter{0};
    auto begin = availableInputFeatures.begin();
    auto end = availableInputFeatures.end();
    for (const auto& inputFeature : cfgInputFeatures) {
      auto itr = std::find(begin, end, inputFeature);
      if (itr != end) {
        // link feature's index in InputFeaturesDplusToPiKPi and cfgInputFeatures
        mCachedIndices.emplace_back((int8_t)std::distance(begin, itr));
      } else {
        LOG(fatal) << "Input feature not available. Please check your configurables.";
      }
      ++counter;
    }
  }

#define FILL_VEC(GETTER) \
  case InputFeaturesDplusToPiKPi::GETTER: \
    inputFeatures.emplace_back(candidate.GETTER()); \
    break

#define FILL_VEC_FULL(OBJECT,FEATURE,GETTER) \
  case InputFeaturesDplusToPiKPi::FEATURE: \
    inputFeatures.emplace_back(OBJECT.GETTER()); \
    break

  /// Method to get the input features vector needed for ML inference
  /// \param candidate is the Dplus candidate
  /// \param prong0 is the candidate's prong0
  /// \param prong1 is the candidate's prong1
  /// \param prong2 is the candidate's prong2
  /// \return inputFeatures vector
  template <typename T1, typename T2>
  std::vector<T> getInputFeatures(T1 const& candidate, T2 const& prong0, T2 const& prong1, T2 const& prong2)
  {
    std::vector<T> inputFeatures;

    for (const auto& idx : mCachedIndices) {
      switch (static_cast<InputFeaturesDplusToPiKPi>(idx)) {
        FILL_VEC(ptProng0);
        FILL_VEC(ptProng1);
        FILL_VEC(ptProng2);
        FILL_VEC(impactParameter0);
        FILL_VEC(impactParameter1);
        FILL_VEC(impactParameter2);
        FILL_VEC(decayLength);
        FILL_VEC(decayLengthXYNormalised);
        FILL_VEC(cpa);
        FILL_VEC(cpaXY);
        FILL_VEC(maxNormalisedDeltaIP);
        FILL_VEC_FULL(prong0, tpcTofNSigmaPi0, tpcTofNSigmaPi);
        FILL_VEC_FULL(prong1, tpcTofNSigmaPi1, tpcTofNSigmaPi);
        FILL_VEC_FULL(prong2, tpcTofNSigmaPi2, tpcTofNSigmaPi);
        FILL_VEC_FULL(prong0, tpcTofNSigmaKa0, tpcTofNSigmaKa);
        FILL_VEC_FULL(prong1, tpcTofNSigmaKa1, tpcTofNSigmaKa);
        FILL_VEC_FULL(prong2, tpcTofNSigmaKa2, tpcTofNSigmaKa);
      }
    }

    return inputFeatures;
  }

 private:
  const std::array<std::string, static_cast<int>(InputFeaturesDplusToPiKPi::nFeatures)> availableInputFeatures = {
    "ptProng0",
    "ptProng1",
    "ptProng2",
    "impactParameter0",
    "impactParameter1",
    "impactParameter2",
    "decayLength",
    "decayLengthXYNormalised",
    "cpa",
    "cpaXY",
    "maxNormalisedDeltaIP",
    "tpcTofNSigmaPi0",
    "tpcTofNSigmaPi1",
    "tpcTofNSigmaPi2",
    "tpcTofNSigmaKa0",
    "tpcTofNSigmaKa1",
    "tpcTofNSigmaKa2"
  };
  std::vector<int8_t> mCachedIndices;
};

} // namespace o2::analysis

#endif // PWGHF_CORE_HFMLRESPONSEDPLUSTOPIKPI_H_
