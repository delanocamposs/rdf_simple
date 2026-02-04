#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
using namespace ROOT;
using namespace ROOT::VecOps;

RVecF getEnergy(RVecF pt, RVecF eta, RVecF phi, RVecF mass){
  RVecF out;
  out.reserve(pt.size());
  for (size_t i = 0; i < pt.size(); i++){
    ROOT::Math::PtEtaPhiMVector p(pt[i], eta[i], phi[i], mass[i]);
    out.emplace_back(p.energy());
  }
  return out;
}

RVecI decayToSignal(RVecI gen_pdgId, RVecI signal_motherIdx){
  RVecI out(gen_pdgId.size(), 0);
  for (size_t i = 0; i < signal_motherIdx.size(); i++){
    if (signal_motherIdx[i] != -1)
      out[signal_motherIdx[i]]++;
  }
  return out;
}

RVecI matchSignal(RVecF sig_energy, RVecF sig_eta, RVecF sig_phi, RVecF g_energy, RVecF g_eta, RVecF g_phi, const float dr_max){

  RVec<int> out(sig_energy.size(), -1);
  RVec<std::pair<size_t, float>> sig_info;  // Create structure to pair indices to energies (for sorting)
  for (size_t i = 0; i < sig_energy.size(); i++){
    std::pair<size_t, float> info;
    info.first = i;
    info.second = sig_energy[i];
    sig_info.emplace_back(info);
  }
  auto sig_info_sorted = Sort(sig_info, [](std::pair<size_t, float> a, std::pair<size_t, float> b) {return a.second > b.second;});

  std::vector<size_t> idx_matchedReco; // vector of reco photon indices that are matched

  // Start by matching the highest energy
  for (size_t i = 0; i < sig_info_sorted.size(); i++){
    size_t idx = sig_info_sorted[i].first;
    RVec<std::pair<size_t, float>> matches;

    // Get all matching reco photons
    for (size_t j = 0; j < g_energy.size(); j++){
      float dr = DeltaR(sig_eta[idx], g_eta[j], sig_phi[idx], g_phi[j]);
      if ((dr < dr_max) && (sig_energy[idx] > 0.5*g_energy[j]) && (sig_energy[idx] < 2.0*g_energy[j])){
	std::pair<size_t, float> info;
	info.first = j;
	info.second = abs(sig_energy[idx] - g_energy[j]);
	matches.emplace_back(info);
      }
    }

    // Sort matches by closest energy
    auto matches_sorted = Sort(matches, [](std::pair<size_t, float> a, std::pair<size_t, float> b) {return a.second < b.second;});


    // Check if best reco match has already been used
    for (size_t j = 0; j < matches_sorted.size(); j++){
      // matched to unused reco photon
      if (std::find(idx_matchedReco.begin(), idx_matchedReco.end(), matches_sorted[j].first) == idx_matchedReco.end()){
	out[sig_info_sorted[i].first] = matches_sorted[j].first;
	idx_matchedReco.push_back(matches_sorted[j].first);
	break;
      }
      // matched to used reco photon, try next best match if possible
      else{
	continue;
      }
    }
  }
  return out;
}

RVecI pairMatchedToReco(RVecI gen_idx, RVecI sig_recoIdx, RVecI sig_motherIdx){
  RVecI out(gen_idx.size(), 0);
  for (size_t i = 0; i < sig_motherIdx.size(); i++){
    if (sig_recoIdx[i] != -1 && sig_motherIdx[i] != -1)
      out[sig_motherIdx[i]]++;
  }
  return out;
}

RVecF getDecayDxy(RVecF vx, RVecF vy, RVecI motherIdx){
  RVecF out(vx.size(), 0);
  for (size_t i = 0; i < vx.size(); i++){
    if (motherIdx[i] == -1)
      continue;
    float dx = vx[motherIdx[i]] - vx[i];
    float dy = vy[motherIdx[i]] - vy[i];
    float dxy = std::sqrt(vx[motherIdx[i]]*vx[motherIdx[i]] + vy[motherIdx[i]]*vy[motherIdx[i]]);
    out[motherIdx[i]] = dxy;
  }
  return out;
}
//TODO Duplucate this function 
// anti-ID: fail by 
RVecI passIDNoIso(RVecI bitmaps){
  RVecI out;
  out.reserve(bitmaps.size());
  for (size_t i = 0; i < bitmaps.size(); i++){
    int bitmap = bitmaps[i];
    bool passSIEIE = (bitmap>>6&3) >= 1;
    bool passHOE = (bitmap>>4&3) >= 1;
    bool passChIso = (bitmap>>8&3) >= 1;
    bool passNeuIso = (bitmap>>10&3) >= 1;
    out.emplace_back(passSIEIE && passHOE && passChIso && passNeuIso);
  }
  return out;
}

RVecI failIDby_HoE_SIEIE(RVecI bitmaps){
  RVecI out;
  out.reserve(bitmaps.size());
  for (size_t i = 0; i < bitmaps.size(); i++){
    int bitmap = bitmaps[i];
    bool failSIEIE = (bitmap>>6&3) == 0;
    bool failHOE = (bitmap>>4&3) == 0;
    bool passChIso = (bitmap>>8&3) >= 1;
    bool passNeuIso = (bitmap>>10&3) >= 1;
    bool passPhIso = (bitmap>>12&3) >= 1;
    out.emplace_back(failSIEIE && failHOE && passChIso && passNeuIso && passPhIso);
  }
  return out;
}

RVecI passHOE(RVecI bitmaps){
  RVecI out;
  out.reserve(bitmaps.size());
  for (size_t i = 0; i < bitmaps.size(); i++){
    int bitmap = bitmaps[i];
    bool pass = (bitmap>>4&3) >= 1;
    out.emplace_back(pass);
  }
  return out;
}

RVecI passSieie(RVecI bitmaps){
  RVecI out;
  out.reserve(bitmaps.size());
  for (size_t i = 0; i < bitmaps.size(); i++){
    int bitmap = bitmaps[i];
    bool pass = (bitmap>>6&3) >= 1;
    out.emplace_back(pass);
  }
  return out;
}

RVecI passChIso(RVecI bitmaps){
  RVecI out;
  out.reserve(bitmaps.size());
  for (size_t i = 0; i < bitmaps.size(); i++){
    int bitmap = bitmaps[i];
    bool pass = (bitmap>>8&3) >= 1;
    out.emplace_back(pass);
  }
  return out;
}

RVecI passNeuIso(RVecI bitmaps){
  RVecI out;
  out.reserve(bitmaps.size());
  for (size_t i = 0; i < bitmaps.size(); i++){
    int bitmap = bitmaps[i];
    bool pass = (bitmap>>10&3) >= 1;
    out.emplace_back(pass);
  }
  return out;
}

RVecI passPhIso(RVecI bitmaps){
  RVecI out;
  out.reserve(bitmaps.size());
  for (size_t i = 0; i < bitmaps.size(); i++){
    int bitmap = bitmaps[i];
    bool pass = (bitmap>>12&3) >= 1;
    out.emplace_back(pass);
  }
  return out;
}

RVecI passPhIso_Run3(RVecI bitmaps){
  RVecI out;
  out.reserve(bitmaps.size());
  for (size_t i = 0; i < bitmaps.size(); i++){
    int bitmap = bitmaps[i];
    bool pass = (bitmap>>12&3) >= 1;
    out.emplace_back(pass);
  }
  return out;
}
