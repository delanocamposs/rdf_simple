#ifndef CHELPERS
#define CHELPERS

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
using namespace ROOT;
using namespace ROOT::VecOps;
float Z_mass=91;


using RVecB = ROOT::VecOps::RVec<bool>;
using RVecC = ROOT::VecOps::RVec<char>;
using RVecD = ROOT::VecOps::RVec<double>;
using RVecF = ROOT::VecOps::RVec<float>;
using RVecI = ROOT::VecOps::RVec<int>;
using RVecL = ROOT::VecOps::RVec<long int>;
using RVecLL = ROOT::VecOps::RVec<long long int>;
using RVecU = ROOT::VecOps::RVec<unsigned int>;
using RVecUL = ROOT::VecOps::RVec<unsigned long int>;
using RVecULL = ROOT::VecOps::RVec<unsigned long long int>;

// Reconstruct the best Z-> ll candidate
RVec<size_t> best_z(RVecF pt, RVecF eta, RVecF phi, RVecF mass, RVecI charge, RVecB isTight)
{

  RVec<size_t> result;
  result.reserve(2);
  auto idx_cmb = ROOT::VecOps::Combinations(pt, 2);
  auto best_mass = -1;
  size_t best_i1=0;
  size_t best_i2=0;
  for (size_t i = 0; i < idx_cmb[0].size(); i++) {
    const auto i1 = idx_cmb[0][i];
    const auto i2 = idx_cmb[1][i];
    if (!(isTight[i1]&&isTight[i2]))
      continue;
    if (charge[i1] != charge[i2]) {
      ROOT::Math::PtEtaPhiMVector p1(pt[i1], eta[i1], phi[i1], mass[i1]);
      ROOT::Math::PtEtaPhiMVector p2(pt[i2], eta[i2], phi[i2], mass[i2]);
      const auto this_mass = (p1 + p2).M();
      if (std::abs(this_mass-Z_mass)<std::abs(best_mass-Z_mass)) {
	best_mass=Z_mass;
	best_i1=i1;
	best_i2=i2;
      }
    }
  }
  result.emplace_back(best_i1);
  result.emplace_back(best_i2);
  return result;
}

// Picking best leptons for FSR Z events
RVec<size_t> best_zg(RVecF pt, RVecF eta, RVecF phi, RVecF mass, RVecI charge, RVecB isTight, RVecF g_pt, RVecF g_eta, RVecF g_phi)
{
  RVec<size_t> result;
  result.reserve(3);
  auto idx_cmb = ROOT::VecOps::Combinations(pt, 2);
  auto best_mass = -999;
  size_t best_i1 = 0;
  size_t best_i2 = 0;
  size_t best_g = 0;

  for (size_t i = 0; i < idx_cmb[0].size(); i++) { // Loop through pairs of leptons
    const auto i1 = idx_cmb[0][i];
    const auto i2 = idx_cmb[1][i];
    if (!(isTight[i1] && isTight[i2])) // Require both passing ID
      continue;
    if (charge[i1] != charge[i2]) { // Require opposite charge
      ROOT::Math::PtEtaPhiMVector p1(pt[i1], eta[i1], phi[i1], mass[i1]);
      ROOT::Math::PtEtaPhiMVector p2(pt[i2], eta[i2], phi[i2], mass[i2]);
      for (size_t j = 0; j < g_pt.size(); j++) { // Loop through photons to find closest inv mass to 90 GeV
	ROOT::Math::PtEtaPhiMVector g(g_pt[j], g_eta[j], g_phi[j], 0);
	const auto this_mass = (p1 + p2 + g).M();
	if (std::abs(this_mass - Z_mass) < std::abs(best_mass - Z_mass)) {
	  best_mass = Z_mass;
	  best_i1 = i1;
	  best_i2 = i2;
	  best_g = j;
	}
      }
    }
  }
  result.emplace_back(best_i1);
  result.emplace_back(best_i2);
  result.emplace_back(best_g);
  return result;
}

RVecF Zgg_fsr(RVecF g_pt, RVecF g_eta, RVecF g_phi, const int idx1, const int idx2, RVecF l_pt, RVecF l_eta, RVecF l_phi, RVecF l_mass, RVec<size_t> z_idx){
  RVecF out;
  out.reserve(3);
  ROOT::Math::PtEtaPhiMVector g1(g_pt[idx1], g_eta[idx1], g_phi[idx1], 0);
  ROOT::Math::PtEtaPhiMVector g2(g_pt[idx2], g_eta[idx2], g_phi[idx2], 0);
  ROOT::Math::PtEtaPhiMVector l1(l_pt[z_idx[0]], l_eta[z_idx[0]], l_phi[z_idx[0]], l_mass[z_idx[0]]);
  ROOT::Math::PtEtaPhiMVector l2(l_pt[z_idx[1]], l_eta[z_idx[1]], l_phi[z_idx[1]], l_mass[z_idx[1]]);
  out.emplace_back((l1+l2+g1).mass());
  out.emplace_back((l1+l2+g2).mass());
  out.emplace_back((l1+l2+g1+g2).mass());
  return out;
}

RVec<bool> overlapClean(RVecF gphi, RVecF geta, RVecF lphi, RVecF leta){
  RVec<bool> out;
  out.reserve(gphi.size());
  for (size_t i = 0; i < gphi.size(); i++){
    bool overlap = false;
    for (size_t j = 0; j < lphi.size(); j++){
      if (pow(DeltaPhi(gphi[i], lphi[j]),2)/pow(0.5,2) + pow(abs(geta[i] - leta[j]), 2)/pow(0.4,2) < 1){
	overlap = true;
	break;
      }
    }
    out.emplace_back(overlap);
  }
  return out;
}

RVecF best_Z_info(RVecF pt, RVecF eta, RVecF phi, RVecF mass, const RVec<size_t>& idx){
 
  RVecF out;
  out.reserve(6);
  size_t id1 = idx[0];
  size_t id2 = idx[1];
  ROOT::Math::PtEtaPhiMVector p0(pt[id1], eta[id1], phi[id1], mass[id1]);
  ROOT::Math::PtEtaPhiMVector p1(pt[id2], eta[id2], phi[id2], mass[id2]);
  out.emplace_back((p0+p1).pt());
  out.emplace_back((p0+p1).eta());
  out.emplace_back((p0+p1).phi());
  out.emplace_back((p0+p1).mass());
  out.emplace_back(DeltaR(eta[id1], eta[id2], phi[id1], phi[id2]));
  out.emplace_back(DeltaPhi(phi[id1], phi[id2]));
  return out;
}

RVecF best_Zg_info(RVecF pt, RVecF eta, RVecF phi, RVecF mass, RVecF g_pt, RVecF g_eta, RVecF g_phi, const RVec<size_t>& idx){
  RVecF out;
  out.reserve(9);
  size_t idx_l1 = idx[0];
  size_t idx_l2 = idx[1];
  size_t idx_g = idx[2];
  ROOT::Math::PtEtaPhiMVector l1(pt[idx_l1], eta[idx_l1], phi[idx_l1], mass[idx_l1]);
  ROOT::Math::PtEtaPhiMVector l2(pt[idx_l2], eta[idx_l2], phi[idx_l2], mass[idx_l2]);
  ROOT::Math::PtEtaPhiMVector g(g_pt[idx_g], g_eta[idx_g], g_phi[idx_g], 0);
  out.emplace_back((l1+l2+g).pt());
  out.emplace_back((l1+l2+g).eta());
  out.emplace_back((l1+l2+g).phi());
  out.emplace_back((l1+l2+g).mass());
  out.emplace_back(DeltaR(eta[idx_l1], g_eta[idx_g], phi[idx_l1], g_phi[idx_g]));
  out.emplace_back(DeltaR(eta[idx_l2], g_eta[idx_g], phi[idx_l2], g_phi[idx_g]));
  out.emplace_back(DeltaPhi(phi[idx_l1], g_phi[idx_g]));
  out.emplace_back(DeltaPhi(phi[idx_l2], g_phi[idx_g]));
  out.emplace_back((l1+l2).mass());
  return out;
}

RVecF best_W_info(RVecF pt, RVecF eta, RVecF phi, RVecF mass, RVecB isTight, const float met_pt, const float met_phi){
  RVecF out;
  out.reserve(8);
  // Find highest pT tight lepton
  auto idx = Argsort(pt, [](double x, double y) {return x > y;});
  for (size_t i = 0; i < pt.size(); i++){
    if (isTight[idx[i]]){
      ROOT::Math::PtEtaPhiMVector l(pt[idx[i]], eta[idx[i]], phi[idx[i]], mass[idx[i]]);
      ROOT::Math::PtEtaPhiMVector MET(met_pt, 0.0, met_phi, 0.0);
      out.emplace_back((l+MET).pt());
      out.emplace_back((l+MET).eta());
      out.emplace_back((l+MET).phi());
      out.emplace_back((l+MET).mass());
      float met_eta = 0.0;
      out.emplace_back(DeltaR(eta[idx[i]], met_eta, phi[idx[i]], met_phi));
      out.emplace_back(DeltaPhi(phi[idx[i]], met_phi));
      float mt2 = l.mass()*l.mass()+2*(l.Et()*MET.Et() - l.px()*MET.px() - l.py()*MET.py());
      out.emplace_back((mt2 > 0) ? std::sqrt(mt2) : 0.0);
      out.emplace_back(idx[i]);
      return out;
    }
  }
  for (int i = 0; i < 8; i++)
    out.emplace_back(0);
  return out;
}

RVecF correct_gammaIso_for_muons(const RVec<size_t>& idx,RVecF mpt, RVecF meta, RVecF mphi,RVecF gpt, RVecF geta, RVecF gphi,RVecF giso,RVec<bool> fsr) {
  RVecF result;  
  result.reserve(gpt.size());
  //add the FSR
  
  for(size_t i=0;i<gpt.size();++i) {
    float iso=giso[i]*gpt[i];
    if (fsr[i]){
      if (DeltaR(geta[i],meta[idx[0]],gphi[i],mphi[idx[0]])<0.3)
	iso=iso-mpt[idx[0]];
      if (DeltaR(geta[i],meta[idx[1]],gphi[i],mphi[idx[1]])<0.3)
	iso=iso-mpt[idx[1]];
    }
   
    if(iso<0)
      iso=0;
    result.emplace_back(iso/gpt[i]);
  }
  return result;
}

RVecF correct_gammaIso_for_photons(const int idx1, const int idx2, RVecF gpt, RVecF geta, RVecF gphi, RVecF giso) {
  RVecF result = giso;  
  if (DeltaR(geta[idx1], geta[idx2], gphi[idx1], gphi[idx2]) < 0.3){
    float iso1 = giso[idx1]*gpt[idx1];
    float iso2 = giso[idx2]*gpt[idx2];
    result[idx1] = (iso1 - gpt[idx2])/gpt[idx1] > 0.0 ? (iso1 - gpt[idx2])/gpt[idx1] : 0.0;
    result[idx2] = (iso2 - gpt[idx1])/gpt[idx2] > 0.0 ? (iso2 - gpt[idx1])/gpt[idx2] : 0.0;
  }
  return result;
}


RVecF correct_gammaIso(const RVecF& gpt, const RVecF& geta, const RVecF& gphi, const RVecF& giso, const RVec<bool>& gID) {
  RVecF result;
  result.reserve(gpt.size());

  for (size_t i = 0; i < gpt.size(); ++i) {
    float iso = giso[i] * gpt[i];

    // Correction for other photons
    for (size_t j = 0; j < gpt.size(); ++j) {
      if ((i == j) || gID[j] == 0)
        continue;

      if (DeltaR(geta[i], geta[j], gphi[i], gphi[j]) < 0.3)
        iso = iso - gpt[j];
    }

    // Ensure the corrected isolation is non-negative
    if (iso < 0)
      iso = 0;

    result.emplace_back(iso / gpt[i]);
  }

  return result;
}

RVecF correct_muoniso_for_photons(RVecF mpt, RVecF meta, RVecF mphi,RVecF miso,RVecF gpt, RVecF geta, RVecF gphi,RVec<bool> gFSR) {
  RVecF result;  
  result.reserve(mpt.size());
  //add the FSR
  
  for(size_t i=0;i<mpt.size();++i) {
    float iso=miso[i]*mpt[i];
    for(size_t j=0;j<gpt.size();++j) {
      if (gFSR[j]==0)
	continue;
      
      if (DeltaR(meta[i],geta[j],mphi[i],gphi[j])<0.4)
	iso=iso-gpt[j];
    }
    if(iso<0)
      iso=0;
    result.emplace_back(iso/mpt[i]);
  }
  return result;
}




RVec<bool> fsr_recovery(const RVec<size_t>& idx,RVecF mpt, RVecF meta, RVecF mphi, RVecF mmass,RVecF gpt, RVecF geta, RVecF gphi,RVec<bool> gID) {
  RVec<bool> result;  
  ROOT::Math::PtEtaPhiMVector p1(mpt[idx[0]], meta[idx[0]], mphi[idx[0]], mmass[idx[0]]);
  ROOT::Math::PtEtaPhiMVector p2(mpt[idx[1]], meta[idx[1]], mphi[idx[1]], mmass[idx[1]]);
  ROOT::Math::PtEtaPhiMVector Z = p1+p2;
  float m = Z.M();

   //add the FSR
   for(size_t i=0;i<gpt.size();++i) {
     bool keep=false;
     if (gID[i]==1) {
       ROOT::Math::PtEtaPhiMVector gamma(gpt[i], geta[i], gphi[i],0);
       float m_fsr = (Z+gamma).M();
       if(std::abs(m_fsr-Z_mass)<std::abs(m-Z_mass)) {
	 keep=true;
	 m=m_fsr;
	 Z=Z+gamma;
       }
     }
     result.emplace_back(keep);
   }
  return result;
}





float ll_mass(const RVec<size_t>& idx,RVecF mpt, RVecF meta, RVecF mphi, RVecF mmass) {
  RVec<bool> result;  
  ROOT::Math::PtEtaPhiMVector p1(mpt[idx[0]], meta[idx[0]], mphi[idx[0]], mmass[idx[0]]);
  ROOT::Math::PtEtaPhiMVector p2(mpt[idx[1]], meta[idx[1]], mphi[idx[1]], mmass[idx[1]]);
  ROOT::Math::PtEtaPhiMVector Z = p1+p2;
  float m = Z.M();
  return m;
}

float calculate_llgamma_mass(const RVec<size_t>& idx,RVecF mpt, RVecF meta, RVecF mphi, RVecF mmass,RVecF gpt, RVecF geta, RVecF gphi, RVecF gFSR) {
  ROOT::Math::PtEtaPhiMVector p1(mpt[idx[0]], meta[idx[0]], mphi[idx[0]], mmass[idx[0]]);
  ROOT::Math::PtEtaPhiMVector p2(mpt[idx[1]], meta[idx[1]], mphi[idx[1]], mmass[idx[1]]);
  ROOT::Math::PtEtaPhiMVector Z = p1+p2;
  for (size_t i=0;i<gpt.size();++i) {
    ROOT::Math::PtEtaPhiMVector g(gpt[i], geta[i], gphi[i], 0);
    if(gFSR[i]) {
      Z=Z+g;
    }
  }
  return Z.M();
}

// Find index for scale factor given the value and the bin low edges
int getBin(const float val, const std::vector<float> bins){
  if (val <= bins.front())
    return 0;
  auto lower = std::lower_bound(bins.begin(), bins.end(), val); // Finds iterator of bin upper edge
  return (int) (std::distance(bins.begin(), lower) - 1);
}

// Get scale factors as a function of 2 variables
RVec<RVecF> scaleFactors_2d(RVecF x, RVecF y, std::vector <std::vector<std::vector<float>>> SFs, std::vector<float> binsX, std::vector<float> binsY, const bool isMC, RVecB selection){
  RVec<RVecF> out;
  RVecF vals;
  RVecF uncs;
  out.reserve(2);
  vals.reserve(x.size());
  uncs.reserve(x.size());
  for (size_t i = 0; i < x.size(); i++){
    if (!(selection[i]&&isMC)) {
      vals.emplace_back(1.0);
      uncs.emplace_back(0.0);
    }
    else {
      int binX = getBin(x[i], binsX);
      int binY = getBin(y[i], binsY);
      vals.emplace_back(SFs[binX][binY][0]);
      uncs.emplace_back(SFs[binX][binY][1]);
    }
  }
  out.emplace_back(vals);
  out.emplace_back(uncs);
  return out;
}

RVec<RVecF> scaleFactors_3d(RVecF x, RVecF y, RVecF z, std::vector <std::vector<std::vector<std::vector<float>>>> SFs, std::vector<float> binsX, std::vector<float> binsY, std::vector<float> binsZ, const bool isMC, RVecB selection){
  RVec<RVecF> out;
  RVecF vals;
  RVecF uncs;
  out.reserve(2);
  vals.reserve(x.size());
  uncs.reserve(x.size());
  for (size_t i = 0; i < x.size(); i++){
    if (!(selection[i]&&isMC)) {
      vals.emplace_back(1.0);
      uncs.emplace_back(0.0);
    }
    else {
      int binX = getBin(x[i], binsX);
      int binY = getBin(y[i], binsY);
      int binZ = getBin(z[i], binsZ);
      vals.emplace_back(SFs[binX][binY][binZ][0]);
      uncs.emplace_back(SFs[binX][binY][binZ][1]);
    } 
  }
  out.emplace_back(vals);
  out.emplace_back(uncs);
  return out;
}

// Get scale factors for ele reco (pt > 20, pt < 20 separate)
RVec<RVecF> scaleFactors_eleReco(RVecF eta, RVecF pt, std::vector <std::vector<std::vector<float>>> SFs_below, std::vector<float> binsX_below, std::vector<float> binsY_below,std::vector <std::vector<std::vector<float>>> SFs_above, std::vector<float> binsX_above, std::vector<float> binsY_above, const bool isMC, RVecB selection){
  RVec<RVecF> out;
  RVecF vals;
  RVecF uncs;
  out.reserve(2);
  vals.reserve(pt.size());
  uncs.reserve(pt.size());
  for (size_t i = 0; i < pt.size(); i++){
    if (!(selection[i]&&isMC)) {
      vals.emplace_back(1.0);
      uncs.emplace_back(0.0);
    }
    else {
      if (pt[i] < 20){
	int binX = getBin(eta[i], binsX_below);
	int binY = getBin(pt[i], binsY_below);
	vals.emplace_back(SFs_below[binX][binY][0]);
	uncs.emplace_back(SFs_below[binX][binY][1]);
      }
      else{
	int binX = getBin(eta[i], binsX_above);
	int binY = getBin(pt[i], binsY_above);
	vals.emplace_back(SFs_above[binX][binY][0]);
	uncs.emplace_back(SFs_above[binX][binY][1]);
      }
    }
  }
  out.emplace_back(vals);
  out.emplace_back(uncs);
  return out;
}

RVecI matchTrigger(RVecF eta, RVecF phi, RVecI pdgId, RVecF trig_eta, RVecF trig_phi, RVecF trig_pt, RVecI trig_id, RVecI trig_bits, const int matchingBits, const float ptThresh){
  RVecI out;
  out.reserve(eta.size());
  for (size_t i = 0; i < eta.size(); i++){
    bool matched = false;
    for (size_t j = 0; j < trig_eta.size(); j++){
      if (std::abs(pdgId[i]) != std::abs(trig_id[j]))
	continue;
      if (trig_pt[j] < ptThresh)
	continue;
      if (DeltaR(eta[i], trig_eta[j], phi[i], trig_phi[j]) < 0.1){
	if ((trig_bits[j] & matchingBits) == matchingBits){
	  matched = true;
	  break;
	}
      }
    }
    out.push_back(matched);
  }
  return out;
}


float getPUweight(const int truePU, std::vector<float> weights, const bool isMC){
  float out = 1.0;
  if (!isMC)
    return out;
  if (truePU > 99 || truePU < 0)
    return out;
  return weights[truePU];
}

RVec<RVecF> getPixelSeedSF(RVecB isEB, RVecB isEE, std::vector< std::vector<float>> weights, const bool isMC, RVecB selection){
  RVec<RVecF> out;
  RVecF vals;
  RVecF uncs;
  vals.reserve(isEB.size());
  uncs.reserve(isEB.size());
  for (size_t i = 0; i < isEB.size(); i++){
    if (!isMC){
      vals.emplace_back(1.0);
      uncs.emplace_back(0.0);
    }
    else if (!selection[i]){
      vals.emplace_back(1.0);
      uncs.emplace_back(0.0);
    }
    else if (isEB[i]){
      vals.emplace_back(weights[0][0]);
      uncs.emplace_back(weights[0][1]);
    }
    else if (isEE[i]){
      vals.emplace_back(weights[1][0]);
      uncs.emplace_back(weights[1][1]);
    }
    else{
      vals.emplace_back(1.0);
      uncs.emplace_back(0.0);
    }
  }
  out.emplace_back(vals);
  out.emplace_back(uncs);
  return out;
}

RVecF check_W_misID(const float e_pt, const float e_eta, const float e_phi, const float e_mass, RVecF g_pt, RVecF g_eta, RVecF g_phi, const int idx1, const int idx2){
  RVecF out;
  out.reserve(3);
  ROOT::Math::PtEtaPhiMVector e(e_pt, e_eta, e_phi, e_mass);
  ROOT::Math::PtEtaPhiMVector g1(g_pt[idx1], g_eta[idx1], g_phi[idx1], 0.0);
  ROOT::Math::PtEtaPhiMVector g2(g_pt[idx2], g_eta[idx2], g_phi[idx2], 0.0);
  out.emplace_back((e+g1).mass());
  out.emplace_back((e+g2).mass());
  out.emplace_back((e+g1+g2).mass());
  return out;
}

RVecF getctau(RVecF dx, RVecF dy, RVecF dz, RVecF pt, RVecF eta, RVecF phi, RVecF mass){
  RVecF out;
  out.reserve(pt.size());
  
  for (size_t i = 0; i < pt.size(); i++){
    float ip3d = std::sqrt(dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i]);
    ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], mass[i]);
    float ctau = ip3d/p4.Gamma()/p4.Beta();
    out.emplace_back(ctau);
  }

  return out;
}

RVecF getGenScEta(RVecF vx, RVecF vy, RVecF vz, RVecF pt, RVecF eta, RVecF phi, RVecF mass){
  RVecF out;
  out.reserve(vx.size());
  for (size_t i = 0; i < vx.size(); i++){
    float lxy = std::sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
    float theta = 2*std::atan(std::exp(-eta[i]));
    float ecal_z = (137.0 - lxy)/std::tan(theta) + vz[i];
    if (std::abs(ecal_z) < 300){
      float ecal_r = 137.;
      float ecal_theta = std::atan2(ecal_r, ecal_z);
      float ecal_eta = -std::log(std::tan(ecal_theta/2.));
      out.emplace_back(ecal_eta);
    }
    else{
      if (ecal_z > 0)
	ecal_z = 324;
      else
	ecal_z = -324;
      ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], mass[i]);
      float dx = p4.px()/p4.pz()*(ecal_z - vz[i]);
      float dy = p4.py()/p4.pz()*(ecal_z - vz[i]);
      float ecal_x = vx[i] + dx;
      float ecal_y = vy[i] + dy;
      float ecal_r = std::sqrt(ecal_x*ecal_x+ecal_y*ecal_y);
      float ecal_theta = std::atan2(ecal_r, ecal_z);
      float ecal_eta = -std::log(std::tan(ecal_theta/2.));
      out.emplace_back(ecal_eta);
    }
  }
  return out;
}

RVecF getGenScPhi(RVecF vx, RVecF vy, RVecF vz, RVecF pt, RVecF eta, RVecF phi, RVecF mass){
  RVecF out;
  out.reserve(vx.size());
  for (size_t i = 0; i < vx.size(); i++){
    float lxy = std::sqrt(vx[i]*vx[i]+vy[i]*vy[i]);
    float theta = 2*std::atan(std::exp(-eta[i]));
    ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], mass[i]);
    float ecal_z = (137. - lxy)/std::tan(theta) + vz[i];
    if (std::abs(ecal_z) > 300){
      if (ecal_z > 0)
	ecal_z = 310;
      else
	ecal_z = -310;
      float ecal_x = vx[i] + p4.px()/p4.pz()*(ecal_z-vz[i]);
      float ecal_y = vy[i] + p4.py()/p4.pz()*(ecal_z-vz[i]);
      float ecal_phi = std::atan2(ecal_y, ecal_x);
      out.emplace_back(ecal_phi);
    }
    else{
      float dx = p4.px()/pt[i];
      float dy = p4.py()/pt[i];
      float dr = std::sqrt(dx*dx+dy*dy);
      float D = vx[i]*(vy[i]+p4.py()/pt[i]) - (vx[i]+p4.px()/pt[i])*vy[i];
      float delta = (137*137.)*(dr*dr)-D*D;
      if (delta < 0){
	out.emplace_back(std::atan2(vy[i], vx[i]));
	continue;
      }
      float xP = (D*dy + std::copysign(1, dy)*dx*std::sqrt(delta))/(dr*dr);
      float xM = (D*dy - std::copysign(1, dy)*dx*std::sqrt(delta))/(dr*dr);
      float yP = (-D*dx + abs(dy)*std::sqrt(delta))/(dr*dr);
      float yM = (-D*dx - abs(dy)*std::sqrt(delta))/(dr*dr);
      float phi1 = std::atan2(yP, xP);
      float phi2 = std::atan2(yM, xM);
      float prod1 = p4.px()*(xP-vx[i])+p4.py()*(yP-vy[i]);
      float prod2 = p4.px()*(xM-vx[i])+p4.py()*(yM-vy[i]);
      if (prod1 > prod2)
	out.emplace_back(phi1);
      else
	out.emplace_back(phi2);
    }
  }
  return out;
}


RVecI isGenSignal(RVecI pdgId, RVecI motherIdx){
  RVecI out;
  out.reserve(pdgId.size());
  for (size_t i = 0; i < pdgId.size(); i++){
    if (motherIdx[i] == -1){
      out.emplace_back(0);
      continue;
    }
    if (pdgId[i] == 22 && std::abs(pdgId[motherIdx[i]]) == 9000006)
      out.emplace_back(1);
    else
      out.emplace_back(0);
  }
  return out;
}

RVecF photonEnergyScale(RVecF eta, RVecU seedgain, std::vector<std::pair<unsigned int, std::vector<float>>> values, std::vector<std::pair<unsigned int, std::vector<float>>> binning, const bool isMC){
  RVecF out;
  out.reserve(eta.size());
  for (size_t i = 0; i < eta.size(); i++){
    if (!isMC){
      out.emplace_back(1.0);
      continue;
    }
    unsigned int gain = seedgain[i];

    // Find iterator for values/bins corresponding to the seedGain
    auto it_vals = std::find_if(values.begin(), values.end(),
			   [gain](const std::pair<unsigned int, std::vector<float>>& element) { return element.first == gain;}
			   );
    auto it_bins = std::find_if(binning.begin(), binning.end(),
			   [gain](const std::pair<unsigned int, std::vector<float>>& element) { return element.first == gain;}
			   );
    if (it_vals == values.end() || it_bins == binning.end())
      out.emplace_back(1.0);
    else{
      int bin = getBin(eta[i], it_bins->second);
      out.emplace_back(it_vals->second[bin]);
    }
  }
  return out;
}


RVecF Photon_corrRelIso_EGM_loose_np(const RVecF& gpt, const RVecB& isEE, const RVecB& isEB, const RVecF& giso){
  RVecB result;
  result.reserve(gpt.size());

  for (size_t i = 0; i < gpt.size(); ++i) {
    const float p = gpt[i];
	const float iso = giso[i];

	float neutral;
	float photon;

    if (isEE[i]) {
      neutral = 19.722 + 0.0117*p + 2.3e-5*p*p;
      photon  = 4.162 + 0.0037*p;
    }
    else if (isEB[i]) {
      neutral = 24.032 + 0.01512*p + 2.259e-5*p*p;
      photon  = 2.876 + 0.004017*p;
    } 

	float iso_cut = min(photon/p, neutral/p);
	
	if (iso<iso_cut){
		result.emplace_back(true);
	} else{
		result.emplace_back(false);
	}
  }

  return result;
}


RVecF Photon_corrRelIso_EGM_loose_chg(const RVecB& isEE, const RVecB& isEB, const RVecF& giso){
  RVecB result;
  result.reserve(giso.size());

  for (size_t i = 0; i < giso.size(); ++i) {
	const float iso = giso[i];

	float charged;

    if (isEE[i]) {
		charged=2.089;
    }
    else if (isEB[i]) {
		charged=1.694;
    } 

	float iso_cut = charged;
	
	if (iso<iso_cut){
		result.emplace_back(true);
	} else{
		result.emplace_back(false);
	}
  }

  return result;
}


RVecF Photon_corrRelIso_custom_loose_np(const RVecF& gpt, const RVecB& isEE, const RVecB& isEB, const RVecF& giso){
  RVecB result;
  result.reserve(gpt.size());

  for (size_t i = 0; i < gpt.size(); ++i) {
    const float p = gpt[i];
	const float iso = giso[i];

	float neutral;
	float photon;

    if (isEE[i]) {
      neutral = 19.722 + 0.0117*p + 2.3e-5*p*p;
      photon  = 4.162 + 0.0037*p;
    }
    else if (isEB[i]) {
      neutral = 24.032 + 0.01512*p + 2.259e-5*p*p;
      photon  = 2.876 + 0.004017*p;
    } 

	float iso_cut = min(1.5*(photon/p), 1.5*(neutral/p));
	
	if (iso<(iso_cut)){
		result.emplace_back(true);
	} else{
		result.emplace_back(false);
	}
  }

  return result;
}


RVecF Photon_corrRelIso_custom_loose_chg(const RVecB& isEE, const RVecB& isEB, const RVecF& giso){
  RVecB result;
  result.reserve(giso.size());

  for (size_t i = 0; i < giso.size(); ++i) {
	const float iso = giso[i];

	float charged;

    if (isEE[i]) {
		charged=2.089;
    }
    else if (isEB[i]) {
		charged=1.694;
    } 

	float iso_cut = 1.5*charged;
	
	if (iso<iso_cut){
		result.emplace_back(true);
	} else{
		result.emplace_back(false);
	}
  }

  return result;
}


RVecF Photon_corrRelIso_EGM_loose_cut_np(const RVecF& gpt, const RVecB& isEE, const RVecB& isEB){
  RVecF result;
  result.reserve(gpt.size());

  for (size_t i = 0; i < gpt.size(); ++i) {
    const float p = gpt[i];

	float neutral;
	float photon;

    if (isEE[i]) {
      neutral = 19.722 + 0.0117*p + 2.3e-5*p*p;
      photon  = 4.162 + 0.0037*p;
    }
    else if (isEB[i]) {
      neutral = 24.032 + 0.01512*p + 2.259e-5*p*p;
      photon  = 2.876 + 0.004017*p;
    } 

	float iso_cut = min(photon/p, neutral/p);
	result.emplace_back(iso_cut);
	
  }

  return result;
}


RVecF Photon_corrRelIso_EGM_loose_cut_chg(const RVecB& isEE, const RVecB& isEB, const RVecF& giso){
  RVecF result;
  result.reserve(giso.size());

  for (size_t i = 0; i < giso.size(); ++i) {
	const float iso = giso[i];

	float charged;

    if (isEE[i]) {
		charged=2.089;
    }
    else if (isEB[i]) {
		charged=1.694;
    } 

	float iso_cut = charged;
	result.emplace_back(iso_cut);
  }

  return result;
}

RVecF Photon_corrRelIso_custom_loose_cut_np(const RVecF& gpt, const RVecB& isEE, const RVecB& isEB){
  RVecF result;
  result.reserve(gpt.size());

  for (size_t i = 0; i < gpt.size(); ++i) {
    const float p = gpt[i];

    float neutral;
    float photon;

    if (isEE[i]) {
      neutral = 19.722 + 0.0117*p + 2.3e-5*p*p;
      photon  = 4.162 + 0.0037*p;
    }
    else if (isEB[i]) {
      neutral = 24.032 + 0.01512*p + 2.259e-5*p*p;
      photon  = 2.876 + 0.004017*p;
    }

    float iso_cut = min(1.5*(photon/p), 1.5*(neutral/p));
    result.emplace_back(iso_cut);

  }

  return result;
}


RVecF Photon_corrRelIso_custom_loose_cut_chg(const RVecB& isEE, const RVecB& isEB, const RVecF& giso){
  RVecF result;
  result.reserve(giso.size());

  for (size_t i = 0; i < giso.size(); ++i) {
	const float iso = giso[i];

	float charged;

    if (isEE[i]) {
		charged=2.089;
    }
    else if (isEB[i]) {
		charged=1.694;
    } 

	float iso_cut = 1.5*charged;
	result.emplace_back(iso_cut);
  }

  return result;
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

#endif
