#include "math.h"
#include <Math/Vector4Dfwd.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <Math/GenVector/PxPyPzM4D.h>
#include "TVector3.h"
#include "TMath.h"
#include <algorithm>
#include <stdio.h>
class VertexCalculator {
 private:
  TVector3 vertex_;
  float pt_;
  float phi_;
  float d0_;
  float ip3d_;
  bool valid_;
  bool verbose_ = false;
  
 public:
  VertexCalculator() {
    vertex_= TVector3(0.0,0.0,0.0);
    pt_=0;
    phi_=0;
    d0_=0;
    ip3d_=0;
    valid_=false;

  }
  ~VertexCalculator() {} ;
  
  void setVertex(const TVector3& vertex){vertex_ = vertex;}
  void setPt(double pt){pt_ = pt;}  
  void setPhi(double phi){phi_ = phi;}
  void setD0(double d0){d0_ = d0;}
  void setIp3d(double ip3d){ip3d_ = ip3d;}
  void setValid(bool valid){valid_ = valid;}
  void setVerbose(bool verbose){verbose_ = verbose;}
  TVector3 vertex(){return vertex_;}    
  float pt(){return pt_;}
  float phi(){return phi_;}
  float d0(){return d0_;}
  float ip3d(){return ip3d_;}
  bool valid(){return valid_;}
  bool verbose(){return verbose_;}
    

    // Propagate prompt photon to ecal
TVector3 propToEcal(const double eta, const double phi, const int isEE, const int isEB){
  TVector3 out;
  float theta = 2*atan(exp(-eta));
  float ecal_z;
  if (theta == 0)
    ecal_z = 0;
  else
    ecal_z = 137./tan(theta);
  if (isEE){
    ecal_z = 324*ecal_z/abs(ecal_z);
    float r = ecal_z * tan(theta);
    float dx = cos(phi) * r;
    float dy = sin(phi) * r;
    out.SetXYZ(dx, dy, ecal_z);
  }
  else if (isEB){
    float r = 137.;
    float x = cos(phi) * r;
    float y = sin(phi) * r;
    out.SetXYZ(x, y, ecal_z);
  }
  else{ // Cases where both isEE an isEB==false
    if (abs(eta) < 1.45){ // Treat as barrel
      float r = 137;
      float x = cos(phi) * r;
      float y = sin(phi) * r;
      out.SetXYZ(x, y, ecal_z);
    }
    else{ // Treat as endcap
      ecal_z = 324*ecal_z/abs(ecal_z);
      float r = ecal_z * tan(theta);
      float dx = cos(phi) * r;
      float dy = sin(phi) * r;
      out.SetXYZ(dx, dy, ecal_z);
    }
  }
  return out;
}

float getRotAngle(const TVector3& v1, const TVector3& v2){
  TVector3 v3 = v1.Cross(v2);
  v3.SetMag(1.0);
  return acos(v3[2]);
}

TVector3 getRotAxis(const TVector3& v1, const TVector3& v2){
  if (v1[2]==0 && v2[2]==0)
    return TVector3(0,0,1.0);
  TVector3 v3 = v1.Cross(v2);
  v3.SetMag(1.0);
  TVector3 ez(0,0,1.0);
  return v3.Cross(ez);
}

double getPhi(const double mass, const double e1, const double e2){
  double cosphi = 1.0-mass*mass/(2.0*e1*e2);
  if (cosphi>=1 || cosphi<=-1)
    return M_PI;
  return acos(cosphi);
}

double getRadius(const TVector3& v1, const TVector3& v2, const double phi){
  double x1 = v1[0];
  double y1 = v1[1];
  double x2 = v2[0];
  double y2 = v2[1];
  return sqrt(pow((x1-x2)/2.,2)+pow((y1-y2)/2.,2))/sin(phi);
}

std::vector<std::pair<double, double>> getCenters(const TVector3& v1, const TVector3& v2, const double phi){
  double rad = getRadius(v1, v2, phi);
  double x1 = v1[0];
  double y1 = v1[1];
  double x2 = v2[0];
  double y2 = v2[1];
  double q = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  double x3 = (x1+x2)/2.;
  double y3 = (y1+y2)/2.;
  double xp = x3 + sqrt(rad*rad-(q*q/4.))*(y1-y2)/q;
  double yp = y3 + sqrt(rad*rad-(q*q/4.))*(x2-x1)/q;
  double xm = x3 - sqrt(rad*rad-(q*q/4.))*(y1-y2)/q;
  double ym = y3 - sqrt(rad*rad-(q*q/4.))*(x2-x1)/q;
  std::pair<double, double> plus, minus;
  plus.first = xp;
  plus.second = yp;
  minus.first = xm;
  minus.second = ym;
  std::vector<std::pair<double, double>> out;
  out.push_back(plus);
  out.push_back(minus);
  return out;
}

TVector3 findIntersection(const double theta, const double radius, const double x0, const double y0, const TVector3& v1, const TVector3& v2, const double phi){
  TVector3 out(0,0,0);
  double a = 1.0+tan(theta)*tan(theta);
  double b = -2.*(x0+tan(theta)*y0);
  double c = pow(x0, 2) + pow(y0, 2) - pow(radius, 2);
  if ((pow(b, 2) - 4*a*c) < 0)
    return out;
  double xp = (-b + sqrt(pow(b, 2) - 4*a*c))/(2*a);
  double xm = (-b - sqrt(pow(b, 2) - 4*a*c))/(2*a);
  double yp = xp * tan(theta);
  double ym = xm * tan(theta);
  
  // Pick point that gives correct angle
  TVector3 p(xp, yp, 0.0);
  TVector3 m(xm, ym, 0.0);

  TVector3 p1 = p - v1;
  TVector3 p2 = p - v2;
  TVector3 m1 = m - v1;
  TVector3 m2 = m - v2;

  double phi1 = acos(p1*p2/(p1.Mag()*p2.Mag()));
  double phi2 = acos(m1*m2/(m1.Mag()*m2.Mag()));

  double dphi1 = fabs(phi1 - phi);
  double dphi2 = fabs(phi2 - phi);

  if (dphi1 <= dphi2){
    out.SetXYZ(xp, yp, 0);
  }else{
    out.SetXYZ(xm, ym, 0);
  }
  return out;
}


double getPt(const TVector3& vertex, const TVector3& v1, const TVector3& v2, const double e1, const double e2){
  TVector3 temp1 = v1-vertex;
  TVector3 temp2 = v2-vertex;
  double sinTheta1 = (vertex.Cross(temp1).Mag())/(vertex.Mag()*temp1.Mag())*copysign(1.0, vertex.Cross(temp1)[2]);
  double sinTheta2 = (vertex.Cross(temp2).Mag())/(vertex.Mag()*temp2.Mag())*copysign(1.0, vertex.Cross(temp2)[2]);
  double pt1 = e1 * sinTheta1;
  double pt2 = e2 * sinTheta2;
  return pt1 + pt2;
  
}


bool tooFar(const TVector3& vertex, const TVector3& v1, const TVector3& v2){
  double x,y;
  // Special cases for vertex at (0, y) or ecal clusters at (x, y1) (x, y2)
  if (vertex[0] == 0){
    if (v1[0]==v2[0])
      return fmax(v1[1], v2[1]) > vertex[1] && fmin(v1[1], v2[1]) < vertex[1];
    double m = (v2[1]-v1[1])/(v2[0]-v1[0]);
    y = v1[1]-m*v1[0];
    x = 0;
  } else {
    double mv = vertex[1]/vertex[0];
    if (v1[0]==v2[0]){
      x = v1[0];
      y = mv*x;
    } else {
      double mp = (v2[1]-v1[1])/(v2[0]-v1[0]);
      x = (v1[1]-mp*v1[0])/(mv-mp);
      y = mv*x;
    }
  }
  return (fmax(0, vertex[0]) > x && fmin(0, vertex[0]) < x && fmax(0, vertex[1]) > y && fmin(0, vertex[1]) < y);
}

int backwards(const TVector3& vertex, const TVector3& v1, const TVector3& v2){
  int out = 1;
  TVector3 midpoint;
  //  Find closest point to line connecting two ecal clusters
  if (v1[0]==v2[0])
    midpoint.SetXYZ(v1[0], 0, 0);
  else{
    double m = (v2[1]-v1[1])/(v2[0]-v1[0]);
    midpoint.SetX((v1[1]-m*v1[0])/(-1./m-m));
    midpoint.SetY(-1./m*(v1[1]-m*v1[0])/(-1./m-m));
  }

  if (midpoint[0]*vertex[0]+midpoint[1]*vertex[1] < 0)
    out = -1;
  if (verbose()){
    printf("Checking backwards: midpoint=(%f, %f) vertex=(%f,%f) \n", midpoint[0], midpoint[1],vertex[0],vertex[1]);
    printf("Midpoint * vertex = %f\n", midpoint*vertex);
    printf("scale=%d\n", out);
  }
  return out;
}

// Set the vertex, pt, and phi
void run(const TVector3& v1, const TVector3& v2, const double e1, const double e2, const double mass){


  double phi = M_PI;
  double pt = 999;
  TVector3 vertex(-999,-999,-999);

  setVertex(vertex);
  setPhi(phi);
  setPt(pt);
  setD0(-sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]));
  setIp3d(-sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]+vertex[2]*vertex[2]));
  setValid(false);

  TVector3 axis = getRotAxis(v1, v2);
  double theta = getRotAngle(v1, v2);
  phi = getPhi(mass, e1, e2);
  if (phi >= M_PI){
    if (verbose_)
      printf("Opening angle > pi, return false\n");
    return;
  }

  TVector3 temp1 = v1;
  TVector3 temp2 = v2;
  temp1.Rotate(theta, axis);
  temp2.Rotate(theta, axis);
  double radius = getRadius(temp1, temp2, phi);
  std::vector<std::pair<double, double>> centers = getCenters(temp1, temp2, phi);
  if (verbose()){
    printf("Ecal cluster 1 energy=%f (%f, %f)\nEcal cluster 2 energy=%f (%f, %f)\n",e1, temp1[0], temp1[1], e2, temp2[0], temp2[1]);
    printf("Photon phi=%f radius=%f\n", phi, radius);
  }

  for(uint i = 0; i < centers.size(); i++){
    pt = 999;
    std::pair<double, double> center = centers[i];
    if (verbose()){
      printf("   Checking circle with center (%f, %f)\n", center.first, center.second);
    }
    
    //Find ranges of phi values to check on locus of points
    double phi1 = atan2(temp1[1], temp1[0]);
    if (phi1 < 0)
      phi1 = M_PI*2 + phi1;
    double phi2 = atan2(temp2[1], temp2[0]);
    if (phi2 < 0)
      phi2 = M_PI*2 + phi2;
    double deltaPhi = phi2-phi1;
    //Bound deltaPhi between -pi,pi
    if (deltaPhi > M_PI){
      deltaPhi = M_PI*2 - deltaPhi;
      phi1 = atan2(temp2[1], temp2[0]);
      phi2 = atan2(temp1[1], temp1[0]);
    }
    if (deltaPhi < -M_PI){
      deltaPhi = M_PI*2 + deltaPhi;
      phi1 = atan2(temp1[1], temp1[0]);
      phi2 = atan2(temp2[1], temp2[0]);
    }
    
    int stepsPhi = 100;
    double stepSize = deltaPhi/stepsPhi;   
    for (int i=0; i<stepsPhi; i++){
      double angle = phi1 + i*stepSize;
      TVector3 temp = findIntersection(angle, radius, center.first, center.second, temp1, temp2, phi);
      double tempPt = getPt(temp, temp1, temp2, e1, e2);
      if (abs(tempPt) <= abs(pt)){
	pt = tempPt;
	vertex = temp;
      }
    }
    
    if (verbose())
      printf ("   Best fit vertex (%f, %f), pt = %f\n", vertex[0], vertex[1], pt);

    // Scale factor <0 means vertex is behind beamline
    int scale = backwards(vertex, temp1, temp2);
    bool unphysical = tooFar(vertex, temp1, temp2);
    if (verbose())
      printf("   checks: scale=%d toofar=%d\n", scale, unphysical);
    // Always pick physical vertex, otherwise only pick backwards if no other vertex found
    // The two circles cannot produce BOTH a physical and a backwards vertex
    if (!unphysical){
      if (verbose())
	printf("   Setting (%f, %f) as vertex\n", vertex[0], vertex[1]);
      vertex.Rotate(-theta, axis);
      if (verbose())
	printf("   Rotated vector around axis (%f, %f, %f) by angle %f\n", axis[0], axis[1], axis[2], -theta);
      setVertex(vertex);
      setPt(pt);
      setPhi(phi);
      setD0(scale*sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]));
      setIp3d(scale*sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]+vertex[2]*vertex[2]));
      setValid(true);
      if (verbose())
	printf("   New Vertex set at (%f, %f, %f), d0=%f, phi=%f\n", vertex[0], vertex[1], vertex[2], d0(), phi);
    }
  }
  return;
}

// With propagation, get the vertex within 3% of using precision calo clusters
void runPtEtaPhi(const double pt1, const double eta1, const double phi1, const int isEE1, const int isEB1, const double pt2, const double eta2, const double phi2, const int isEE2, const int isEB2, const double mass){
    TVector3 v1 = propToEcal(eta1, phi1, isEE1, isEB1);
    TVector3 v2 = propToEcal(eta2, phi2, isEE2, isEB2);
    ROOT::Math::PtEtaPhiMVector p4_1(pt1, eta1, phi1, 0);
    ROOT::Math::PtEtaPhiMVector p4_2(pt2, eta2, phi2, 0);
    run(v1, v2, p4_1.energy(), p4_2.energy(), mass);
    return;
}


// Call only after running to set p4 using vertex coords
ROOT::Math::PxPyPzMVector correctP4(const double pt, const double eta, const double phi, const int isEE, const int isEB){
  TVector3 v = vertex();
  TVector3 calo = propToEcal(eta, phi, isEE, isEB);
  TVector3 dR = calo - v;
  
  dR.SetMag(1.0);
  ROOT::Math::PtEtaPhiMVector p4(pt, eta, phi, 0);
  float e = p4.energy();
  ROOT::Math::PxPyPzMVector out(e*dR.x(), e*dR.y(), e*dR.z(), 0);
  return out;   
}

// Returns vector of [pt1, eta1, phi1, pt2, eta2, phi2, lxy, valid]
std::vector<float> getVertexInfo(const double pt1, const double eta1, const double phi1, const int isEE1, const int isEB1, const double pt2, const double eta2, const double phi2, const int isEE2, const int isEB2, const double mass){
  std::vector<float> out;
  runPtEtaPhi(pt1, eta1, phi1, isEE1, isEB1, pt2, eta2, phi2, isEE2, isEB2, mass);
  ROOT::Math::PxPyPzMVector p4_1 = correctP4(pt1, eta1, phi1, isEE1, isEB1);
  ROOT::Math::PxPyPzMVector p4_2 = correctP4(pt2, eta2, phi2, isEE2, isEB2);
  out.push_back(p4_1.pt());
  out.push_back(p4_1.eta());
  out.push_back(p4_1.phi());
  out.push_back(p4_2.pt());
  out.push_back(p4_2.eta());
  out.push_back(p4_2.phi());
  out.push_back(d0());
  out.push_back(valid());
  return out;
}

};


bool compare_pair(const RVecF& q1,const RVecF& q2) {
  //first invariant mass of pairs
  unsigned int mass1 = q1[8]<65. ? 1:0;
  unsigned int mass2 = q2[8]<65 ? 1:0;

  unsigned int valid1 = q1[7]>0 ? 1:0;
  unsigned int valid2 = q2[7]>0 ? 1:0;

  unsigned int dxy1   = q1[6]>-8 ? 1:0;
  unsigned int dxy2   = q2[6]>-8 ? 1:0;
  ROOT::Math::PtEtaPhiMVector p1_0(q1[0],q1[1],q1[2],0.0);
  ROOT::Math::PtEtaPhiMVector p1_1(q1[3],q1[4],q1[5],0.0);
  ROOT::Math::PtEtaPhiMVector p1 = p1_0+p1_1;

  ROOT::Math::PtEtaPhiMVector p2_0(q2[0],q2[1],q2[2],0.0);
  ROOT::Math::PtEtaPhiMVector p2_1(q2[3],q2[4],q2[5],0.0);
  ROOT::Math::PtEtaPhiMVector p2 = p2_0+p2_1;
  
  if (valid1>valid2)
    return true;
  else if (valid1<valid2)
    return false;
  else {
    if (mass1>mass2)
      return true;
    else if (mass1<mass2)
      return false;
    else {
      if (valid1==1) {
	if(dxy1>dxy2)
	  return true;
	else if (dxy1<dxy2)
	  return false;
	else {
	  if(p1.pt()>p2.pt())
	    return true;
	  else
	    return false;
	}
      }
      else {
	if(p1.pt()>p2.pt())
	  return true;
	else
	  return false;
      }
    }
  }



}

bool compare_pair_withId(const RVecF& q1,const RVecF& q2) {
  //first invariant mass of pairs
  unsigned int mass1 = q1[8]<65. ? 1:0;
  unsigned int mass2 = q2[8]<65 ? 1:0;

  unsigned int valid1 = q1[7]>0 ? 1:0;
  unsigned int valid2 = q2[7]>0 ? 1:0;

  unsigned int dxy1   = q1[6]>-8 ? 1:0;
  unsigned int dxy2   = q2[6]>-8 ? 1:0;

  unsigned int id1 = q1[13] + q1[14];
  unsigned int id2 = q2[13] + q2[14];

  ROOT::Math::PtEtaPhiMVector p1_0(q1[0],q1[1],q1[2],0.0);
  ROOT::Math::PtEtaPhiMVector p1_1(q1[3],q1[4],q1[5],0.0);
  ROOT::Math::PtEtaPhiMVector p1 = p1_0+p1_1;

  ROOT::Math::PtEtaPhiMVector p2_0(q2[0],q2[1],q2[2],0.0);
  ROOT::Math::PtEtaPhiMVector p2_1(q2[3],q2[4],q2[5],0.0);
  ROOT::Math::PtEtaPhiMVector p2 = p2_0+p2_1;
  
  if (valid1>valid2)
    return true;
  else if (valid1<valid2)
    return false;
  else {
    if (id1 > id2)
      return true;
    else if (id2 > id1)
      return false;
    else {
      if (mass1>mass2)
	return true;
      else if (mass1<mass2)
	return false;
      else {
	if (valid1==1) {
	  if(dxy1>dxy2)
	    return true;
	  else if (dxy1<dxy2)
	    return false;
	  else {
	    if(p1.pt()>p2.pt())
	      return true;
	    else
	      return false;
	  }
	}
	else {
	  if(p1.pt()>p2.pt())
	    return true;
	  else
	    return false;
	}
      }
    }
  }
}

RVecF best_2gamma(RVecF pt,RVecF eta, RVecF phi,RVec<bool> EB, RVec<bool> EE, RVecI isLoose, RVecI gID, float mass) {
  RVec<RVecF> all_combos;
  auto idx_cmb = ROOT::VecOps::Combinations(pt, 2);
  for (size_t i = 0; i < idx_cmb[0].size(); i++) {
    const auto i1 = idx_cmb[0][i];
    const auto i2 = idx_cmb[1][i];
    if ( (isLoose[i1] != 1) || (isLoose[i2] != 1) )
      continue;

    RVecF result;
    result.reserve(18);
    VertexCalculator *calc = new VertexCalculator();
    //four vector
    ROOT::Math::PtEtaPhiMVector p0(pt[i1],eta[i1],phi[i1],0.0);
    ROOT::Math::PtEtaPhiMVector p1(pt[i2],eta[i2],phi[i2],0.0);
    float raw_m = (p0+p1).M();
    std::vector<float> kin_fit= calc->getVertexInfo(pt[i1],eta[i1],phi[i1],EE[i1],EB[i1],pt[i2],eta[i2],phi[i2],EE[i2],EB[i2],mass); 
    delete calc;
    result.emplace_back(kin_fit[0]);
    result.emplace_back(kin_fit[1]);
    result.emplace_back(kin_fit[2]);
    result.emplace_back(kin_fit[3]);
    result.emplace_back(kin_fit[4]);
    result.emplace_back(kin_fit[5]);
    result.emplace_back(kin_fit[6]);
    result.emplace_back(kin_fit[7]);
    result.emplace_back(raw_m);
    result.emplace_back(i1);
    result.emplace_back(i2);
    result.emplace_back(DeltaPhi(phi[i1], phi[i2]));
    result.emplace_back(DeltaR(eta[i1], eta[i2], phi[i1], phi[i2]));
    result.emplace_back(gID[i1]>0); //photon cut based ID
    result.emplace_back(gID[i2]>0);
    result.emplace_back((p0+p1).pt());
    result.emplace_back((p0+p1).eta());
    result.emplace_back((p0+p1).phi());
    all_combos.emplace_back(result);
  }
  auto sortedIndices = ROOT::VecOps::Argsort(all_combos,compare_pair_withId);
  return all_combos[sortedIndices[0]];
}

RVecF best_3gamma(RVecF pt,RVecF eta, RVecF phi,RVec<bool> EB, RVec<bool> EE,RVecI isLoose ,RVecI gID, RVecF gIso,float mass) {  
  RVec<RVecF> all_combos;
  auto idx_cmb = ROOT::VecOps::Combinations(pt, 3);
  for (size_t i = 0; i < idx_cmb[0].size(); i++) {
    const auto i1 = idx_cmb[0][i];
    const auto i2 = idx_cmb[1][i];
    const auto i3 = idx_cmb[2][i];
      
    if (!(isLoose[i1] && isLoose[i2] && isLoose[i3])) // only use preselected photons
      continue;

    RVecF result;
    result.reserve(21);
    VertexCalculator *calc = new VertexCalculator();
    std::vector<float> best_pair;
    best_pair.reserve(14);
    float raw_m;
    int g1;
    int g2;
    int g3;
    //four vector
    ROOT::Math::PtEtaPhiMVector p0(pt[i1],eta[i1],phi[i1],0.0);
    ROOT::Math::PtEtaPhiMVector p1(pt[i2],eta[i2],phi[i2],0.0);
    ROOT::Math::PtEtaPhiMVector p2(pt[i3],eta[i3],phi[i3],0.0);

    //First Pair
    std::vector<float> pair_1 = calc->getVertexInfo(pt[i1],eta[i1],phi[i1],EE[i1],EB[i1],pt[i2],eta[i2],phi[i2],EE[i2],EB[i2],mass);

    //second pair
    std::vector<float> pair_2 = calc->getVertexInfo(pt[i1],eta[i1],phi[i1],EE[i1],EB[i1],pt[i3],eta[i3],phi[i3],EE[i3],EB[i3],mass);

    //third pair
    std::vector<float> pair_3 = calc->getVertexInfo(pt[i2],eta[i2],phi[i2],EE[i2],EB[i2],pt[i3],eta[i3],phi[i3],EE[i3],EB[i3],mass);


    if (compare_pair(pair_1,pair_2)) {
      if (compare_pair(pair_1,pair_3)){
        best_pair.insert(best_pair.end(),pair_1.begin(),pair_1.end());
        g1 = i1;
        g2 = i2;
        g3 = i3;
        raw_m = (p0+p1).M();

        

      } else {
        best_pair.insert(best_pair.end(),pair_3.begin(),pair_3.end());
        g1 = i2;
        g2 = i3;
        g3 = i1;
        raw_m = (p1+p2).M();
        
      }
      
    } else if (compare_pair(pair_2,pair_3)) {
      best_pair.insert(best_pair.end(),pair_2.begin(),pair_2.end());
      g1 = i1;
      g2 = i3;
      g3 = i2;
      raw_m = (p0+p2).M();
      

      } else {
        best_pair.insert(best_pair.end(),pair_3.begin(),pair_3.end());
        g1 = i2;
        g2 = i3;
        g3 = i1;
        raw_m = (p1+p2).M();
        

      }

    ROOT::Math::PtEtaPhiMVector pc0(best_pair[0],best_pair[1],best_pair[2],0.0);
    ROOT::Math::PtEtaPhiMVector pc1(best_pair[3],best_pair[4],best_pair[5],0.0);
    ROOT::Math::PtEtaPhiMVector pc2(pt[g3],eta[g3],phi[g3],0.0);

    delete calc;   

    result.emplace_back(best_pair[0]);
    result.emplace_back(best_pair[1]);
    result.emplace_back(best_pair[2]);
    result.emplace_back(best_pair[3]);
    result.emplace_back(best_pair[4]);
    result.emplace_back(best_pair[5]);
    result.emplace_back(best_pair[6]);
    result.emplace_back(best_pair[7]);
    result.emplace_back(raw_m);
    result.emplace_back(gID[g1] && gIso[g1]<0.1);
    result.emplace_back(gID[g2] && gIso[g2]<0.1);
    result.emplace_back(gID[g3] && gIso[g3]<0.1);
    result.emplace_back(pt[g3]);
    result.emplace_back(eta[g3]);
    result.emplace_back(phi[g3]);
    result.emplace_back((p0+p1+p2).M());
    result.emplace_back((pc0+pc1+pc2).M());
    result.emplace_back(DeltaPhi(phi[g1],phi[g2]));
    result.emplace_back(DeltaR(eta[g1],eta[g2],phi[g1],phi[g2]));
	result.emplace_back(g1);
    result.emplace_back(g2);
    result.emplace_back(g3);
    all_combos.emplace_back(result);
  }
  auto sortedIndices = ROOT::VecOps::Argsort(all_combos,compare_pair);
  return all_combos[sortedIndices[0]];
}

bool compare_quad(const RVecF& q1,const RVecF& q2) {
      ROOT::Math::PtEtaPhiMVector p1_A1(q1[0],q1[1],q1[2],0.0);
      ROOT::Math::PtEtaPhiMVector p1_A2(q1[3],q1[4],q1[5],0.0);
      ROOT::Math::PtEtaPhiMVector p1_B1(q1[8],q1[9],q1[10],0.0);
      ROOT::Math::PtEtaPhiMVector p1_B2(q1[11],q1[12],q1[13],0.0);
      float sumpt1 = (p1_A1+p1_A2).pt()+(p1_B1+p1_B2).pt();

      ROOT::Math::PtEtaPhiMVector p2_A1(q2[0],q2[1],q2[2],0.0);
      ROOT::Math::PtEtaPhiMVector p2_A2(q2[3],q2[4],q2[5],0.0);
      ROOT::Math::PtEtaPhiMVector p2_B1(q2[8],q2[9],q2[10],0.0);
      ROOT::Math::PtEtaPhiMVector p2_B2(q2[11],q2[12],q2[13],0.0);
      float sumpt2 = (p2_A1+p2_A2).pt()+(p2_B1+p2_B2).pt();


  //first invariant mass of pairs
  unsigned int mass1_1 = q1[16]<62.5 ? 1:0;
  unsigned int mass1_2 = q1[17]<62.5 ? 1:0;
  unsigned int mass2_1 = q2[16]<62.5 ? 1:0;
  unsigned int mass2_2 = q2[17]<62.5 ? 1:0;

  unsigned int valid1_1 = q1[7]>0 ? 1:0;
  unsigned int valid1_2 = q1[15]>0 ? 1:0;
  unsigned int valid2_1 = q2[7]>0 ? 1:0;
  unsigned int valid2_2 = q2[15]>0 ? 1:0;

  unsigned int dxy1_1   = q1[6]>-8 ? 1:0;
  unsigned int dxy1_2   = q1[14]>-8 ? 1:0;
  unsigned int dxy2_1   = q2[6]>-8 ? 1:0;
  unsigned int dxy2_2    = q2[14]>-8 ? 1:0;



  
  if ((mass1_1+mass1_2)>(mass2_1+mass2_2)) {
    return true;
  }
  else if ((mass1_1+mass1_2)<(mass2_1+mass2_2)) {
    return false;
  }
  else {
    if ((valid1_1+valid1_2)>(valid2_1+valid2_2))
      return true;
    else if ((valid1_1+valid1_2)<(valid2_1+valid2_2))
      return false;
    else {
      //if both valid look at dxy
      if ((valid1_1+valid1_2)==2) {
	if ((dxy1_1+dxy1_2)>(dxy2_1+dxy2_2))
	  return true;
	else if ((dxy1_1+dxy1_2)<(dxy2_1+dxy2_2))
	  return false;
	else {
	  if (sumpt1>sumpt2)
	    return true;
	  else
	    return false;
	}
      }
      else  {
	  if (sumpt1>sumpt2)
	    return true;
	  else
	    return false;
      }
    }
  }

}

bool compare_quad_pairing(const std::vector<float>p1_A,const std::vector<float>p1_B,const std::vector<float>p2_A,const std::vector<float>p2_B) {
      ROOT::Math::PtEtaPhiMVector p1_A1(p1_A[0],p1_A[1],p1_A[2],0.0);
      ROOT::Math::PtEtaPhiMVector p1_A2(p1_A[3],p1_A[4],p1_A[5],0.0);
      ROOT::Math::PtEtaPhiMVector p1_B1(p1_B[0],p1_B[1],p1_B[2],0.0);
      ROOT::Math::PtEtaPhiMVector p1_B2(p1_B[3],p1_B[4],p1_B[5],0.0);

      ROOT::Math::PtEtaPhiMVector p2_A1(p2_A[0],p2_A[1],p2_A[2],0.0);
      ROOT::Math::PtEtaPhiMVector p2_A2(p2_A[3],p2_A[4],p2_A[5],0.0);
      ROOT::Math::PtEtaPhiMVector p2_B1(p2_B[0],p2_B[1],p2_B[2],0.0);
      ROOT::Math::PtEtaPhiMVector p2_B2(p2_B[3],p2_B[4],p2_B[5],0.0);
      //first we see if the mass is everywhere below 125/2.
      unsigned int p1_masscut_A = (p1_A1+p1_A2).M()<65. ? 1:0; 
      unsigned int p1_masscut_B = (p1_B1+p1_B2).M()<65. ? 1:0; 
      unsigned int p2_masscut_A = (p2_A1+p2_A2).M()<65. ? 1:0; 
      unsigned int p2_masscut_B = (p2_B1+p2_B2).M()<65. ? 1:0; 

      if ((p1_masscut_A+p1_masscut_B)>(p2_masscut_A+p2_masscut_B)) {
	return true;
      }
      else if ((p1_masscut_A+p1_masscut_B)<(p2_masscut_A+p2_masscut_B)) {
	return false;
      } 
      else {

	if((p1_A[7]+p1_B[7])>(p2_A[7]+p2_B[7])) {
	  return true;
	}
	else if ((p1_A[7]+p1_B[7])<(p2_A[7]+p2_B[7])) {
	  return false; //should this be false? - changed to false by delano aug. 8/29/25
	} 
	else if ((p1_A[7]+p1_B[7])==(p2_A[7]+p2_B[7])) {
	  float p1_pt = (p1_A1+p1_A2).pt()+(p1_B1+p1_B2).pt();
	  float p2_pt = (p2_A1+p2_A2).pt()+(p2_B1+p2_B2).pt();

	  
	  //if both invalid or half valid just group by higher SUM pt
	  if ((p1_A[7]+p1_B[7])<2) {
	    if (p1_pt>=p2_pt)
	      return true;
	    else
	      return false;
	  }
	  else {
	    
	    //both valid if one of them is negative dxy pick the other
	    unsigned int p1_A_dxy =p1_A[6]>-8.0 ? 1 :0; 
	    unsigned int p1_B_dxy =p1_B[6]>-8.0 ? 1 :0; 
	    unsigned int p2_A_dxy =p2_A[6]>-8.0 ? 1 :0; 
	    unsigned int p2_B_dxy =p2_B[6]>-8.0 ? 1 :0; 
	    if((p1_A_dxy+p1_B_dxy)>(p2_A_dxy+p2_B_dxy))
	      return true;
	    else if ((p1_A_dxy+p1_B_dxy)<(p2_A_dxy+p2_B_dxy))
	      return false;
	    else {
	      //sum pt
	      if (p1_pt>=p2_pt)
		return true;
	      else
		return false;
	    }
	  }
	} 
      }
      return true;
} 


RVecF best_4gamma(RVecF pt,RVecF eta, RVecF phi,RVec<bool> EB, RVec<bool> EE,RVecI isLoose, RVecI gID, RVecF gIso, float mass) {

  RVec<RVecF> all_combos;
  RVecF result;
  result.reserve(28);
  auto idx_cmb = ROOT::VecOps::Combinations(pt, 4);
  VertexCalculator *calc = new VertexCalculator();
  std::vector<float> best23_A;
  std::vector<float> best23_B;
  best23_A.reserve(20);
  best23_B.reserve(20);
  std::vector<float> best_A;
  std::vector<float> best_B;
  best_A.reserve(20);
  best_B.reserve(20);
  int g1;
  int g2;
  int g3;
  int g4;

  for (size_t i = 0; i < idx_cmb[0].size(); i++) {
    const auto i1 = idx_cmb[0][i];
    const auto i2 = idx_cmb[1][i];
    const auto i3 = idx_cmb[2][i];
    const auto i4 = idx_cmb[3][i];

    if (!(isLoose[i1] && isLoose[i2] && isLoose[i3] && isLoose[i4])) // only use preselected photons
      continue;

    result.clear();
    best23_A.clear();
    best23_B.clear();
    best_A.clear();
    best_B.clear();

    //four vector
    ROOT::Math::PtEtaPhiMVector p0(pt[i1],eta[i1],phi[i1],0.0);
    ROOT::Math::PtEtaPhiMVector p1(pt[i2],eta[i2],phi[i2],0.0);
    ROOT::Math::PtEtaPhiMVector p2(pt[i3],eta[i3],phi[i3],0.0);
    ROOT::Math::PtEtaPhiMVector p3(pt[i4],eta[i4],phi[i4],0.0);
    float raw_m = (p0+p1+p2+p3).M();
    ROOT::Math::PtEtaPhiMVector pA = p0+p1;
    ROOT::Math::PtEtaPhiMVector pB = p2+p3;
  

    //first pairing
    std::vector<float> pairing1_A = calc->getVertexInfo(pt[i1],eta[i1],phi[i1],EE[i1],EB[i1],pt[i2],eta[i2],phi[i2],EE[i2],EB[i2],mass); 
    std::vector<float> pairing1_B = calc->getVertexInfo(pt[i3],eta[i3],phi[i3],EE[i3],EB[i3],pt[i4],eta[i4],phi[i4],EE[i4],EB[i4],mass); 
    
    //second pairing
    std::vector<float> pairing2_A = calc->getVertexInfo(pt[i1],eta[i1],phi[i1],EE[i1],EB[i1],pt[i3],eta[i3],phi[i3],EE[i3],EB[i3],mass); 
    std::vector<float> pairing2_B = calc->getVertexInfo(pt[i2],eta[i2],phi[i2],EE[i2],EB[i2],pt[i4],eta[i4],phi[i4],EE[i4],EB[i4],mass); 
    
    //third pairing
    std::vector<float> pairing3_A = calc->getVertexInfo(pt[i1],eta[i1],phi[i1],EE[i1],EB[i1],pt[i4],eta[i4],phi[i4],EE[i4],EB[i4],mass); 
    std::vector<float> pairing3_B = calc->getVertexInfo(pt[i2],eta[i2],phi[i2],EE[i2],EB[i2],pt[i3],eta[i3],phi[i3],EE[i3],EB[i3],mass); 
    
    
    
    if (compare_quad_pairing(pairing2_A,pairing2_B,pairing3_A,pairing3_B)) {
      best23_A.insert(best23_A.end(),pairing2_A.begin(),pairing2_A.end());
      best23_B.insert(best23_B.end(),pairing2_B.begin(),pairing2_B.end());
      g1 = i1;
      g2 = i3;
      g3 = i2;
      g4 = i4;
      pA = p0+p2;
      pB = p1+p3;
      
    }
    else {
      best23_A.insert(best23_A.end(),pairing3_A.begin(),pairing3_A.end());
      best23_B.insert(best23_B.end(),pairing3_B.begin(),pairing3_B.end());
      g1 = i1;
      g2 = i4;
      g3 = i2;
      g4 = i3;
      pA=p0+p3;
      pB=p1+p2;
    }
    if (compare_quad_pairing(pairing1_A,pairing1_B,best23_A,best23_B)) {
      best_A.insert(best_A.end(),pairing1_A.begin(),pairing1_A.end());
      best_B.insert(best_B.end(),pairing1_B.begin(),pairing1_B.end());
      g1 = i1;
      g2 = i2;
      g3 = i3;
      g4 = i4;
      pA=p0+p1;
      pB=p2+p3;
    }
    else {
      best_A.insert(best_A.end(),best23_A.begin(),best23_A.end());
      best_B.insert(best_B.end(),best23_B.begin(),best23_B.end());
    }

    result.emplace_back(best_A[0]);
    result.emplace_back(best_A[1]);
    result.emplace_back(best_A[2]);
    result.emplace_back(best_A[3]);
    result.emplace_back(best_A[4]);
    result.emplace_back(best_A[5]);
    result.emplace_back(best_A[6]);
    result.emplace_back(best_A[7]);
    result.emplace_back(best_B[0]);
    result.emplace_back(best_B[1]);
    result.emplace_back(best_B[2]);
    result.emplace_back(best_B[3]);
    result.emplace_back(best_B[4]);
    result.emplace_back(best_B[5]);
    result.emplace_back(best_B[6]);
    result.emplace_back(best_B[7]);
    //raw masses
    result.emplace_back(pA.M());
    result.emplace_back(pB.M());
    result.emplace_back(raw_m);
    //corrected mass
    ROOT::Math::PtEtaPhiMVector pc0(best_A[0],best_A[1],best_A[2],0.0);
    ROOT::Math::PtEtaPhiMVector pc1(best_A[3],best_A[4],best_A[5],0.0);
    ROOT::Math::PtEtaPhiMVector pc2(best_B[0],best_B[1],best_B[2],0.0);
    ROOT::Math::PtEtaPhiMVector pc3(best_B[3],best_B[4],best_B[5],0.0);
    result.emplace_back((pc0+pc1+pc2+pc3).M());
    //indices of used photons
    result.emplace_back(gID[g1] && gIso[g1]<0.1);
    result.emplace_back(gID[g2] && gIso[g2]<0.1);
    result.emplace_back(gID[g3] && gIso[g3]<0.1);
    result.emplace_back(gID[g4] && gIso[g4]<0.1);
	result.emplace_back(g1);
    result.emplace_back(g2);
    result.emplace_back(g3);
    result.emplace_back(g4);
    all_combos.emplace_back(result);
  }
  delete calc;  
  if (all_combos.size()>1) {
    auto sortedIndices = ROOT::VecOps::Argsort(all_combos,compare_quad);
    return all_combos[sortedIndices[0]];
  }
  else {
    return all_combos[0];
  }


}


RVecF best_4gamma_1bad(RVecF pt,RVecF eta, RVecF phi,RVec<bool> EB, RVec<bool> EE,RVecI isLoose, RVecI gID, RVecF gIso, float mass) {

  RVec<RVecF> all_combos;
  RVecF result;
  result.reserve(28);
  auto idx_cmb = ROOT::VecOps::Combinations(pt, 4);
  VertexCalculator *calc = new VertexCalculator();
  std::vector<float> best23_A;
  std::vector<float> best23_B;
  best23_A.reserve(20);
  best23_B.reserve(20);
  std::vector<float> best_A;
  std::vector<float> best_B;
  best_A.reserve(20);
  best_B.reserve(20);
  int g1;
  int g2;
  int g3;
  int g4;

  for (size_t i = 0; i < idx_cmb[0].size(); i++) {
    const auto i1 = idx_cmb[0][i];
    const auto i2 = idx_cmb[1][i];
    const auto i3 = idx_cmb[2][i];
    const auto i4 = idx_cmb[3][i];

    int loose_number = isLoose[i1]+isLoose[i2]+isLoose[i3]+isLoose[i4];
    if (loose_number != 3)
      continue;


    result.clear();
    best23_A.clear();
    best23_B.clear();
    best_A.clear();
    best_B.clear();

    //four vector
    ROOT::Math::PtEtaPhiMVector p0(pt[i1],eta[i1],phi[i1],0.0);
    ROOT::Math::PtEtaPhiMVector p1(pt[i2],eta[i2],phi[i2],0.0);
    ROOT::Math::PtEtaPhiMVector p2(pt[i3],eta[i3],phi[i3],0.0);
    ROOT::Math::PtEtaPhiMVector p3(pt[i4],eta[i4],phi[i4],0.0);
    float raw_m = (p0+p1+p2+p3).M();
    ROOT::Math::PtEtaPhiMVector pA = p0+p1;
    ROOT::Math::PtEtaPhiMVector pB = p2+p3;
  

    //first pairing
    std::vector<float> pairing1_A = calc->getVertexInfo(pt[i1],eta[i1],phi[i1],EE[i1],EB[i1],pt[i2],eta[i2],phi[i2],EE[i2],EB[i2],mass); 
    std::vector<float> pairing1_B = calc->getVertexInfo(pt[i3],eta[i3],phi[i3],EE[i3],EB[i3],pt[i4],eta[i4],phi[i4],EE[i4],EB[i4],mass); 
    
    //second pairing
    std::vector<float> pairing2_A = calc->getVertexInfo(pt[i1],eta[i1],phi[i1],EE[i1],EB[i1],pt[i3],eta[i3],phi[i3],EE[i3],EB[i3],mass); 
    std::vector<float> pairing2_B = calc->getVertexInfo(pt[i2],eta[i2],phi[i2],EE[i2],EB[i2],pt[i4],eta[i4],phi[i4],EE[i4],EB[i4],mass); 
    
    //third pairing
    std::vector<float> pairing3_A = calc->getVertexInfo(pt[i1],eta[i1],phi[i1],EE[i1],EB[i1],pt[i4],eta[i4],phi[i4],EE[i4],EB[i4],mass); 
    std::vector<float> pairing3_B = calc->getVertexInfo(pt[i2],eta[i2],phi[i2],EE[i2],EB[i2],pt[i3],eta[i3],phi[i3],EE[i3],EB[i3],mass); 
    
    
    
    if (compare_quad_pairing(pairing2_A,pairing2_B,pairing3_A,pairing3_B)) {
      best23_A.insert(best23_A.end(),pairing2_A.begin(),pairing2_A.end());
      best23_B.insert(best23_B.end(),pairing2_B.begin(),pairing2_B.end());
      g1 = i1;
      g2 = i3;
      g3 = i2;
      g4 = i4;
      pA = p0+p2;
      pB = p1+p3;
      
    }
    else {
      best23_A.insert(best23_A.end(),pairing3_A.begin(),pairing3_A.end());
      best23_B.insert(best23_B.end(),pairing3_B.begin(),pairing3_B.end());
      g1 = i1;
      g2 = i4;
      g3 = i2;
      g4 = i3;
      pA=p0+p3;
      pB=p1+p2;
    }
    if (compare_quad_pairing(pairing1_A,pairing1_B,best23_A,best23_B)) {
      best_A.insert(best_A.end(),pairing1_A.begin(),pairing1_A.end());
      best_B.insert(best_B.end(),pairing1_B.begin(),pairing1_B.end());
      g1 = i1;
      g2 = i2;
      g3 = i3;
      g4 = i4;
      pA=p0+p1;
      pB=p2+p3;
    }
    else {
      best_A.insert(best_A.end(),best23_A.begin(),best23_A.end());
      best_B.insert(best_B.end(),best23_B.begin(),best23_B.end());
    }

    result.emplace_back(best_A[0]);
    result.emplace_back(best_A[1]);
    result.emplace_back(best_A[2]);
    result.emplace_back(best_A[3]);
    result.emplace_back(best_A[4]);
    result.emplace_back(best_A[5]);
    result.emplace_back(best_A[6]);
    result.emplace_back(best_A[7]);
    result.emplace_back(best_B[0]);
    result.emplace_back(best_B[1]);
    result.emplace_back(best_B[2]);
    result.emplace_back(best_B[3]);
    result.emplace_back(best_B[4]);
    result.emplace_back(best_B[5]);
    result.emplace_back(best_B[6]);
    result.emplace_back(best_B[7]);
    //raw masses
    result.emplace_back(pA.M());
    result.emplace_back(pB.M());
    result.emplace_back(raw_m);
    //corrected mass
    ROOT::Math::PtEtaPhiMVector pc0(best_A[0],best_A[1],best_A[2],0.0);
    ROOT::Math::PtEtaPhiMVector pc1(best_A[3],best_A[4],best_A[5],0.0);
    ROOT::Math::PtEtaPhiMVector pc2(best_B[0],best_B[1],best_B[2],0.0);
    ROOT::Math::PtEtaPhiMVector pc3(best_B[3],best_B[4],best_B[5],0.0);
    result.emplace_back((pc0+pc1+pc2+pc3).M());
	
    //indices of used photons
    result.emplace_back(gID[g1] && gIso[g1]<0.1);
    result.emplace_back(gID[g2] && gIso[g2]<0.1);
    result.emplace_back(gID[g3] && gIso[g3]<0.1);
    result.emplace_back(gID[g4] && gIso[g4]<0.1);
    result.emplace_back(g1);
    result.emplace_back(g2);
    result.emplace_back(g3);
    result.emplace_back(g4);
    all_combos.emplace_back(result);
  }
  delete calc;  
  if (all_combos.size()>1) {
    auto sortedIndices = ROOT::VecOps::Argsort(all_combos,compare_quad);
    return all_combos[sortedIndices[0]];
  }
  else {
    return all_combos[0];
  }


}



RVecF minMatchDR(const RVecF& photonEta,const RVecF& photonPhi,const RVecF& genEta,const RVecF& genPhi) {
  RVecF output;
  for (int i=0;i<photonEta.size();++i) {
    float minDR=1000000;
    for (int j=0;j<genEta.size();++j) {
      float dr = ROOT::VecOps::DeltaR(photonEta[i],genEta[j],photonPhi[i],genPhi[j]);
      if (dr<minDR)
	minDR=dr;
    }
    output.emplace_back(minDR);
  }
  return output;
}





//add methods to deal with the systematics of the photon scale
RVecF ptSmearUp(RVecF pt, RVecF eta, RVecF phi, RVecF scale){
   RVecF out;
   for (size_t i = 0; i < pt.size(); i++){
      ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], 0);
      float eCorr = p4.energy() + std::abs(scale[i]);
      float scale = eCorr / p4.energy();
      out.emplace_back(p4.pt() * scale);
   }
   return out;
}

RVecF ptSmearDown(RVecF pt, RVecF eta, RVecF phi, RVecF scale){
   RVecF out;
   for (size_t i = 0; i < pt.size(); i++){
      ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], 0);
      float eCorr = p4.energy() - std::abs(scale[i]);
      float scale = eCorr / p4.energy();
      out.emplace_back(p4.pt() * scale);
   }
   return out;
}

RVecF kinfit_systematics(RVecF pt, RVecF eta, RVecF phi, RVecF isEB, RVecF isEE, const float idx1, const float idx2, const float mass){
  RVecF out;
  out.reserve(2);
  VertexCalculator *calc = new VertexCalculator();
  std::vector<float> kin_fit = calc->getVertexInfo(pt[idx1], eta[idx1], phi[idx1], isEE[idx1], isEB[idx1], pt[idx2], eta[idx2], phi[idx2], isEE[idx2], isEB[idx2], mass);
  delete calc;
  ROOT::Math::PtEtaPhiMVector g1(pt[idx1], eta[idx1], phi[idx1], 0);
  ROOT::Math::PtEtaPhiMVector g2(pt[idx2], eta[idx2], phi[idx2], 0);
  out.emplace_back(kin_fit[6]);
  out.emplace_back((g1+g2).mass());
  return out;
}
