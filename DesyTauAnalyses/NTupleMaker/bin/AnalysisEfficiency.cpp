#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

int binNumber(float x, int nbins, float * bins) {

  int binN = 0;

  for (int iB=0; iB<nbins; ++iB) {
    if (x>=bins[iB]&&x<bins[iB+1]) {
      binN = iB;
      break;
    }
  }

  return binN;

}

float effBin(float x, int nbins, float * bins, float * eff) {

  int bin = binNumber(x, nbins, bins);

  return eff[bin];

}

double cosRestFrame(TLorentzVector boost, TLorentzVector vect) {

  double bx = -boost.Px()/boost.E();
  double by = -boost.Py()/boost.E();
  double bz = -boost.Pz()/boost.E();

  vect.Boost(bx,by,bz);
  double prod = -vect.Px()*bx-vect.Py()*by-vect.Pz()*bz;
  double modBeta = TMath::Sqrt(bx*bx+by*by+bz*bz); 
  double modVect = TMath::Sqrt(vect.Px()*vect.Px()+vect.Py()*vect.Py()+vect.Pz()*vect.Pz());
  
  double cosinus = prod/(modBeta*modVect);

  return cosinus;

}

double QToEta(double Q) {
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;
}

double EtaToQ(double Eta) {
  double Q = 2.0*TMath::ATan(TMath::Exp(-Eta));
  if (Q<0.0) Q += TMath::Pi();
  return Q;
}

double PtoEta(double Px, double Py, double Pz) {

  double P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
  double cosQ = Pz/P;
  double Q = TMath::ACos(cosQ);
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;

}

double PtoPhi(double Px, double Py) {
  return TMath::ATan2(Py,Px);
}

double PtoPt(double Px, double Py) {
  return TMath::Sqrt(Px*Px+Py*Py);
}

double dPhiFrom2P(double Px1, double Py1,
		  double Px2, double Py2) {


  double prod = Px1*Px2 + Py1*Py2;
  double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);
  
  double cosDPhi = prod/(mod1*mod2);
  
  return TMath::ACos(cosDPhi);

}

double deltaEta(double Px1, double Py1, double Pz1,
		double Px2, double Py2, double Pz2) {

  double eta1 = PtoEta(Px1,Py1,Pz1);
  double eta2 = PtoEta(Px2,Py2,Pz2);

  double dEta = eta1 - eta2;

  return dEta;

}

double deltaR(double Eta1, double Phi1,
	      double Eta2, double Phi2) {

  double Px1 = TMath::Cos(Phi1);
  double Py1 = TMath::Sin(Phi1);

  double Px2 = TMath::Cos(Phi2);
  double Py2 = TMath::Sin(Phi2);

  double dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
  double dEta = Eta1 - Eta2;

  double dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

  return dR;

}

double PtEtaToP(double Pt, double Eta) {

  //  double Q = EtaToQ(Eta);

  //double P = Pt/TMath::Sin(Q);
  double P = Pt*TMath::CosH(Eta);

  return P;
}
double Px(double Pt, double Phi){

  double Px=Pt*TMath::Cos(Phi);
  return Px;
}
double Py(double Pt, double Phi){

  double Py=Pt*TMath::Sin(Phi);
  return Py;
}
double Pz(double Pt, double Eta){

  double Pz=Pt*TMath::SinH(Eta);
  return Pz;
}
double InvariantMass(double energy,double Px,double Py, double Pz){

  double M_2=energy*energy-Px*Px-Py*Py-Pz*Pz;
  double M=TMath::Sqrt(M_2);
  return M;


}
double EFromPandM0(double M0,double Pt,double Eta){

  double E_2=M0*M0+PtEtaToP(Pt,Eta)*PtEtaToP(Pt,Eta);
  double E =TMath::Sqrt(E_2);
  return E;

}
bool electronMvaIdTight(float eta, float mva) {

  float absEta = fabs(eta);

  bool passed = false;
  if (absEta<0.8) {
    if (mva>0.73) passed = true;
  }
  else if (absEta<1.479) {
    if (mva>0.57) passed = true;
  }
  else {
    if (mva>0.05) passed = true;
  }

  return passed;

}

const float electronMass = 0;
const float muonMass = 0.10565837;
const float pionMass = 0.1396;

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // kinematic cuts on electrons
  const float ptElectronLowCut   = cfg.get<float>("ptElectronLowCut");
  const float ptElectronHighCut  = cfg.get<float>("ptElectronHighCut");
  const float etaElectronCut     = cfg.get<float>("etaElectronCut");
  const float dxyElectronCut     = cfg.get<float>("dxyElectronCut");
  const float dzElectronCut      = cfg.get<float>("dzElectronCut");
  const float isoElectronLowCut  = cfg.get<float>("isoElectronLowCut");
  const float isoElectronHighCut = cfg.get<float>("isoElectronHighCut");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonCut     = cfg.get<float>("etaMuonCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonLowCut  = cfg.get<float>("isoMuonLowCut");
  const float isoMuonHighCut = cfg.get<float>("isoMuonHighCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dZetaCut       = cfg.get<float>("dZetaCut");
  const bool oppositeSign    = cfg.get<bool>("oppositeSign");

  const bool applyTriggerMatch = cfg.get<bool>("ApplyTriggerMatch");
  const float DRTrigMatch      = cfg.get<float>("DRTrigMatch");

  // **** end of configuration

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");
  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);

  int nEtaBins = 3;
  float etaBins[4] = {-0.1, 0.8, 1.479, 2.3};

  int nPtBins = 7;
  float ptBins[8] = {10,15,20,30,40,50,70,100};

  TString EtaBinName[3] = {"0To0p8",
			   "0p8To1p5",
			   "1p5To2p3"};

  TString PtBinName[7] = {"pt10to15",
			  "pt15to20",
			  "pt20to30",
			  "pt30to40",
			  "pt40to50",
			  "pt50to70",
			  "pt70to100"};


  TString PassLevel[13] = {"True","TrueAll","TruePassed","RecoAll","RecoPassed",
			   "TrigDenHigh","TrigNumHigh",
			   "TrigDenLow","TrigNumLow",
			   "TPAll","TPPassed",
			   "RecoAllTag","RecoPassedTag"};
  TString Type[2] = {"Prompt","NonPrompt"};
  TString LepQ[2] = {"Pos","Neg"};

  // histograms for lepton Id
  TH1F * MuonPtH[2][2][13][3];
  TH1F * ElectronPtH[2][2][13][3];
  for (int iLepQ=0; iLepQ<2; ++iLepQ) {
    for (int iType=0; iType<2; ++iType) {
      for (int iPassLevel=0; iPassLevel<13; ++iPassLevel) {
	for (int iEta=0; iEta<3; ++iEta) {
	  MuonPtH[iLepQ][iType][iPassLevel][iEta] = new TH1F("MuonPt_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0.,100);
	  ElectronPtH[iLepQ][iType][iPassLevel][iEta] = new TH1F("ElectronPt_"+LepQ[iLepQ]+Type[iType]+PassLevel[iPassLevel]+EtaBinName[iEta],"",100,0.,100);
	}
      }
    }
  }

  TString passName[2] = {"Pass","Fail"};

  TH1F * massMuMuId[2][3][7][2];
  TH1F * massEEId[2][3][7][2];
  
  for (int iQ=0; iQ<2; ++iQ) {
    for (int iEta=0; iEta<3; ++iEta) {
      for (int iPt=0; iPt<7; ++iPt) {
	for (int iPass=0; iPass<2; ++iPass) {
	  massMuMuId[iQ][iEta][iPt][iPass] = new TH1F("massMuMuId_"+LepQ[iQ]+EtaBinName[iEta]+PtBinName[iPt]+passName[iPass],"",30,60,120);
	  massEEId[iQ][iEta][iPt][iPass] = new TH1F("massEEId_"+LepQ[iQ]+EtaBinName[iEta]+PtBinName[iPt]+passName[iPass],"",30,60,120);
	}
      }
    }
  }



  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  
    if (_tree==NULL) continue;
    
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    = " << NE << std::endl;
    
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
    
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      float weight = 1;

      //      std::cout << "Entry : " << iEntry << std::endl;
      //      std::cout << "Number of gen particles = " << analysisTree.genparticles_count << std::endl;
      //      std::cout << "Number of taus  = " << analysisTree.tau_count << std::endl;
      //      std::cout << "Number of jets  = " << analysisTree.pfjet_count << std::endl;
      //      std::cout << "Number of muons = " << analysisTree.muon_count << std::endl;
      
      // **** Analysis of generator info
      vector<unsigned int> indexZ; indexZ.clear();
      vector<unsigned int> indexMu; indexMu.clear();
      vector<unsigned int> indexE; indexE.clear();
      vector<unsigned int> indexTau; indexTau.clear();

      for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

       	float pxGen = analysisTree.genparticles_px[igen];
       	float pyGen = analysisTree.genparticles_py[igen];
       	float pzGen = analysisTree.genparticles_pz[igen];
       	float etaGen = PtoEta(pxGen,pyGen,pzGen);
       	float ptGen  = PtoPt(pxGen,pyGen);
	float phiGen = TMath::ATan2(pyGen,pxGen);

       	if (fabs(analysisTree.genparticles_pdgid[igen])==23 && 
	    (analysisTree.genparticles_status[igen]==62||analysisTree.genparticles_status[igen]==52)) { 
       	  indexZ.push_back(igen);
	  // std::cout << "Z boson :  pt = " << ptGen << "  eta = " << etaGen << "   phi = " << phiGen 
	  //   	    << "   status = " << analysisTree.genparticles_status[igen] <<  std::endl; 
	}
       	if (fabs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_status[igen==1]) {
	  // std::cout << "muon Id = " <<  analysisTree.genparticles_pdgid[igen]
	  //  	    << "   pt = "<< ptGen << "  eta = " << etaGen << "   phi = " << phiGen 
	  //  	    << "   decay = " << analysisTree.genparticles_info[igen] << std::endl; 
	  if (analysisTree.genparticles_info[igen]==1||analysisTree.genparticles_info[igen]==5)
	    indexMu.push_back(igen);
       	}
       	if (fabs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_status[igen==1]) {
	  // std::cout << "electron Id = " <<  analysisTree.genparticles_pdgid[igen]
	  //   	    << "   pt = "<< ptGen << "  eta = " << etaGen << "   phi = " << phiGen 
	  //   	    << "   decay = " << analysisTree.genparticles_info[igen] << std::endl; 
	  if (analysisTree.genparticles_info[igen]==1||analysisTree.genparticles_info[igen]==5)
	    indexE.push_back(igen);
       	}
      }
      for (unsigned int itau = 0; itau < analysisTree.gentau_count; ++ itau) {
       	float pxGen = analysisTree.gentau_px[itau];
	float pyGen = analysisTree.gentau_py[itau];
	float pzGen = analysisTree.gentau_pz[itau];
	float etaGen = PtoEta(pxGen,pyGen,pzGen);
	float phiGen = TMath::ATan2(pyGen,pxGen);
	float ptGen  = PtoPt(pxGen,pyGen);
	indexTau.push_back(itau);
       	// std::cout << "tau    pt = "<< ptGen << "  eta = " << etaGen << "   phi = " << phiGen
       	// 	  << "   decay = " <<  analysisTree.gentau_decayMode[itau] << std::endl; 
      }

      bool isZToEE = false;
      bool isZToEEopen = false;
      bool isZToMM = false;
      bool isZToMMopen = false;

      if (indexMu.size()==2) {
	unsigned int index1 = indexMu.at(0);
	unsigned int index2 = indexMu.at(1);
	int q1 = analysisTree.genparticles_pdgid[index1];
	int q2 = analysisTree.genparticles_pdgid[index2];
	if ( (q1*q2<0) && 
	     (analysisTree.genparticles_info[index1]==1) &&
	     (analysisTree.genparticles_info[index2]==1) ) { 
	  isZToMM = true;
	  
	  float pxGen1 = analysisTree.genparticles_px[index1];
	  float pyGen1 = analysisTree.genparticles_py[index1];
	  float pzGen1 = analysisTree.genparticles_pz[index1];
	  float etaGen1 = PtoEta(pxGen1,pyGen1,pzGen1);
	  float phiGen1 = TMath::ATan2(pyGen1,pxGen1);
	  
	  float pxGen2 = analysisTree.genparticles_px[index2];
	  float pyGen2 = analysisTree.genparticles_py[index2];
	  float pzGen2 = analysisTree.genparticles_pz[index2];
	  float etaGen2 = PtoEta(pxGen2,pyGen2,pzGen2);
	  float phiGen2 = TMath::ATan2(pyGen2,pxGen2);

	  float DR = deltaR(etaGen1,phiGen1,
			    etaGen2,phiGen2);
	  if (DR>dRleptonsCut)
	    isZToMMopen = true;

	}
      }

      if (indexE.size()==2) {
	unsigned int index1 = indexE.at(0);
	unsigned int index2 = indexE.at(1);
	int q1 = analysisTree.genparticles_pdgid[index1];
	int q2 = analysisTree.genparticles_pdgid[index2];
	if ( (q1*q2<0) && 
	     (analysisTree.genparticles_info[index1]==1) &&
	     (analysisTree.genparticles_info[index2]==1) ) { 
	  isZToEE = true;
	  
	  float pxGen1 = analysisTree.genparticles_px[index1];
	  float pyGen1 = analysisTree.genparticles_py[index1];
	  float pzGen1 = analysisTree.genparticles_pz[index1];
	  float etaGen1 = PtoEta(pxGen1,pyGen1,pzGen1);
	  float phiGen1 = TMath::ATan2(pyGen1,pxGen1);
	  
	  float pxGen2 = analysisTree.genparticles_px[index2];
	  float pyGen2 = analysisTree.genparticles_py[index2];
	  float pzGen2 = analysisTree.genparticles_pz[index2];
	  float etaGen2 = PtoEta(pxGen2,pyGen2,pzGen2);
	  float phiGen2 = TMath::ATan2(pyGen2,pxGen2);

	  float DR = deltaR(etaGen1,phiGen1,
			    etaGen2,phiGen2);

	  if (DR>dRleptonsCut)
	    isZToEEopen = true;
	}
      }
      bool isZToTauTau = (!isZToEE) && (!isZToMM);

      // ****************
      // trigger matching
      // ****************

      bool isMu23fired  = false;
      bool isMu8fired   = false;
      bool isEle23fired = false;
      bool isEle12fired = false;

      for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {

	bool trigFilter = false;
	TString filterName;
	float etaTrig = analysisTree.trigobject_eta[iT];
	float phiTrig = analysisTree.trigobject_phi[iT];
	float ptTrig  = analysisTree.trigobject_pt[iT];

	if (analysisTree.trigobject_filters[iT][2]) { // Mu23 leg
	  trigFilter = true;
	  filterName = "Mu23";
	  for (unsigned int J=0; J<indexMu.size(); ++J) {
	    int indexJ = indexMu.at(J);
	    if (analysisTree.genparticles_info[indexJ]!=5) continue;
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isMu23fired = true;
	      break;
	    }
	  }
	}

	if (analysisTree.trigobject_filters[iT][6]) { // Ele12 leg
	  trigFilter = true;
	  filterName = "Ele12";
	  for (unsigned int J=0; J<indexE.size(); ++J) {
	    int indexJ = indexE.at(J);
	    if (analysisTree.genparticles_info[indexJ]!=5) continue;
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isEle12fired = true;
	      break;
	    }
	  }
	}

	if (analysisTree.trigobject_filters[iT][3]) { // Mu8 leg
	  trigFilter = true;
	  filterName = "Mu8";
	  for (unsigned int J=0; J<indexMu.size(); ++J) {
	    int indexJ = indexMu.at(J);
	    if (analysisTree.genparticles_info[indexJ]!=5) continue;
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isMu8fired = true;
	      break;
	    }
	  }
	}

	if (analysisTree.trigobject_filters[iT][7]) { // Ele23 leg
	  trigFilter = true;
	  filterName = "Ele23";
	  for (unsigned int J=0; J<indexE.size(); ++J) {
	    int indexJ = indexE.at(J);
	    if (analysisTree.genparticles_info[indexJ]!=5) continue;
	    float pxGen = analysisTree.genparticles_px[indexJ];
	    float pyGen = analysisTree.genparticles_py[indexJ];
	    float pzGen = analysisTree.genparticles_pz[indexJ];
	    float etaGen = PtoEta(pxGen,pyGen,pzGen);
	    float phiGen = TMath::ATan2(pyGen,pxGen);
	    float deltaRtrig = deltaR(etaGen,phiGen,
				      etaTrig,phiTrig);
	    if (deltaRtrig<DRTrigMatch) {
	      isEle23fired = true;
	      break;
	    }
	  }
	}

	// if (trigFilter) {
	//   std::cout << filterName 
	// 	    << "  pt = " << ptTrig
	// 	    << "  eta = " << etaTrig
	// 	    << "  phi = " << phiTrig << std::endl;
	// }

      }
      
      // std::cout << "Z Decay : " << isZToEE << isZToMM << isZToTauTau << std::endl; 
      // std::cout << std::endl;


      // *********************************************
      // selecting tag and probe electron (highest pt)
      // *********************************************

      int indexTagPosE = -1;
      int indexTagNegE = -1;

      float ptMinTagPosE = 0;
      float ptMinTagNegE = 0;

      int indexProbePosE = -1;
      int indexProbeNegE = -1;

      float ptMinProbePosE = 0;
      float ptMinProbeNegE = 0;
     
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]>ptElectronLowCut&&
	    analysisTree.electron_pt[ie]<99.99&&
	    fabs(analysisTree.electron_eta[ie])<etaElectronCut) {
	  if (analysisTree.electron_pt[ie]>ptMinProbePosE&&analysisTree.electron_charge[ie]>0.5) {
	    indexProbePosE = ie;
	    ptMinProbePosE = analysisTree.electron_pt[ie];
	  }
	  if (analysisTree.electron_pt[ie]>ptMinProbeNegE&&analysisTree.electron_charge[ie]<-0.5) {
	    indexProbeNegE = ie;
	    ptMinProbeNegE = analysisTree.electron_pt[ie];
	  }
	}
	if (analysisTree.electron_pt[ie]<ptElectronHighCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	float neutralIso = 
	  analysisTree.electron_neutralHadIso[ie] + 
	  analysisTree.electron_photonIso[ie] - 
	  0.5*analysisTree.electron_puIso[ie];
	neutralIso = TMath::Max(float(0),neutralIso); 
	float absIso = analysisTree.electron_chargedHadIso[ie] + neutralIso;
	float relIso = absIso/analysisTree.electron_pt[ie];
	if (relIso>isoElectronHighCut) continue;
	if (relIso<isoElectronLowCut) continue;
	bool electronMvaId = electronMvaIdTight(analysisTree.electron_superclusterEta[ie],
						analysisTree.electron_mva_id_nontrigPhys14[ie]);
	if (!electronMvaId&&applyElectronId) continue;
	if (!analysisTree.electron_pass_conversion[ie]&&applyElectronId) continue;
	bool trigMatch = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][4]) { // Single Ele27 Leg
	    float dRtrig = deltaR(analysisTree.electron_eta[ie],analysisTree.electron_phi[ie],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<DRTrigMatch) {
	      trigMatch = true;
	      break;
	    }
	  }
	}
	if (!trigMatch) continue;
	if (analysisTree.electron_pt[ie]>ptMinTagPosE&&analysisTree.electron_charge[ie]>0.5) {
	  indexTagPosE = ie;
	  ptMinTagPosE = analysisTree.electron_pt[ie];
	}
	if (analysisTree.electron_pt[ie]>ptMinTagNegE&&analysisTree.electron_charge[ie]<-0.5) {
	  indexTagNegE = ie;
	  ptMinTagNegE = analysisTree.electron_pt[ie];
	}
      }

      // *****************************
      // Loop over generator electrons 
      // *****************************
      for (unsigned int iE=0; iE<indexE.size(); ++iE) {
	unsigned int index = indexE.at(iE);
	int electronQ = 0; 
	if (analysisTree.genparticles_pdgid[index]>0) 
	  electronQ = 1;
	TLorentzVector electronGen4P; electronGen4P.SetXYZM(analysisTree.genparticles_px[index],
							    analysisTree.genparticles_py[index],
							    analysisTree.genparticles_pz[index],
							    electronMass);
	float electronGenEta = electronGen4P.Eta();
	float electronGenPt  = electronGen4P.Pt();
	int indexReco = -1;

	int etaBinGen = binNumber(fabs(electronGenEta),nEtaBins,etaBins);

	float dPMin = 0.05;
	int matchedElectron = -1;
	if (isZToEE) {
	  ElectronPtH[electronQ][0][0][etaBinGen]->Fill(electronGenPt,weight);
	}
	else {
	  ElectronPtH[electronQ][1][0][etaBinGen]->Fill(electronGenPt,weight);
	}

	for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	  TLorentzVector electronReco4P; electronReco4P.SetXYZM(analysisTree.electron_px[ie],
								analysisTree.electron_py[ie],
								analysisTree.electron_pz[ie],
								electronMass);
	  float dP = (electronGen4P-electronReco4P).P();
	  float dPoverP = dP/electronGen4P.P();
	  if (dPoverP<dPMin && 
	      fabs(analysisTree.electron_eta[ie])<etaElectronCut &&
	      analysisTree.electron_pt[ie]>ptElectronLowCut &&
	      analysisTree.electron_pt[ie]<99.99) {
	    dPMin = dPoverP;
	    matchedElectron = ie;
	  }
	}
	if (matchedElectron>=0) {
	  float electronRecoPt = analysisTree.electron_pt[matchedElectron];
	  float electronRecoEta = analysisTree.electron_eta[matchedElectron];
	  int etaBinReco = binNumber(fabs(electronRecoEta),nEtaBins,etaBins);
	  if (isZToEE) {
	    ElectronPtH[electronQ][0][1][etaBinGen]->Fill(electronGenPt,weight);
	    if (isZToEEopen) {
	      ElectronPtH[electronQ][0][3][etaBinReco]->Fill(electronRecoPt,weight);
	      if (electronQ==0&&indexTagNegE>=0) 
		ElectronPtH[electronQ][0][11][etaBinReco]->Fill(electronRecoPt,weight);
	      if (electronQ==1&&indexTagPosE>=0) 
		ElectronPtH[electronQ][0][11][etaBinReco]->Fill(electronRecoPt,weight);
	    }
	  }
	  else {
	    ElectronPtH[electronQ][1][1][etaBinGen]->Fill(electronGenPt,weight);
	    ElectronPtH[electronQ][1][3][etaBinReco]->Fill(electronRecoPt,weight);
	  }
	  if (fabs(analysisTree.electron_dxy[matchedElectron])>dxyElectronCut) continue;
	  if (fabs(analysisTree.electron_dz[matchedElectron])>dzElectronCut) continue;
	  float neutralIso = 
	    analysisTree.electron_neutralHadIso[matchedElectron] + 
	    analysisTree.electron_photonIso[matchedElectron] - 
	    0.5*analysisTree.electron_puIso[matchedElectron];
	  neutralIso = TMath::Max(float(0),neutralIso); 
	  float absIso = analysisTree.electron_chargedHadIso[matchedElectron] + neutralIso;
	  float relIso = absIso/analysisTree.electron_pt[matchedElectron];
	  if (relIso>isoElectronHighCut) continue;
	  if (relIso<isoElectronLowCut) continue;
	  bool electronMvaId = electronMvaIdTight(analysisTree.electron_superclusterEta[matchedElectron],
						  analysisTree.electron_mva_id_nontrigPhys14[matchedElectron]);
	  if (!electronMvaId&&applyElectronId) continue;
	  if (!analysisTree.electron_pass_conversion[matchedElectron]&&applyElectronId) continue;
	  if (isZToEE) {
	    ElectronPtH[electronQ][0][2][etaBinGen]->Fill(electronGenPt,weight);
	    if (isZToEEopen) { 
	      ElectronPtH[electronQ][0][4][etaBinReco]->Fill(electronRecoPt,weight);
	      if (electronQ==0&&indexTagNegE>=0) 
		ElectronPtH[electronQ][0][12][etaBinReco]->Fill(electronRecoPt,weight);
	      if (electronQ==1&&indexTagPosE>=0) 
		ElectronPtH[electronQ][0][12][etaBinReco]->Fill(electronRecoPt,weight);
	    }
	  }
	  else {
	    ElectronPtH[electronQ][1][2][etaBinGen]->Fill(electronGenPt,weight);
	    ElectronPtH[electronQ][1][4][etaBinReco]->Fill(electronRecoPt,weight);
	    if (isMu8fired) {
	      ElectronPtH[electronQ][1][5][etaBinReco]->Fill(electronRecoPt,weight);
	      if (isEle23fired)
		ElectronPtH[electronQ][1][6][etaBinReco]->Fill(electronRecoPt,weight);
	    }
	    if (isMu23fired) {
	      ElectronPtH[electronQ][1][7][etaBinReco]->Fill(electronRecoPt,weight);
	      if (isEle12fired)
		ElectronPtH[electronQ][1][8][etaBinReco]->Fill(electronRecoPt,weight);
	    }
	  }
	}
      }


      // ******************************************
      // selecting tag and probe muons (highest pt)
      // ******************************************

      int indexTagPosMu = -1;
      int indexTagNegMu = -1;

      float ptMinTagPosMu = 0;
      float ptMinTagNegMu = 0;

      int indexProbePosMu = -1;
      int indexProbeNegMu = -1;

      float ptMinProbePosMu = 0;
      float ptMinProbeNegMu = 0;

      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]>ptMuonLowCut&&
	    analysisTree.muon_pt[im]<99.99&&
	    fabs(analysisTree.muon_eta[im])<etaMuonCut) {
	  if (analysisTree.muon_pt[im]>ptMinProbePosMu&&analysisTree.muon_charge[im]>0.5 ) {
	    indexProbePosMu = im;
	    ptMinProbePosMu = analysisTree.muon_pt[im];
	  }
	  if (analysisTree.muon_pt[im]>ptMinProbeNegMu&&analysisTree.muon_charge[im]<-0.5) {
	    indexProbeNegMu = im;
	    ptMinProbeNegMu = analysisTree.muon_pt[im];
	  }
	}
	if (analysisTree.muon_pt[im]<ptMuonHighCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	float neutralIso = 
	  analysisTree.muon_neutralHadIso[im] + 
	  analysisTree.muon_photonIso[im] - 
	  0.5*analysisTree.muon_puIso[im];
	neutralIso = TMath::Max(float(0),neutralIso); 
	float absIso = analysisTree.muon_chargedHadIso[im] + neutralIso;
	float relIso = absIso/analysisTree.muon_pt[im];
	if (relIso>isoMuonHighCut) continue;
	if (relIso<isoMuonLowCut) continue;
	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
	bool trigMatch = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][0]||analysisTree.trigobject_filters[iT][1]) { // Single IsoMu24 or IsoTkMu24 Leg
	    float dRtrig = deltaR(analysisTree.muon_eta[im],analysisTree.muon_phi[im],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<DRTrigMatch) {
	      trigMatch = true;
	      break;
	    }
	  }
	}
	if (!trigMatch) continue;
	if (analysisTree.muon_pt[im]>ptMinTagPosMu&&analysisTree.muon_charge[im]>0.5) {
	  indexTagPosMu = im;
	  ptMinTagPosMu = analysisTree.muon_pt[im];
	}
	if (analysisTree.muon_pt[im]>ptMinTagNegMu&&analysisTree.muon_charge[im]<-0.5) {
	  indexTagNegMu = im;
	  ptMinTagNegMu = analysisTree.muon_pt[im];
	}
      }

      // *************************
      // Loop over generator muons
      // *************************
      for (unsigned int iM=0; iM<indexMu.size(); ++iM) {
	unsigned int index = indexMu.at(iM);
	int muonQ = 0; 
	if (analysisTree.genparticles_pdgid[index]>0) 
	  muonQ = 1;
	TLorentzVector muonGen4P; muonGen4P.SetXYZM(analysisTree.genparticles_px[index],
						    analysisTree.genparticles_py[index],
						    analysisTree.genparticles_pz[index],
						    muonMass);
	float muonGenEta = muonGen4P.Eta();
	float muonGenPt  = muonGen4P.Pt();
	int indexReco = -1;

	int etaBinGen = binNumber(fabs(muonGenEta),nEtaBins,etaBins);

	float dPMin = 0.05;
	int matchedMuon = -1;
	if (isZToMM) {
	  MuonPtH[muonQ][0][0][etaBinGen]->Fill(muonGenPt,weight);
	}
	else {
	  MuonPtH[muonQ][1][0][etaBinGen]->Fill(muonGenPt,weight);
	}

	for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	  TLorentzVector muonReco4P; muonReco4P.SetXYZM(analysisTree.muon_px[im],
							analysisTree.muon_py[im],
							analysisTree.muon_pz[im],
							muonMass);
	  float dP = (muonGen4P-muonReco4P).P();
	  float dPoverP = dP/muonGen4P.P();
	  // std::cout << "dPoverP = " << dPoverP 
	  // 	    << "  id = " << analysisTree.genparticles_pdgid[index]
	  // 	    << "   q = " << analysisTree.muon_charge[im] << std::endl;
	  if (dPoverP<dPMin && 
	      fabs(analysisTree.muon_eta[im])<etaMuonCut &&
	      analysisTree.muon_pt[im]>ptMuonLowCut &&
	      analysisTree.muon_pt[im]<99.99) {
	    dPMin = dPoverP;
	    matchedMuon = im;
	  }
	}
	if (matchedMuon>=0) {
	  float muonRecoPt = analysisTree.muon_pt[matchedMuon];
	  float muonRecoEta = analysisTree.muon_eta[matchedMuon];
	  int etaBinReco = binNumber(fabs(muonRecoEta),nEtaBins,etaBins);
	  if (isZToMM) {
	    MuonPtH[muonQ][0][1][etaBinGen]->Fill(muonGenPt,weight);
	    if (isZToMMopen) { 
	      MuonPtH[muonQ][0][3][etaBinReco]->Fill(muonRecoPt,weight);
	      if (muonQ==0&&indexTagNegMu>=0) 
		MuonPtH[muonQ][0][11][etaBinReco]->Fill(muonRecoPt,weight);
	      if (muonQ==1&&indexTagPosMu>=0) 
		MuonPtH[muonQ][0][11][etaBinReco]->Fill(muonRecoPt,weight);
	    }
	  }
	  else {
	    MuonPtH[muonQ][1][1][etaBinGen]->Fill(muonGenPt,weight);
	    MuonPtH[muonQ][1][3][etaBinReco]->Fill(muonRecoPt,weight);
	  }
	  if (fabs(analysisTree.muon_dxy[matchedMuon])>dxyMuonCut) continue;
	  if (fabs(analysisTree.muon_dz[matchedMuon])>dzMuonCut) continue;
	  float neutralIso = 
	    analysisTree.muon_neutralHadIso[matchedMuon] + 
	    analysisTree.muon_photonIso[matchedMuon] - 
	    0.5*analysisTree.muon_puIso[matchedMuon];
	  neutralIso = TMath::Max(float(0),neutralIso); 
	  float absIso = analysisTree.muon_chargedHadIso[matchedMuon] + neutralIso;
	  float relIso = absIso/analysisTree.muon_pt[matchedMuon];
	  if (relIso>isoMuonHighCut) continue;
	  if (relIso<isoMuonLowCut) continue;
	  if (applyMuonId && !analysisTree.muon_isMedium[matchedMuon]) continue;
	  if (isZToMM) {
	    MuonPtH[muonQ][0][2][etaBinGen]->Fill(muonGenPt,weight);
	    if (isZToMMopen) { 
	      MuonPtH[muonQ][0][4][etaBinReco]->Fill(muonRecoPt,weight);
	      if (muonQ==0&&indexTagNegMu>=0) 
		MuonPtH[muonQ][0][12][etaBinReco]->Fill(muonRecoPt,weight);
	      if (muonQ==1&&indexTagPosMu>=0) 
		MuonPtH[muonQ][0][12][etaBinReco]->Fill(muonRecoPt,weight);

	    }
	  }
	  else {
	    MuonPtH[muonQ][1][2][etaBinGen]->Fill(muonGenPt,weight);
	    MuonPtH[muonQ][1][4][etaBinReco]->Fill(muonRecoPt,weight);
	    if (isEle12fired) {
	      MuonPtH[muonQ][1][5][etaBinReco]->Fill(muonRecoPt,weight);
	      if (isMu23fired)
		MuonPtH[muonQ][1][6][etaBinReco]->Fill(muonRecoPt,weight);
	    }
	    if (isEle23fired) {
	      MuonPtH[muonQ][1][7][etaBinReco]->Fill(muonRecoPt,weight);
	      if (isMu8fired)
		MuonPtH[muonQ][1][8][etaBinReco]->Fill(muonRecoPt,weight);
	    }
	  }
	}
      }


      // ******************************
      // electron Id with tag-and-probe
      // ******************************
      if (analysisTree.hltriggerresults_second[0]==1) { // Single Ele27 trigger
	for (int iQ=0; iQ<2; ++iQ) {

	  int indexTag = indexTagPosE;
	  int indexProbe = indexProbeNegE;
	  if (iQ==1) {
	    indexTag = indexTagNegE;
	    indexProbe = indexProbePosE;
	  }
	  
	  if (indexTag>=0&&indexProbe>=0) { 
	    TLorentzVector tagP4; tagP4.SetXYZM(analysisTree.electron_px[indexTag],
						analysisTree.electron_py[indexTag],
						analysisTree.electron_pz[indexTag],
						electronMass);
	    TLorentzVector probeP4; probeP4.SetXYZM(analysisTree.electron_px[indexProbe],
						    analysisTree.electron_py[indexProbe],
						    analysisTree.electron_pz[indexProbe],
						    electronMass);
	    
	    // artificial matching
	    if (!isZToEE) continue;
	    float dPMin = 1;
	    for (unsigned int genE=0; genE<indexE.size(); ++genE) {
	      int index = indexE.at(genE);
	      TLorentzVector genP4; genP4.SetXYZM(analysisTree.genparticles_px[index],
						  analysisTree.genparticles_py[index],
						  analysisTree.genparticles_pz[index],
						  muonMass);

	      float dP = (genP4-probeP4).P();
	      float dPoverP = dP/genP4.P();
	      if (dPoverP<dPMin) {
		dPMin = dPoverP;
	      }
	    }
	    if (dPMin>0.05) continue;

	    float deltaREE = deltaR(analysisTree.electron_eta[indexTag],
				    analysisTree.electron_phi[indexTag],
				    analysisTree.electron_eta[indexProbe],
				    analysisTree.electron_phi[indexProbe]);
	    if (deltaREE<dRleptonsCut) continue;

	    float mass = (tagP4+probeP4).M();
	    bool probePassed = true;
	    if (fabs(analysisTree.electron_eta[indexProbe])>etaElectronCut) probePassed = false;
	    if (fabs(analysisTree.electron_dxy[indexProbe])>dxyElectronCut) probePassed = false;
	    if (fabs(analysisTree.electron_dz[indexProbe])>dzElectronCut) probePassed = false;
	    float neutralIso = 
	      analysisTree.electron_neutralHadIso[indexProbe] + 
	      analysisTree.electron_photonIso[indexProbe] - 
	      0.5*analysisTree.electron_puIso[indexProbe];
	    neutralIso = TMath::Max(float(0),neutralIso); 
	    float absIso = analysisTree.electron_chargedHadIso[indexProbe] + neutralIso;
	    float relIso = absIso/analysisTree.electron_pt[indexProbe];
	    if (relIso>isoElectronHighCut) probePassed = false;
	    if (relIso<isoElectronLowCut) probePassed = false;
	    bool electronMvaId = electronMvaIdTight(analysisTree.electron_superclusterEta[indexProbe],
						    analysisTree.electron_mva_id_nontrigPhys14[indexProbe]);
	    if (!electronMvaId&&applyElectronId) probePassed = false;
	    if (!analysisTree.electron_pass_conversion[indexProbe]&&applyElectronId) probePassed = false;
	    int etaBin = binNumber(fabs(analysisTree.electron_eta[indexProbe]),nEtaBins,etaBins);
	    int ptBin = binNumber(analysisTree.electron_pt[indexProbe],nPtBins,ptBins);
	    ElectronPtH[iQ][0][9][etaBin]->Fill(analysisTree.electron_pt[indexProbe],weight);
	    if (probePassed) {
	      massEEId[iQ][etaBin][ptBin][0]->Fill(mass,weight);
	      ElectronPtH[iQ][0][10][etaBin]->Fill(analysisTree.electron_pt[indexProbe],weight);
	    }
	    else
	      massEEId[iQ][etaBin][ptBin][1]->Fill(mass,weight);
	  }
	}
      }

      // **************************
      // muon Id with tag-and-probe
      // **************************
      if (analysisTree.hltriggerresults_second[2]==1||analysisTree.hltriggerresults_second[3]==1) { // Single IsoMu24 || IsoTkMu24 trigger
	for (int iQ=0; iQ<2; ++iQ) {

	  int indexTag = indexTagPosMu;
	  int indexProbe = indexProbeNegMu;
	  if (iQ==1) {
	    indexTag = indexTagNegMu;
	    indexProbe = indexProbePosMu;
	  }
	  
	  if (indexTag>=0&&indexProbe>=0) { 
	    TLorentzVector tagP4; tagP4.SetXYZM(analysisTree.muon_px[indexTag],
						analysisTree.muon_py[indexTag],
						analysisTree.muon_pz[indexTag],
						muonMass);
	    TLorentzVector probeP4; probeP4.SetXYZM(analysisTree.muon_px[indexProbe],
						    analysisTree.muon_py[indexProbe],
						    analysisTree.muon_pz[indexProbe],
						    muonMass);
	    // artificial matching
	    if (!isZToMM) continue;
	    float dPMin = 1;
	    for (unsigned int genMu=0; genMu<indexMu.size(); ++genMu) {
	      int index = indexMu.at(genMu);
	      TLorentzVector genP4; genP4.SetXYZM(analysisTree.genparticles_px[index],
						  analysisTree.genparticles_py[index],
						  analysisTree.genparticles_pz[index],
						  muonMass);

	      float dP = (genP4-probeP4).P();
	      float dPoverP = dP/genP4.P();
	      if (dPoverP<dPMin) {
		dPMin = dPoverP;
	      }
	    }
	    if (dPMin>0.05) continue;

	    float deltaRMuMu = deltaR(analysisTree.muon_eta[indexTag],
				      analysisTree.muon_phi[indexTag],
				      analysisTree.muon_eta[indexProbe],
				      analysisTree.muon_phi[indexProbe]);
	    if (deltaRMuMu<dRleptonsCut) continue;

	    float mass = (tagP4+probeP4).M();
	    bool probePassed = true;
	    if (fabs(analysisTree.muon_dxy[indexProbe])>dxyMuonCut) probePassed = false;
	    if (fabs(analysisTree.muon_dz[indexProbe])>dzMuonCut) probePassed = false;
	    float neutralIso = 
	      analysisTree.muon_neutralHadIso[indexProbe] + 
	      analysisTree.muon_photonIso[indexProbe] - 
	      0.5*analysisTree.muon_puIso[indexProbe];
	    neutralIso = TMath::Max(float(0),neutralIso); 
	    float absIso = analysisTree.muon_chargedHadIso[indexProbe] + neutralIso;
	    float relIso = absIso/analysisTree.muon_pt[indexProbe];
	    if (relIso>isoMuonHighCut) probePassed = false;
	    if (relIso<isoMuonLowCut) probePassed = false;
	    if (applyMuonId && !analysisTree.muon_isMedium[indexProbe]) probePassed = false;
	    int etaBin = binNumber(fabs(analysisTree.muon_eta[indexProbe]),nEtaBins,etaBins);
	    int ptBin = binNumber(analysisTree.muon_pt[indexProbe],nPtBins,ptBins);
	    MuonPtH[iQ][0][9][etaBin]->Fill(analysisTree.muon_pt[indexProbe],weight);
	    if (probePassed) {
	      massMuMuId[iQ][etaBin][ptBin][0]->Fill(mass,weight);
	      MuonPtH[iQ][0][10][etaBin]->Fill(analysisTree.muon_pt[indexProbe],weight);
	    }
	    else
	      massMuMuId[iQ][etaBin][ptBin][1]->Fill(mass,weight);
	  }
	}
      }


      selEvents++;
      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;

  file->cd("");
  file->Write();
  file->Close();
  delete file;
  
}



