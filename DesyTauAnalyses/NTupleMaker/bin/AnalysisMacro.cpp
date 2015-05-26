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

const float MuMass = 0.105658367;

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

  TH1F * muonPtAllH = new TH1F("muonPtAllH","",40,0,200);
  TH1F * electronPtAllH = new TH1F("electronPtAllH","",40,0,200);

  // histograms (dilepton selection)
  TH1F * electronPtH  = new TH1F("electronPtH","",40,0,200);
  TH1F * electronEtaH = new TH1F("electronEtaH","",50,-2.5,2.5); 
  TH1F * muonPtH  = new TH1F("muonPtH","",40,0,200);
  TH1F * muonEtaH = new TH1F("muonEtaH","",50,-2.5,2.5); 

  TH1F * dileptonMassH = new TH1F("dileptonMassH","",40,0,200);
  TH1F * dileptonPtH = new TH1F("dileptonPtH","",40,0,200);
  TH1F * dileptonEtaH = new TH1F("dileptonEtaH","",100,-5,5);
  TH1F * dileptondRH = new TH1F("dileptondRH","",60,0,6);
  TH1F * ETmissH = new TH1F("ETmissH","",40,0,200);
  TH1F * MtH = new TH1F("MtH","",40,0,200);
  TH1F * DZetaH = new TH1F("DZetaH","",60,-400,200);

  // histograms (dilepton selection + DZeta cut DZeta)
  TH1F * electronPtSelH  = new TH1F("electronPtSelH","",40,0,200);
  TH1F * electronEtaSelH = new TH1F("electronEtaSelH","",50,-2.5,2.5); 
  TH1F * muonPtSelH  = new TH1F("muonPtSelH","",40,0,200);
  TH1F * muonEtaSelH = new TH1F("muonEtaSelH","",50,-2.5,2.5); 

  TH1F * dileptonMassSelH = new TH1F("dileptonMassSelH","",40,0,200);
  TH1F * dileptonPtSelH = new TH1F("dileptonPtSelH","",40,0,200);
  TH1F * dileptonEtaSelH = new TH1F("dileptonEtaSelH","",100,-5,5);
  TH1F * dileptondRSelH = new TH1F("dileptondRSelH","",60,0,6);
  TH1F * ETmissSelH = new TH1F("ETmissSelH","",40,0,200);
  TH1F * MtSelH = new TH1F("MtSelH","",40,0,200);
  TH1F * DZetaSelH = new TH1F("DZetaSelH","",60,-400,200);

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
      // int indexW  = -1;
      // int indexNu = -1; 
      // int indexMu = -1;
      // int indexE  = -1;
      // int nGenMuons = 0;
      // int nGenElectrons = 0;
      // for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

      // 	float pxGen = analysisTree.genparticles_px[igen];
      // 	float pyGen = analysisTree.genparticles_py[igen];
      // 	float pzGen = analysisTree.genparticles_pz[igen];
      // 	float etaGen = PtoEta(pxGen,pyGen,pzGen);
      // 	float ptGen  = PtoPt(pxGen,pyGen);

      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==24 && analysisTree.genparticles_status[igen]==62) 
      // 	  indexW = igen;
      // 	if ((fabs(analysisTree.genparticles_pdgid[igen])==12 
      // 	     ||fabs(analysisTree.genparticles_pdgid[igen])==14
      // 	     ||fabs(analysisTree.genparticles_pdgid[igen])==16) 
      // 	    && analysisTree.genparticles_info[igen]== (1<<1) )
      // 	  indexNu = igen;

      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==13) {
      // 	  if ( analysisTree.genparticles_info[igen]== (1<<1) ) {
      // 	    indexMu = igen;
      // 	    if (fabs(etaGen)<2.3 && ptGen>10.)
      // 	      nGenMuons++;
      // 	  }
      // 	}
      // 	if (fabs(analysisTree.genparticles_pdgid[igen])==11) {
      // 	  if ( analysisTree.genparticles_info[igen]== (1<<1) ) {
      // 	    indexE = igen;
      // 	    if (fabs(etaGen)<2.3 && ptGen>10.)
      // 	      nGenElectrons++;
      // 	  }
      // 	}
      // }

      // trigger selection
       
      bool trigAccept = false;

      for (int i=0; i<kMaxhltriggerresults; ++i) {
	if ((i==5||i==6)&&analysisTree.hltriggerresults_second[i]==1) {
	  //	  std::cout << analysisTree.run_hltnames->at(i) << " : " << analysisTree.hltriggerresults_second[i] << std::endl;
	  trigAccept = true;
	}
      }

      if (!trigAccept) continue;
      
      // electron selection

      vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	electronPtAllH->Fill(analysisTree.electron_pt[ie],weight);
	if (analysisTree.electron_pt[ie]<ptElectronLowCut) continue;
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
	electrons.push_back(ie);
      }

      // muon selection

      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	muonPtAllH->Fill(analysisTree.muon_pt[im],weight);
	if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
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
	muons.push_back(im);
      }

      if (electrons.size()==0) continue;
      if (muons.size()==0) continue;

      // selecting muon and electron pair (OS or SS);
      float ptScalarSum = -1;
      float dRleptons = -1;
      int electronIndex = -1;
      int muonIndex = -1;
      for (unsigned int ie=0; ie<electrons.size(); ++ie) {
	int eIndex = electrons[ie];
	for (unsigned int im=0; im<muons.size(); ++im) {
	  int mIndex = muons[im];
	  float qProd = analysisTree.electron_charge[eIndex]*analysisTree.muon_charge[mIndex];
	  if (oppositeSign && qProd>0) continue;
	  if (!oppositeSign && qProd<0) continue;
	  float dR = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
			    analysisTree.muon_eta[mIndex],analysisTree.muon_phi[mIndex]);

	  if (dR<dRleptonsCut) continue;

	  bool kinematics = 
	    (analysisTree.electron_pt[eIndex]>ptElectronLowCut && analysisTree.muon_pt[mIndex]>ptMuonHighCut) || 
	    (analysisTree.electron_pt[eIndex]>ptElectronHighCut && analysisTree.muon_pt[mIndex]>ptMuonLowCut);

	  if (!kinematics) continue;

	  float sumPt = analysisTree.electron_pt[eIndex] + analysisTree.muon_pt[mIndex];
	  if (sumPt>ptScalarSum) {
	    ptScalarSum = sumPt;
	    dRleptons = dR;
	    electronIndex = ie;
	    muonIndex = im;
	  }
	  
	}
      }

      if (ptScalarSum<0) continue;

      // computations of kinematic variables

      TLorentzVector muonLV; muonLV.SetXYZM(analysisTree.muon_px[muonIndex],
					    analysisTree.muon_py[muonIndex],
					    analysisTree.muon_pz[muonIndex],
					    muonMass);

      TLorentzVector electronLV; electronLV.SetXYZM(analysisTree.electron_px[electronIndex],
						    analysisTree.electron_py[electronIndex],
						    analysisTree.electron_pz[electronIndex],
						    electronMass);

      TLorentzVector dileptonLV = muonLV + electronLV;
      float dileptonMass = dileptonLV.M();
      float dileptonPt = dileptonLV.Pt();
      float dileptonEta = dileptonLV.Eta();

      float ETmiss = TMath::Sqrt(analysisTree.pfmet_ex*analysisTree.pfmet_ex + analysisTree.pfmet_ey*analysisTree.pfmet_ey);

      // bisector of electron and muon transverse momenta
      float electronUnitX = electronLV.Px()/electronLV.Pt();
      float electronUnitY = electronLV.Py()/electronLV.Pt();
	
      float muonUnitX = muonLV.Px()/muonLV.Pt();
      float muonUnitY = muonLV.Py()/muonLV.Pt();

      float zetaX = electronUnitX + muonUnitX;
      float zetaY = electronUnitY + muonUnitY;
      
      float normZeta = TMath::Sqrt(zetaX*zetaX+zetaY*zetaY);

      zetaX = zetaX/normZeta;
      zetaY = zetaY/normZeta;

      float vectorX = analysisTree.pfmet_ex + muonLV.Px() + electronLV.Px();
      float vectorY = analysisTree.pfmet_ey + muonLV.Py() + electronLV.Py();
      
      float vectorVisX = muonLV.Px() + electronLV.Px();
      float vectorVisY = muonLV.Py() + electronLV.Py();

      // computation of DZeta variable
      float PZeta = vectorX*zetaX + vectorY*zetaY;
      float PVisZeta = vectorVisX*zetaX + vectorVisY*zetaY;
      float DZeta = PZeta - 1.85*PVisZeta;

      // computation of MT variable
      float dPhi = dPhiFrom2P( dileptonLV.Px(), dileptonLV.Py(),
			       analysisTree.pfmet_ex,  analysisTree.pfmet_ey );

      float MT = TMath::Sqrt(2*dileptonPt*ETmiss*(1-TMath::Cos(dPhi)));

      // filling histograms after dilepton selection

      electronPtH->Fill(electronLV.Pt(),weight);
      electronEtaH->Fill(electronLV.Eta(),weight);
      
      muonPtH->Fill(muonLV.Pt(),weight);
      muonEtaH->Fill(muonLV.Eta(),weight);
      
      dileptonMassH->Fill(dileptonMass,weight);
      dileptonPtH->Fill(dileptonPt,weight);
      dileptonEtaH->Fill(dileptonEta,weight);
      dileptondRH->Fill(dRleptons,weight);
      
      ETmissH->Fill(ETmiss,weight);
      MtH->Fill(MT,weight);
      DZetaH->Fill(DZeta,weight);

      // topological cut
      if (DZeta<dZetaCut) continue;
      
      electronPtSelH->Fill(electronLV.Pt(),weight);
      electronEtaSelH->Fill(electronLV.Eta(),weight);
      
      muonPtSelH->Fill(muonLV.Pt(),weight);
      muonEtaSelH->Fill(muonLV.Eta(),weight);
      
      dileptonMassSelH->Fill(dileptonMass,weight);
      dileptonPtSelH->Fill(dileptonPt,weight);
      dileptonEtaSelH->Fill(dileptonEta,weight);
      dileptondRSelH->Fill(dRleptons,weight);
      
      ETmissSelH->Fill(ETmiss,weight);
      MtSelH->Fill(MT,weight);
      DZetaSelH->Fill(DZeta,weight);

      //      std::cout << std::endl;
      
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



