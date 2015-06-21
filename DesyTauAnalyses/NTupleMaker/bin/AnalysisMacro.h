
#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

using namespace std;
const Float_t MuMass = 0.105658367;

const Float_t electronMass = 0;
const Float_t muonMass = 0.10565837;
const Float_t pionMass = 0.1396;
const  int CutN=10;
   Float_t XSec=-1;
   Float_t xs,fact,fact2;
   //string CutList[10];
    vector<string> CutList;

TH1D *hHT[CutN];
TH1D *hST[CutN];
TH1D *h0JetpT[CutN];
TH1D *hnJet[CutN];
TH1D *hnBJet[CutN];
TH1D *hnEl[CutN];
TH1D *hnMu[CutN];
TH1D *hLeppt[CutN];
TH1D *hElpt[CutN];
TH1D *hMupt[CutN];
TH1D *hLepeta[CutN];
TH1D *hEleta[CutN];
TH1D *hMueta[CutN];
TH1D *hMET[CutN];
TH1D *hnOver[CutN];
TH1D *hdPhiMETLep[CutN];
TH1D *hdPhiJMET[CutN];
TH1D *hToppT[CutN];
TH1D *hnTop[CutN];
TH1D *hnW[CutN];
TH1D *hWTagpT[CutN];
TH1D *hWmassTagpT[CutN];
TH1D *hWTagMass[CutN];
TH1D *hWmassTagMass[CutN];
TH1D *hnWmass[CutN];
TH1D *hMT[CutN];
TH1D *hDZeta[CutN];
TH1D *hel_miniISO[CutN];
TH1D *hmu_miniISO[CutN];
TH1D *hel_relISO[CutN];
TH1D *hmu_relISO[CutN];
TH2D *hmet_MT[CutN];
TH2D *hmet_dPhi[CutN];
TH2D *hMT_dPhi[CutN];
  
  TH1D *CutFlow= new TH1D("CutFlow","Cut Flow",CutN,0.5,CutN+0.5);

  TH1D * inputEventsH = new TH1D("inputEventsH","",1,-0.5,0.5);
  TH1D * hxsec = new TH1D("xsec","",1,0,10e+06);

  TH1D * muonPtAllH = new TH1D("muonPtAllH","",40,0,200);
  TH1D * electronPtAllH = new TH1D("electronPtAllH","",40,0,200);

  // histograms (dilepton selection)
  TH1D * electronPtH  = new TH1D("electronPtH","",40,0,200);
  TH1D * electronEtaH = new TH1D("electronEtaH","",50,-2.5,2.5); 
  TH1D * muonPtH  = new TH1D("muonPtH","",40,0,200);
  TH1D * muonEtaH = new TH1D("muonEtaH","",50,-2.5,2.5); 

  TH1D * dileptonMassH = new TH1D("dileptonMassH","",40,0,200);
  TH1D * dileptonPtH = new TH1D("dileptonPtH","",40,0,200);
  TH1D * dileptonEtaH = new TH1D("dileptonEtaH","",100,-5,5);
  TH1D * dileptondRH = new TH1D("dileptondRH","",60,0,6);
  TH1D * ETmissH = new TH1D("ETmissH","",40,0,200);
  TH1D * MtH = new TH1D("MtH_2l","",40,0,200);
  TH1D * DZetaH = new TH1D("DZetaH","",60,-400,200);

  // histograms (dilepton selection + DZeta cut DZeta)
  TH1D * electronPtSelH  = new TH1D("electronPtSelH","",40,0,200);
  TH1D * electronEtaSelH = new TH1D("electronEtaSelH","",50,-2.5,2.5); 
  TH1D * muonPtSelH  = new TH1D("muonPtSelH","",40,0,200);
  TH1D * muonEtaSelH = new TH1D("muonEtaSelH","",50,-2.5,2.5); 

  TH1D * dileptonMassSelH = new TH1D("dileptonMassSelH","",40,0,200);
  TH1D * dileptonPtSelH = new TH1D("dileptonPtSelH","",40,0,200);
  TH1D * dileptonEtaSelH = new TH1D("dileptonEtaSelH","",100,-5,5);
  TH1D * dileptondRSelH = new TH1D("dileptondRSelH","",60,0,6);
  TH1D * ETmissSelH = new TH1D("ETmissSelH","",40,0,200);
  TH1D * MtSelH = new TH1D("MtSelH_2l","",40,0,200);
  TH1D * DZetaSelH = new TH1D("DZetaSelH","",60,-400,200);

     TLorentzVector ElV, MuV, JetsV, METV;

     vector<TLorentzVector> JetsMV;
     vector<TLorentzVector>  ElMV;
     vector<TLorentzVector>  MuMV;
     vector<TLorentzVector>  LeptMV;

bool ComparePt(TLorentzVector a, TLorentzVector b) { return a.Pt() > b.Pt(); }


int binNumber(Float_t x, int nbins, Float_t * bins) {

  int binN = 0;

  for (int iB=0; iB<nbins; ++iB) {
    if (x>=bins[iB]&&x<bins[iB+1]) {
      binN = iB;
      break;
    }
  }

  return binN;

}

Float_t effBin(Float_t x, int nbins, Float_t * bins, Float_t * eff) {

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


double DeltaPhi(TLorentzVector METV, TLorentzVector LeptonV){
  TLorentzVector Ws = METV + LeptonV;

        //Delta phi between W and Lep
        //standard root defintion (+ fabs)takes care of getting values between 0 and pi
        double DelPhiWLep = fabs(Ws.DeltaPhi(LeptonV));
        //alternative definiton with the same result, if you want to cross check
      	 Double_t DelPhiWLepAlt = (Ws.Phi() - LeptonV.Phi());
        if (DelPhiWLepAlt > TMath::Pi()) DelPhiWLepAlt -= 2*TMath::Pi();
        if (DelPhiWLepAlt <= -TMath::Pi()) DelPhiWLepAlt += 2*TMath::Pi();
        DelPhiWLepAlt = fabs(DelPhiWLepAlt);
		
	return DelPhiWLep;

}



bool electronMvaIdTight(Float_t eta, Float_t mva) {

  Float_t absEta = fabs(eta);

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

bool electronVetoTight(Float_t SuperClusterEta, Float_t eta, Float_t phi, Float_t full5, Float_t hOverE, Float_t d0, Float_t dZ, Float_t ooE, Float_t pfISO, Float_t nMissing, bool convVeto) {


  bool passed = false;


if (fabs(SuperClusterEta) <= 1.479 ){
	
	if (  fabs(eta)<0.013625 && fabs(phi)< 0.230374 && fabs(full5) < 0.011586 && hOverE < 0.181130 && fabs(d0) < 0.094095 && fabs(dZ) < 0.713070 && ooE < 0.295751 && pfISO < 0.158721 && nMissing <= 2 && convVeto) 
	       passed = true;
}

else if (      1.479  < fabs(SuperClusterEta)  && fabs(SuperClusterEta) < 2.5 ){
 
  	if ( fabs(eta)<0.011932 && fabs(phi)< 0.255450 && fabs(full5) < 0.031849 && hOverE < 0.223870 && fabs(d0) < 0.342293 && fabs(dZ) < 0.953461 && ooE < 0.155501 && pfISO < 0.177032 && nMissing <= 3 && convVeto) 

	       passed = true;
}
       else  
	       passed = false;

  return passed;
}









//string CutList[CutN];// ={"No cut","Trigger","2- l", "dR < "};
void SetupHists(int CutNer){
    for(int cj = 0; cj < CutNer; cj++)
    {
        CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj].c_str());
        TString cutName=CutList[cj];
        TString nCut;
        nCut.Form("%d",cj);

        hHT[cj] = new TH1D ("HT_"+nCut,"HT "+cutName,400,0.0,4000.0);
        hHT[cj]->Sumw2();
        hST[cj] = new TH1D ("ST_"+nCut,"ST "+cutName,400,0.0,4000.0);
        hST[cj]->Sumw2();
        h0JetpT[cj] = new TH1D ("0JetpT_"+nCut,"0JetpT "+cutName,200,0.0,2000.0);
        h0JetpT[cj]->Sumw2();
        hnJet[cj] = new TH1D ("nJet_"+nCut,"nJet "+cutName,20,0,20);
        hnJet[cj]->Sumw2();
        hnBJet[cj] = new TH1D ("nBJet_"+nCut,"nBJet "+cutName,20,0,20);
        hnBJet[cj]->Sumw2();
        hnEl[cj] = new TH1D ("nEl_"+nCut,"nEl "+cutName,10,0,10);
        hnEl[cj]->Sumw2();
        hLeppt[cj] = new TH1D ("LeppT_"+nCut,"Lep pT "+cutName,100,0,1000);
        hLeppt[cj]->Sumw2();
        hElpt[cj] = new TH1D ("ElpT_"+nCut,"El pT "+cutName,100,0,1000);
        hElpt[cj]->Sumw2();
        hnMu[cj] = new TH1D ("nMu_"+nCut,"nMu "+cutName,10,0,10);
        hnMu[cj]->Sumw2();
        hMupt[cj] = new TH1D ("MupT_"+nCut,"Mu pT "+cutName,100,0,1000);
        hMupt[cj]->Sumw2();
        hnOver[cj] = new TH1D ("nOver_"+nCut,"nOver "+cutName,2,0,2);
        hEleta[cj] = new TH1D ("Eleta_"+nCut,"El eta "+cutName,100,-4,4);
        hEleta[cj]->Sumw2();
        hMueta[cj] = new TH1D ("Mueta_"+nCut,"Mu eta "+cutName,100,-4,4);
        hMueta[cj]->Sumw2();
        hLepeta[cj] = new TH1D ("Lepeta_"+nCut,"Lep eta "+cutName,100,-4,4);
        hLepeta[cj]->Sumw2();
        hMET[cj] = new TH1D("MET_"+nCut,"MET "+cutName,200.0,0.0,2000.0);
        hMET[cj]->Sumw2();
        hdPhiMETLep[cj] = new TH1D("dPhiMETLep_"+nCut,"dPhiMETLep "+cutName,64,0.0,3.2);
        hdPhiMETLep[cj]->Sumw2();
        hdPhiJMET[cj] = new TH1D("dPhiJMET_"+nCut,"dPhiJMET "+cutName,64,0.0,3.2);
        hdPhiJMET[cj]->Sumw2();
        hToppT[cj] = new TH1D ("ToppT_"+nCut,"ToppT "+cutName,200,0.0,2000.0);
        hToppT[cj]->Sumw2();
        hnTop[cj] = new TH1D ("nTop_"+nCut,"nTop "+cutName,20,0,20);
        hnTop[cj]->Sumw2();
        hnW[cj] = new TH1D ("nW_"+nCut,"nW "+cutName,20,0,20);
        hnW[cj]->Sumw2();
        hWTagpT[cj] = new TH1D ("WpT_"+nCut,"WpT "+cutName,200,0.0,2000.0);
        hWTagpT[cj]->Sumw2();
        hWmassTagpT[cj] = new TH1D ("WmasspT_"+nCut,"WmasspT "+cutName,200,0.0,2000.0);
        hWmassTagpT[cj]->Sumw2();
        hWTagMass[cj] = new TH1D ("WMass_"+nCut,"WMass "+cutName,200,0.0,2000.0);
        hWTagMass[cj]->Sumw2();
        hWmassTagMass[cj] = new TH1D ("WmassMass_"+nCut,"WmassMass "+cutName,200,0.0,2000.0);
        hWmassTagMass[cj]->Sumw2();
        hnWmass[cj] = new TH1D ("nWmass_"+nCut,"nWmass "+cutName,20,0,20);
        hnWmass[cj]->Sumw2();
        hMT[cj] = new TH1D ("MT_"+nCut,"MT "+cutName,40,0,200);
        hMT[cj]->Sumw2();
 /*       hDZeta[CutN]= new TH1D ("DZeta_"+nCut,"DZeta "+cutName,300,-400,200);
        hDZeta[cj]->Sumw2();
*/
        hel_miniISO[cj]= new TH1D ("elminiISO_"+nCut,"elminiISO "+cutName,50,0,5);;
        hel_miniISO[cj]->Sumw2();
        hmu_miniISO[cj]= new TH1D ("muminiISO_"+nCut,"muminiISO "+cutName,50,0,5);;
        hmu_miniISO[cj]->Sumw2();
        hel_relISO[cj]= new TH1D ("elrelISO_"+nCut,"elrelISO "+cutName,50,0,5);;
        hel_relISO[cj]->Sumw2();
        hmu_relISO[cj]= new TH1D ("murelISO_"+nCut,"murelISO "+cutName,50,0,5);;
        hmu_relISO[cj]->Sumw2();
        hmet_dPhi[cj] = new TH2D ("met_dPhi_"+nCut,"met_dPhi "+cutName,200.0,0.0,2000.0,64,0.0,3.2);
	hmet_dPhi[cj]->Sumw2();
        hmet_MT[cj] = new TH2D ("met_MT_"+nCut,"met_MT "+cutName,200.0,0.0,2000.0,40,0,200);
	hmet_MT[cj]->Sumw2();
        hMT_dPhi[cj]= new TH2D ("MT_dPhi"+nCut,"MT_dPhi "+cutName,40,0,200,64,0.0,3.2);
	hMT_dPhi[cj]->Sumw2();

    }
}


void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  ElV, vector<TLorentzVector>  MuV,vector<TLorentzVector>  JetsV, TLorentzVector  MetV, AC1B &tree_){
	//void FillMainHists(int CutIndex, Double_t EvWeight, vector<TLorentzVector>  *JetsV){
	hnJet[CutIndex]->Fill(JetsV.size(),EvWeight);
        hnEl[CutIndex]->Fill(ElV.size(),EvWeight);
        hnMu[CutIndex]->Fill(MuV.size(),EvWeight);
   //     if (JetsV.size() > 0) h0JetpT[CutIndex]->Fill(JetsV.at(0).Pt(),EvWeight);
  //  if(FillBJets){
  //      hnBJet[CutIndex]->Fill(Obj.nBJetGood,EvWeight);
  //  }
    if (ElV.size() > 0)
    {
     hElpt[CutIndex]->Fill(ElV.at(0).Pt(),EvWeight);
     hLeppt[CutIndex]->Fill(ElV.at(0).Pt(),EvWeight);
     
 //    for (unsigned int nel = 0; nel<ElV.size();nel++){
 //    hel_miniISO[CutIndex]->Fill(analysisTree.electron_miniISO[nel],EvWeight);
 //    }
     Float_t dPhi=dPhiFrom2P( ElV.at(0).Px(), ElV.at(0).Py(), MetV.Px(),  MetV.Py() );
     hdPhiMETLep[CutIndex]->Fill(dPhi,EvWeight);
      		
      Float_t MT = TMath::Sqrt(2*ElV.at(0).Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
      hMT[CutIndex]->Fill(MT,EvWeight);
      hmet_dPhi[CutIndex]->Fill(MetV.Pt(),dPhi,EvWeight);
      hmet_MT[CutIndex]->Fill(MetV.Pt(),MT,EvWeight);
      hMT_dPhi[CutIndex]->Fill(MT,dPhi,EvWeight);

   for (unsigned int ie=0;ie<tree_.electron_count;ie++){
	   if (tree_.electron_pt[ie] == ElV.at(0).Pt() ) 
	   {
		   hel_miniISO[CutIndex]->Fill(tree_.electron_miniISO[ie],EvWeight);

	Float_t neutralIso = 
	  tree_.electron_neutralHadIso[ie] + 
	  tree_.electron_photonIso[ie] - 
	  0.5*tree_.electron_puIso[ie];
	neutralIso = TMath::Max(Float_t(0),neutralIso); 
	Float_t absIso = tree_.electron_chargedHadIso[ie] + neutralIso;
	Float_t relIso = absIso/tree_.electron_pt[ie];
        hel_relISO[CutIndex]->Fill(relIso,EvWeight);
	}
   	}
    }


    if (MuV.size() > 0)
    {
     TLorentzVector WBos = MetV + MuV.at(0);
     hMupt[CutIndex]->Fill(MuV.at(0).Pt(),EvWeight);
     hLeppt[CutIndex]->Fill(MuV.at(0).Pt(),EvWeight);
  //   for (unsigned int nmu = 0; nmu<MuV.size();nmu++){
  //   hmu_miniISO[CutIndex]->Fill(analysisTree.muon_miniISO[nmu],EvWeight);
  //  }
     Float_t dPhi=dPhiFrom2P( MuV.at(0).Px(), MuV.at(0).Py(), MetV.Px(),  MetV.Py() );
     hdPhiMETLep[CutIndex]->Fill(dPhi,EvWeight);
      Float_t MT = TMath::Sqrt(2*MuV.at(0).Pt()*MetV.Pt()*(1-TMath::Cos(dPhi)));
      hMT[CutIndex]->Fill(MT,EvWeight);
      hmet_dPhi[CutIndex]->Fill(MetV.Pt(),dPhi,EvWeight);
      hmet_MT[CutIndex]->Fill(MetV.Pt(),MT,EvWeight);
      hMT_dPhi[CutIndex]->Fill(MT,dPhi,EvWeight);
     //cout<<"  "<<dPhi<<"  "<<fabs(WBos.DeltaPhi(MuV.at(0)))<<endl;
  
      for (unsigned int im=0;im<tree_.muon_count;im++){
      if (tree_.muon_pt[im] == MuV.at(0).Pt()) 
	      
	{hmu_miniISO[CutIndex]->Fill(tree_.muon_miniISO[im],EvWeight);
   
	Float_t neutralIso = 
	  tree_.muon_neutralHadIso[im] + 
	  tree_.muon_photonIso[im] - 
	  0.5*tree_.muon_puIso[im];
	neutralIso = TMath::Max(Float_t(0),neutralIso); 
	Float_t absIso = tree_.muon_chargedHadIso[im] + neutralIso;
	Float_t relIso = absIso/tree_.muon_pt[im];
        hmu_relISO[CutIndex]->Fill(relIso,EvWeight);
			}
   		}
    }


    Float_t sumpT=0;
    int bjets=0;
    
    if (JetsV.size()>0){
    for (unsigned int ij=0;ij<JetsV.size();ij++){
         sumpT+=JetsV.at(ij).Pt();
         Float_t dPhiJ=dPhiFrom2P( JetsV.at(ij).Px(), JetsV.at(ij).Py(), MetV.Px(),  MetV.Py() );
    
     hdPhiJMET[CutIndex]->Fill(dPhiJ,EvWeight);
	 
     for (unsigned int ib = 0; ib <tree_.pfjet_count;ib++){
        if (tree_.pfjet_pt[ib] == JetsV.at(ij).Pt()  &&  tree_.pfjet_btag[ib][6]  > 0.814) bjets++;
      //  cout<<tree_.pfjet_pt[ib] <<"  "<<JetsV.at(ij).Pt()<<"  "<< tree_.pfjet_btag[ib][6]  <<endl;
     }
      hnBJet[CutIndex]->Fill(bjets,EvWeight);
	    
		    }
     hHT[CutIndex]->Fill(sumpT,EvWeight);
    }
      hMET[CutIndex]->Fill(MetV.Pt(),EvWeight);

		                      //cout<<" pfjet_b "<<ib<<"  "<<analysisTree.pfjet_btag[ib][6]<<endl;
				      //          }
				      //                    if (btagged) continue;


     // hMT[CutIndex]->Fill(MT,EvWeight);
    //  hDZeta->Fill(DZeta,EvWeight);

   /*
    hST[CutIndex]->Fill(Obj.ST,EvWeight);
    hnTop[CutIndex]->Fill(Obj.nTopTagJetGood,EvWeight);
    hdPhiJMET[CutIndex]->Fill(Obj.minDelPhiJMet,EvWeight);
    if(Obj.nTopTagJetGood>0)hToppT[CutIndex]->Fill(Obj.goodTopTagJet[0].Pt(),EvWeight);
    hnW[CutIndex]->Fill(Obj.nWTagJetGood,EvWeight);
    if(Obj.nWTagJetGood>0){
    hWTagpT[CutIndex]->Fill(Obj.goodWTagJet[0].Pt(),EvWeight);
    hWTagMass[CutIndex]->Fill(Obj.goodWTagJet[0].M(),EvWeight);

    }
    hnWmass[CutIndex]->Fill(Obj.nWmassTagJetGood,EvWeight);
    if(Obj.nWmassTagJetGood>0){
     hWmassTagpT[CutIndex]->Fill(Obj.goodWmassTagJet[0].Pt(),EvWeight);
    hWmassTagMass[CutIndex]->Fill(Obj.goodWmassTagJet[0].M(),EvWeight);
   }*/
	
}









