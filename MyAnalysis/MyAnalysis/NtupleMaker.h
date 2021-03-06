#ifndef MyAnalysis_NtupleMaker_H
#define MyAnalysis_NtupleMaker_H

#include <EventLoop/Algorithm.h>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"

#include "GoodRunsLists/GoodRunsListSelectionTool.h"

#include "SUSYTools/SUSYToolsDict.h"
#include "SUSYTools/ISUSYObjDef_xAODTool.h"
#include "SUSYTools/SUSYCrossSection.h"

#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"

#include <memory>
#include <vector>
#include <string>
#include <map>


class NtupleMaker : public EL::Algorithm
{
    // put your configuration variables here as public variables.
    // that way they can be set directly from CINT and python.

    public:

        int m_eventCounter; //!
        int m_numCleanEvents; //!
        bool AddMuToJets_; //!
        bool isMC;
        bool doPRW_; //!
        
        int debug_;
        std::string outputfile_;
        
        bool calculateTrueMET_; //!
        float rebalancedJetPt_; //!

        // GRL tool
        GoodRunsListSelectionTool *m_grl; //!

   		ST::SUSYObjDef_xAOD *objTool; //!
   		SUSY::CrossSectionDB *my_XsecDB; //!
		std::string prw_file_;
		std::string ilumicalc_file_;

        // variables that don't get filled at submission time should be
        // protected from being send from the submission node to the worker
        // node (done by the //!)

        // this is a standard constructor
        NtupleMaker ();

    private:

        typedef TLorentzVector LorentzVector;
        typedef std::vector<std::string>::const_iterator StrIter;
        typedef struct myJet {
            LorentzVector momentum;
            bool btag;
            bool good;
            float jvt;
            bool fjvt;
            float tw;
            int ntracks;
            int vtx;
            float sumpt;
            bool OR;
			float FracSamplingMax;
			float HECFrac;
			float EMFrac;

        } myJet;

        typedef struct myPhoton {
            LorentzVector momentum;
            bool isSignal;
            bool OR;
        } myPhoton;

        typedef struct myElectron {
            LorentzVector momentum;
            int charge;
            bool isSignal;
            bool OR;
        } myElectron;

        typedef struct myMuon {
            LorentzVector momentum;
            bool isSignal;
            bool isBad;
            int charge;
            bool OR;
        } myMuon;

        typedef struct myTau {
            LorentzVector momentum;
            int charge;
            bool isSignal;
            bool OR;
        } myTau;

        // these are the functions inherited from Algorithm
        virtual EL::StatusCode setupJob (EL::Job& job);
        virtual EL::StatusCode fileExecute ();
        virtual EL::StatusCode histInitialize ();
        virtual EL::StatusCode changeInput (bool firstFile);
        virtual EL::StatusCode initialize ();
        virtual EL::StatusCode execute ();
        virtual EL::StatusCode postExecute ();
        virtual EL::StatusCode finalize ();
        virtual EL::StatusCode histFinalize ();

        TTree *EventTree; //!
        UShort_t NVtx_; //!
        Float_t weight_; //!
		UInt_t dsid_;  //!
		UInt_t evtno_;  //!
		bool pvtx_;  //!
		bool j400triggered_;  //!
		bool xe70triggered_;  //!
		bool xe90triggered_;  //!
		bool xe110triggered_;  //!

        std::vector<Float_t>  JetPt_; //!
		std::vector<Float_t> * JetPt_n = &JetPt_; //!
        std::vector<Float_t>  JetEta_; //!
		std::vector<Float_t> * JetEta_n = &JetEta_; //!
        std::vector<Float_t>  JetPhi_; //!
		std::vector<Float_t> * JetPhi_n = &JetPhi_; //!
        std::vector<Float_t>  JetM_; //!
		std::vector<Float_t> * JetM_n = &JetM_; //!
        std::vector<bool>  JetBtag_; //!
		std::vector<bool> * JetBtag_n = &JetBtag_; //!
        std::vector<Float_t>  JetJVT_; //!
		std::vector<Float_t> * JetJVT_n = &JetJVT_; //!
        std::vector<bool>  JetFJVT_; //!
		std::vector<bool> * JetFJVT_n = &JetFJVT_; //!
        std::vector<bool>  JetGood_; //!
		std::vector<bool> * JetGood_n = &JetGood_; //!
        std::vector<bool>  JetPassOR_; //!
		std::vector<bool> * JetPassOR_n = &JetPassOR_; //!
        std::vector<UShort_t>  HighestJVFVtx_; //!
		std::vector<UShort_t> * HighestJVFVtx_n = &HighestJVFVtx_; //!
        std::vector<UShort_t>  JetNTracks_; //!
		std::vector<UShort_t> * JetNTracks_n = &JetNTracks_; //!
        std::vector<Float_t>  JetSumPtTracks_; //!
		std::vector<Float_t> * JetSumPtTracks_n = &JetSumPtTracks_; //!
        std::vector<Float_t>  JetTrackWidth_; //!
		std::vector<Float_t> * JetTrackWidth_n = &JetTrackWidth_; //!
        std::vector<Float_t>  FracSamplingMax_; //!
		std::vector<Float_t> * FracSamplingMax_n = &FracSamplingMax_; //!
        std::vector<Float_t>  HECFrac_; //!
		std::vector<Float_t> * HECFrac_n = &HECFrac_; //!
        std::vector<Float_t>  EMFrac_; //!
		std::vector<Float_t> * EMFrac_n = &EMFrac_; //!

        std::vector<Float_t>  GenJetPt_; //!
		std::vector<Float_t> * GenJetPt_n = &GenJetPt_; //!
        std::vector<Float_t>  GenJetEta_; //!
		std::vector<Float_t> * GenJetEta_n = &GenJetEta_; //!
        std::vector<Float_t>  GenJetPhi_; //!
		std::vector<Float_t> * GenJetPhi_n = &GenJetPhi_; //!
        std::vector<Float_t>  GenJetM_; //!
		std::vector<Float_t> * GenJetM_n = &GenJetM_; //!
        std::vector<bool>  GenJetBtag_; //!
		std::vector<bool> * GenJetBtag_n = &GenJetBtag_; //!

        std::vector<Float_t>  GenJetNoNuPt_; //!
		std::vector<Float_t> * GenJetNoNuPt_n = &GenJetNoNuPt_; //!
        std::vector<Float_t>  GenJetNoNuEta_; //!
		std::vector<Float_t> * GenJetNoNuEta_n = &GenJetNoNuEta_; //!
        std::vector<Float_t>  GenJetNoNuPhi_; //!
		std::vector<Float_t> * GenJetNoNuPhi_n = &GenJetNoNuPhi_; //!
        std::vector<Float_t>  GenJetNoNuM_; //!
		std::vector<Float_t> * GenJetNoNuM_n = &GenJetNoNuM_; //!
        std::vector<bool>  GenJetNoNuBtag_; //!
		std::vector<bool> * GenJetNoNuBtag_n = &GenJetNoNuBtag_; //!

        std::vector<Float_t>  GenJetNoNuMuPt_; //!
		std::vector<Float_t> * GenJetNoNuMuPt_n = &GenJetNoNuMuPt_; //!
        std::vector<Float_t>  GenJetNoNuMuEta_; //!
		std::vector<Float_t> * GenJetNoNuMuEta_n = &GenJetNoNuMuEta_; //!
        std::vector<Float_t>  GenJetNoNuMuPhi_; //!
		std::vector<Float_t> * GenJetNoNuMuPhi_n = &GenJetNoNuMuPhi_; //!
        std::vector<Float_t>  GenJetNoNuMuM_; //!
		std::vector<Float_t> * GenJetNoNuMuM_n = &GenJetNoNuMuM_; //!
        std::vector<bool>  GenJetNoNuMuBtag_; //!
		std::vector<bool> * GenJetNoNuMuBtag_n = &GenJetNoNuMuBtag_; //!

        std::vector<Float_t>  ElePt_; //!
		std::vector<Float_t> * ElePt_n = &ElePt_; //!
        std::vector<Float_t>  EleEta_; //!
		std::vector<Float_t> * EleEta_n = &EleEta_; //!
        std::vector<Float_t>  ElePhi_; //!
		std::vector<Float_t> * ElePhi_n = &ElePhi_; //!
        std::vector<bool>  EleIsSignal_; //!
		std::vector<bool> * EleIsSignal_n = &EleIsSignal_; //!
		std::vector<Int_t> EleCharge_; //!
		std::vector<Int_t> * EleCharge_n = &EleCharge_; //!
		std::vector<bool>  ElePassOR_; //!
		std::vector<bool> * ElePassOR_n = &ElePassOR_; //!


        std::vector<Float_t>  MuonPt_; //!
		std::vector<Float_t> * MuonPt_n = &MuonPt_; //!
        std::vector<Float_t>  MuonEta_; //!
		std::vector<Float_t> * MuonEta_n = &MuonEta_; //!
        std::vector<Float_t>  MuonPhi_; //!
		std::vector<Float_t> * MuonPhi_n = &MuonPhi_; //!
        std::vector<bool>  MuonIsBad_; //!
		std::vector<bool> * MuonIsBad_n = &MuonIsBad_; //!
        std::vector<bool>  MuonIsSignal_; //!
		std::vector<bool> * MuonIsSignal_n = &MuonIsSignal_; //!
		std::vector<Int_t> MuonCharge_; //!
		std::vector<Int_t> * MuonCharge_n = &MuonCharge_; //!
		std::vector<bool>  MuonPassOR_; //!
		std::vector<bool> * MuonPassOR_n = &MuonPassOR_; //!

        std::vector<Float_t>  TauPt_; //!
		std::vector<Float_t> * TauPt_n = &TauPt_; //!
        std::vector<Float_t>  TauEta_; //!
		std::vector<Float_t> * TauEta_n = &TauEta_; //!
        std::vector<Float_t>  TauPhi_; //!
		std::vector<Float_t> * TauPhi_n = &TauPhi_; //!
        std::vector<bool>  TauIsSignal_; //!
		std::vector<bool> * TauIsSignal_n = &TauIsSignal_; //!
		std::vector<Int_t> TauCharge_; //!
		std::vector<Int_t> * TauCharge_n = &TauCharge_; //!
		std::vector<bool>  TauPassOR_; //!
		std::vector<bool> * TauPassOR_n = &TauPassOR_; //!

        std::vector<Float_t>  PhotonPt_; //!
		std::vector<Float_t> * PhotonPt_n = &PhotonPt_; //!
        std::vector<Float_t>  PhotonEta_; //!
		std::vector<Float_t> * PhotonEta_n = &PhotonEta_; //!
        std::vector<Float_t>  PhotonPhi_; //!
		std::vector<Float_t> * PhotonPhi_n = &PhotonPhi_; //!
        std::vector<bool>  PhotonIsSignal_; //!
		std::vector<bool> * PhotonIsSignal_n = &PhotonIsSignal_; //!
		std::vector<bool>  PhotonPassOR_; //!
		std::vector<bool> * PhotonPassOR_n = &PhotonPassOR_; //!

		float MET_pt_;
		float MET_phi_;
		float METjet_pt_;
		float METjet_phi_;
		float METmu_pt_;
		float METmu_phi_;
		float METele_pt_;
		float METele_phi_;
		float METgamma_pt_;
		float METgamma_phi_;
		float METtrack_pt_;
		float METtrack_phi_;

		float GenMET_pt_;
		float GenMET_phi_;
		float TrueMHT_pt_;
		float TrueMHT_phi_;

        // this is needed to distribute the algorithm to the workers
        ClassDef(NtupleMaker, 1);
};

template<typename T>
struct GreaterByPt {
    typedef T first_argument_type;
    typedef T second_argument_type;
    bool operator()( const T & t1, const T & t2 ) const {
        return t1.momentum.Pt() > t2.momentum.Pt();
    }
};

#endif
