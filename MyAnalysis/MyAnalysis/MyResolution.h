#ifndef MyAnalysis_MyResolution_H
#define MyAnalysis_MyResolution_H

#include <EventLoop/Algorithm.h>
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "PileupReweighting/PileupReweightingTool.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "SUSYTools/SUSYToolsDict.h"
#include "SUSYTools/ISUSYObjDef_xAODTool.h"
#include <TH1.h>

#include <memory>
#include <string>

class MyResolution : public EL::Algorithm
{
    private:
        double m_VetoCone;
        double m_AddActivityCone;
        double m_MatchingCone;
        double m_RelGenActivityVeto;
        double m_RelRecoActivityVeto;
		std::string prw_file_;
		std::string ilumicalc_file_;

        // put your configuration variables here as public variables.
        // that way they can be set directly from CINT and python.
    public:
        // float cutValue;
        int m_eventCounter; //!
        int m_numCleanEvents; //!
        GoodRunsListSelectionTool *m_grl; //!
        JetCleaningTool *m_jetCleaning; //!
		ST::SUSYObjDef_xAOD *objTool; //!

        // variables that don't get filled at submission time should be
        // protected from being send from the submission node to the worker
        // node (done by the //!)

        std::vector<double> PtBinEdges; //!
        std::vector<double> EtaBinEdges; //!

        // Resize histo vectors
        void ResizeHistoVector(std::vector<std::vector<TH1F*> > &histoVector);

        std::vector< std::vector<TH1F*> > PtResolution_tot; //!
        std::vector< std::vector<TH1F*> > EtaResolution_tot; //!
        std::vector< std::vector<TH1F*> > PhiResolution_tot; //!

        std::vector< std::vector<TH1F*> > PtResolution_LF; //!
        std::vector< std::vector<TH1F*> > EtaResolution_LF; //!
        std::vector< std::vector<TH1F*> > PhiResolution_LF; //!

        std::vector< std::vector<TH1F*> > PtResolution_HF; //!
        std::vector< std::vector<TH1F*> > EtaResolution_HF; //!
        std::vector< std::vector<TH1F*> > PhiResolution_HF; //!

        std::vector< std::vector<TH1F*> > PtResolution_b; //!
        std::vector< std::vector<TH1F*> > EtaResolution_b; //!
        std::vector< std::vector<TH1F*> > PhiResolution_b; //!

        std::vector< std::vector<TH1F*> > PtResolution_nob; //!
        std::vector< std::vector<TH1F*> > EtaResolution_nob; //!
        std::vector< std::vector<TH1F*> > PhiResolution_nob; //!

        std::vector<TH1F*> NReco_tot; //!
        std::vector<TH1F*> NReco_b; //!
        std::vector<TH1F*> NReco_nob; //!

        std::vector<TH1F*> NGen_tot; //!
        std::vector<TH1F*> NGen_b; //!
        std::vector<TH1F*> NGen_nob; //!

    public:
        // Tree *myTree; //!
        // TH1 *myHist; //!

        // this is a standard constructor
        MyResolution ();

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

        std::string GetHistName(unsigned int i_pt, unsigned int i_eta, std::string s1, std::string s2);
        std::string GetHistNameEta(unsigned int i_pt, std::string s1);
        unsigned int GetPtBin(double pt);
        unsigned int GetEtaBin(double eta);

        // this is needed to distribute the algorithm to the workers
        ClassDef(MyResolution, 1);
};

#endif
