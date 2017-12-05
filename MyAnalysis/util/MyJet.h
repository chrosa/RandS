/*
 * MyJet.h
 *
 *  Created on: Oct 21, 2016
 *      Author: csander
 */

#ifndef MYJET_H_
#define MYJET_H_

#include <TLorentzVector.h>

class MyJet: public TLorentzVector {

    public:

        MyJet();
        MyJet(double pt, double eta, double phi, double m) {
            SetPtEtaPhiM(pt, eta, phi, m);
            btag = false;
            good = true;
            isPhoton = false;
            isTau = false;
            jvt = -1;
            fjvt = false;
            ntracks = 0;
            trackWidth = 0;
            sumpt = 0;
            OR = true;
			vtx = 0;
			FracSamplingMax = 0.5;
			HECFrac = 0.5;
			EMFrac = 0.5;
        };
        virtual ~MyJet();

        void SetJVT(double j) {
            jvt = j;
        };

        void SetFJVT(bool j) {
            fjvt = j;
        };

        double GetJVT() {
            return jvt;
        };

        bool IsFJVT() {
            return fjvt;
        };

        bool PassOR() {
            return OR;
        };

        bool IsPU(double cut) {
            return (jvt < cut);
        };

        bool IsNoPU(double cut) {
            return (jvt > cut);
        };

        void SetNTracks(unsigned short n) {
            ntracks = n;
        };

        unsigned short GetNTracks() {
            return ntracks;
        };

        void SetTrackWidth(double tw) {
            trackWidth = tw;
        };

        double GetTrackWidth() {
            return trackWidth;
        };

        void SetSumPtTracks(double f) {
            sumpt = f;
        };

        double GetSumPtTracks() {
            return sumpt;
        };

        void SetBTag(bool b) {
            btag = b;
        };

        bool IsB() {
            return btag;
        };

        void SetJetID(bool id) {
            good = id;
        };

        bool IsGood() {
            return good;
        };

        bool IsBad() {
            return !good;
        };

        void SetIsPhoton(bool b) {
            isPhoton = b;
        };

        bool IsPhoton() {
            return !isPhoton;
        };

        void SetIsTau(bool b) {
            isTau = b;
        };

        void SetPassOR(bool b) {
            OR = b;
        };

        bool IsTau() {
            return !isTau;
        };

		int GetVtx() {
			return vtx;
		};
		
		void SetVtx(int i) {
			vtx = i;
		};

		float GetFracSamplingMax(){
			return FracSamplingMax;
		};
		
		void SetFracSamplingMax(double f){
			FracSamplingMax = f;
		};

		float GetHECFrac(){
			return HECFrac;
		};
		
		void SetHECFrac(double f){
			HECFrac = f;
		};

		float GetEMFrac(){
			return EMFrac;
		};

		void SetEMFrac(double f){
			EMFrac = f;
		};

    private:

        bool btag;
        bool good;
        bool isPhoton;
        bool isTau;
        bool OR;
        bool fjvt;
        float jvt;
        float sumpt;
        float trackWidth;
        unsigned short ntracks;
		int vtx;
		float FracSamplingMax;
		float HECFrac;
		float EMFrac;

        
};

#endif /* MYJET_H_ */
