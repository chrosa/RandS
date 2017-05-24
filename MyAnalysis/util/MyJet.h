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
            fjvt = -1;
            ntracks = 0;
            trackWidth = 0;
            sumpt = 0;
        };
        virtual ~MyJet();

        void SetJVT(double j) {
            jvt = j;
        };

        void SetFJVT(double j) {
            fjvt = j;
        };

        double GetJVT() {
            return jvt;
        };

        double GetFJVT() {
            return fjvt;
        };

        bool IsPU(double cut) {
            return (jvt < cut);
        };

        bool IsNoPU(double cut) {
            return (jvt > cut);
        };

        bool IsFPU(double cut) {
            return (fjvt < cut);
        };

        bool IsNoFPU(double cut) {
            return (fjvt > cut);
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

        bool IsTau() {
            return !isTau;
        };

    private:

        bool btag;
        bool good;
        bool isPhoton;
        bool isTau;
        float jvt;
        float fjvt;
        float sumpt;
        float trackWidth;
        unsigned short ntracks;
        

};

#endif /* MYJET_H_ */
