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
            jvt = -1;
        };
        virtual ~MyJet();

        void SetJVT(int j) {
            jvt = j;
        };

        double GetJVT() {
            return jvt;
        };

        void SetBTag(bool b) {
            btag = b;
        };

        void SetJetID(bool id) {
            good = id;
        };

        bool IsPU(double cut) {
            return (jvt < cut);
        };

        bool IsNoPU(double cut) {
            return (jvt > cut);
        };

        bool IsGood() {
            return good;
        };

        bool IsBad() {
            return !good;
        };

        bool IsB() {
            return btag;
        };

    private:

        bool btag;
        bool good;
        float jvt;

};

#endif /* MYJET_H_ */
