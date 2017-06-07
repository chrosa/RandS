/*
 * MyMuon.h
 *
 *  Created on: Oct 21, 2016
 *      Author: csander
 */

#ifndef MYMUON_H_
#define MYMUON_H_

#include <TLorentzVector.h>

class MyMuon: public TLorentzVector {

    public:

        MyMuon();
        MyMuon(double pt, double eta, double phi) {
            SetPtEtaPhiM(pt, eta, phi, 106.e-3);
            OR = true;
            signal = false;
            bad = false;
            charge = 0;
        };
        virtual ~MyMuon();

        bool PassOR() {
            return OR;
        };

        void SetPassOR(bool b) {
            OR = b;
        };

        bool IsSignal() {
            return signal;
        };

        void SetIsSignal(bool b) {
            signal = b;
        };

        bool IsBad() {
            return bad;
        };

        void SetIsBad(bool b) {
            bad = b;
        };

        int GetCharge() {
            return charge;
        };

        void SetCharge(int q) {
            charge = q;
        };


    private:

        bool bad;
        bool signal;
        bool OR;
        int charge;

};

#endif /* MYMUON_H_ */
