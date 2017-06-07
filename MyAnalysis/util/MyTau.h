/*
 * MyTau.h
 *
 *  Created on: Oct 21, 2016
 *      Author: csander
 */

#ifndef MYTAU_H_
#define MYTAU_H_

#include <TLorentzVector.h>

class MyTau: public TLorentzVector {

    public:

        MyTau();
        MyTau(double pt, double eta, double phi) {
            SetPtEtaPhiM(pt, eta, phi, 1.776e0);
            signal = false;
            OR = true;
            charge = 0;
        };
        virtual ~MyTau();

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

        int GetCharge() {
            return charge;
        };

        void SetCharge(int q) {
            charge = q;
        };

    private:

        bool signal;
        bool OR;
        int charge;

};

#endif /* MYTAU_H_ */
