/*
 * MyElectron.h
 *
 *  Created on: Oct 21, 2016
 *      Author: csander
 */

#ifndef MYELECTRON_H_
#define MYELECTRON_H_

#include <TLorentzVector.h>

class MyElectron: public TLorentzVector {

    public:

        MyElectron();
        MyElectron(double pt, double eta, double phi) {
            SetPtEtaPhiM(pt, eta, phi, 511.e-6);
            signal = false;
            OR = true;
            charge = 0;
        };
        virtual ~MyElectron();

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

#endif /* MYELECTRON_H_ */
