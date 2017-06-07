/*
 * MyPhoton.h
 *
 *  Created on: Oct 21, 2016
 *      Author: csander
 */

#ifndef MYPHOTON_H_
#define MYPHOTON_H_

#include <TLorentzVector.h>

class MyPhoton: public TLorentzVector {

    public:

        MyPhoton();
        MyPhoton(double pt, double eta, double phi) {
            SetPtEtaPhiM(pt, eta, phi, 0.);
            signal = false;
            OR = true;
        };
        virtual ~MyPhoton();

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

    private:

        bool signal;
        bool OR;

};

#endif /* MYPHOTON_H_ */
