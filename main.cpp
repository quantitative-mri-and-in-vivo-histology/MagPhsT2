#include <iostream>
#include <stdio.h>
#include "epg.h"
#include "external/cnpy/cnpy.h"
#include <chrono>

using namespace std::chrono;


int main() {

    unsigned long nPulses = 64;
    double dAlpha = 40.;
    double dPhi = 0.0;
    double dTR = 5.0;
    double dM0 = 1.0;
    double dT1 = 800.0;
    double dT2 = 50.0;

    std::vector<double> vSignal(nPulses, 0.0);
    std::vector<double> vFlipArray(nPulses, dAlpha);

    // Create EPG object (M0, T1, T2)
    EPG epg(dM0, dT1, dT2);

    // Output tissue parameters
    printf("\n\nTissue parameters\n");
    printf("=================\n");
    printf("M0/T1/T2: %.1f/%.1f/%.1f [ms]\n\n", epg.GetM0(), epg.GetT1(), epg.GetT2());

    printf("Simple example\n");
    printf("==============\n");
    // Apply RF pulse (flipangle dAlpha, phase dPhi) and evolve one repetition time TR
    epg.Step(dAlpha, dPhi, dTR);
    printf("Step (# applied RF pulses): %d\n", epg.GetStep());
    printf("Mx: %.2f\n", epg.GetMagA());
    printf("Mz: %.2f\n\n", epg.GetMagA(0, false));

    // Apply RF pulse again
    epg.Step(dAlpha, dPhi, dTR);
    printf("Step (# applied RF pulses): %d\n", epg.GetStep());
    printf("Mx: %.2f\n", epg.GetMagA());
    printf("Mz: %.2f\n\n", epg.GetMagA(0, false));

    // Spoil transverse magnetization
    epg.NullTransverse();
    printf("Spoiling of transverse magnetization\n");
    printf("Mx: %.2f\n", epg.GetMagA());
    printf("Mz: %.2f\n\n", epg.GetMagA(0, false));

    // Back to equilibrium
    epg.Equilibrium();

    // RF-spoiled GRE: Steps (number of RF pulses) needed to reach steady state using quadratic phase increment of 50deg
    printf("RF-spoiled GRE\n");
    printf("==============\n");
    printf("Steps to steady state: %.d\n\n", epg.StepsToSS(dAlpha, 50.0, dTR));

    printf("Non-spoiled GRE\n");
    printf("===============\n");

    // Write signal evolution of flipangle train (here: still constant FA) in vector.
    epg.Equilibrium();
    for (int k=0;k<nPulses;k++) {
        epg.Step(vFlipArray[k], 0.0, dTR);
        vSignal[k] = epg.GetMagA(0);
        if(k%10 == 0)
        std::cout << "Signal after pulse " << (k+1) << ": " << vSignal[k] << std::endl;
    }

    // Find flipangle train and measure time
    auto start = high_resolution_clock::now();

    nPulses = 256;
    std::vector<double> vSignal2(nPulses, 1.0);
    std::vector<double> vFlipArray2(nPulses, 0.0);
    epg.Equilibrium();
    epg.FindFlipAngleTrain(nPulses,&vFlipArray2[0], &vSignal2[0], dTR);

    //Save signal array to npy file
    cnpy::npy_save("/tmp/gre.npy",&vFlipArray2[0],{nPulses},"w");

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    cout << "\nFinding and writing flipangle train took " << duration.count() << " milliseconds" << endl;


    return 0;
}
