# Using a .pxd file gives us a separate namespace for
# the C++ declarations. Using a .pxd file also allows
# us to reuse the declaration in multiple .pyx modules.
from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "epg.h":
   
    cppclass EPG:
        EPG(double, double, double) except +

        void SetParameters(double, double, double)
        void DeleteStates()

        int  GetStep()
        double GetM0()
        double GetT1()
        double GetT2()

        bool GetVerbose()
        void SetVerbose(bool)
        void Equilibrium()
        void NullTransverse()
        void SetStep(int)

        double GetMagA(int,bool)
        double GetMagB(int,bool)

        double GetReA(int)
        double GetImA(int)
        double GetReB(int)
        double GetImB(int)

        double GetNextMagA(double, double, int, bool)

        double GetPhase()

        void Step(double, double, double)
        void Steps(double, double, double, int)
        vector[double] GetMagTrain(vector[double] fa, vector[double] phi, double)

        int StepsToSS(double, double, double, double) ;

        bool FindFlipAngleTrain(int, double* fa, double* Ftarget, double, double, int, double)
