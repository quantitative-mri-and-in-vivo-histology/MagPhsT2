/** @file epg.h
 *  @brief Implementation of Extended Phase Graph Signal Simulations
 * 
 * Author: tony.stoecker@dzne.de
 * Date  : May 20, 2012
 */

#ifndef EPG_H
#define EPG_H

#include <iostream>
#include <vector>
using namespace std;
#include "math.h"



#define EPG_PI 3.14159265358979
#define EPG_EPSILON  1e-16
#define EPG_TOL  0.0001
#define EPG_FA_TOL  0.01
#define EPG_MAXITER  100
#define EPG_MAXSIZE  99999


class EPG {

 
 public:

       /** @brief Default Constructor */
       EPG (const double &M0 = 1.0, const double &T1 = 1000.0, const double &T2 = 100.0) ;
      ~EPG (        ) ;      /**< @brief destructor                 */
       EPG       (const EPG &epg); /**< @brief copy constructor           */
  void Copy      (const EPG &epg); /**< @brief copy without delete        */
  EPG& operator= (const EPG &epg); /**< @brief assignment (deletes first) */
  bool operator==(const EPG &epg); /**< @brief compare                    */
  bool operator!=(const EPG &epg); /**< @brief compare                    */

 /* EPG state creation and initialization */
 /** @brief Reset EPG to equilibrium, possibly reassign new parameters */
 void SetParameters ( const double &M0, const double &T1, const double &T2 );	/**< @brief assign parameters (sets equilibrium)  */
 void DeleteStates  (  );	/**< @brief clear state vectors */

 /* get parameters of this epg */
 int     GetStep() const { return m_step; };	/**< @brief Return current step of this EPG  */
 double  GetM0  () const { return m_M0;   };	/**< @brief Return equlibirium magnetization */
 double  GetT1  () const { return m_T1;   };	/**< @brief Return longitudinal relaxation   */
 double  GetT2  () const { return m_T2;   };	/**< @brief Return transverse relaxation     */

 /* other attributes of this epg */
 bool    GetVerbose(             ) const { return m_verbose; };
 void    SetVerbose(bool val=true)       { m_verbose=val;    };

 /* State manipulation */
 void Equilibrium    (       ); /**< @brief Reset EPG to equilibrium          */
 void NullTransverse (       ); /**< @brief Set all transverse states to zero */
 void SetStep        (int val); /**< @brief Set step counter of this EPG      */

 /* get state components of order num */
 void   GetFaState ( double* real, double* imag, const int &num = 0 ) const ; /**< @brief get complex transverse   state after last pulse  */
 void   GetZaState ( double* real, double* imag, const int &num = 0 ) const ; /**< @brief get complex longitudinal state after last pulse  */
 void   GetFbState ( double* real, double* imag, const int &num = 0 ) const ; /**< @brief get complex transverse   state before next pulse */
 void   GetZbState ( double* real, double* imag, const int &num = 0 ) const ; /**< @brief get complex longitudinal state before next pulse */

 /** @brief returns transverse magnitude of state num after last RF pulse; longitudinal if transverse = false */
 double GetMagA  ( const int &num = 0 , const bool &transverse = true ) const;
 /** @brief returns transverse magnitude of state num before next RF pulse; longitudinal if transverse = false */
 double GetMagB  ( const int &num = 0 , const bool &transverse = true ) const;

 /** @brief returns real part of transverse magnetization for state num after last RF pulse */
 double GetReA  ( const int &num = 0 ) const;
 /** @brief returns imag part of transverse magnetization for state num after last RF pulse */
 double GetImA  ( const int &num = 0 ) const;
 /** @brief returns real part of transverse magnetization for state num before next RF pulse */
 double GetReB  ( const int &num = 0 ) const;
 /** @brief returns imag part of transverse magnetization for state num before next RF pulse */
 double GetImB  ( const int &num = 0 ) const;

 /** @brief returns transverse magnitude of next state, i.e. after the next RF pulse (fa,phi); longitudinal if transverse = false */
 double GetNextMagA       ( const double &fa, const double &phi,                   const int &num = 0, const bool &transverse = true  ) ;
 

 /** @brief returns phase of last RF pulse */
 double GetPhase () {return m_phase;}

 /* stepping forward in EPG: RF pulse plus TR intervall */
 void Step     ( const double &fa, const double &phi, const double &TR ) ;			/**< @brief single step forward. input: RF pulse plus TR intervall */
 void Steps    ( const double &fa, const double &phi, const double &TR, const int &steps ) ;	/**< @brief multiple steps, constant flipangle and constant phase  */
 void Steps    ( const double* fa, const double &phi, const double &TR, const int &steps ) ;	/**< @brief multiple steps, variable flipangle and constant phase  */
 void Steps    ( const double* fa, const double* phi, const double &TR, const int &steps ) ;	/**< @brief multiple steps, variable flipangle and variable phase  */

// Overloading not working in Cython correctly, seems like compiler thinks Steps(double, double, double, int) and Steps(double*, double, double, int) have the same signature
// This is a workaround. Additionally, this function allows to retrieve the transverse magnitude of state 0 after every RF pulse
vector<double> GetMagTrain    ( const vector<double> &fa, const vector<double> &phi, const double &TR) ;	/**< @brief multiple steps, variable flipangle and variable phase  */

 /** @brief stepping until steady state is reached. Qphi = quadratic phase incr. returns: number of steps needed */
 int StepsToSS ( const double &fa, const double &Qphi, const double &TR, const double &tol = EPG_TOL ) ;

 /** @brief find train of flip angles to obtain a target signal shaping */
 bool FindFlipAngleTrain     (const int &length, double* fa, const double* Ftarget, const double &TR, const double &reduce=0.0, const int &num =0, const double &tol = EPG_TOL );
 
/** @brief find optimal flip angle using Brent's algorithm */
 double FindFlipAngle     ( const double &ax, const double &bx, const double &Ftarget, const double &phi,                   const int &num, const double &tol);

 protected:
 void   Rotation  ( const double &fa, const double &phi, const int &first =0, int last =-1) ; /**< @brief rotation of states  */
 void   Evolution ( const double &TR                                                      ) ; /**< @brief evolution of states */
 
  //private members
  int      m_step;		/**< @brief current step */
  bool     m_verbose;		/**< @brief Equlibirium magnetization */
  double   m_M0;		/**< @brief Equlibirium magnetization */
  double   m_T1;		/**< @brief Longitudinal relaxation   */
  double   m_T2;		      /**< @brief Transverse relaxation     */
  double   m_phase;		    /**< @brief phase of last RF pulse (might be needed for phase locking)  */
  vector<double>  m_FaRe;	/**< @brief Real part of transverse magnetization after the last pulse    */
  vector<double>  m_FaIm;	/**< @brief Imag part of transverse magnetization after the last pulse    */
  vector<double>  m_FbRe;	/**< @brief Real part of transverse magnetization before the next pulse   */
  vector<double>  m_FbIm;	/**< @brief Imag part of transverse magnetization before the next pulse   */
  vector<double>  m_ZaRe;	/**< @brief Real part of longitudinal magnetization after the last pulse  */
  vector<double>  m_ZaIm;	/**< @brief Imag part of longitudinal magnetization after the last pulse  */
  vector<double>  m_ZbRe;	/**< @brief Real part of longitudinal magnetization before the next pulse */
  vector<double>  m_ZbIm;	/**< @brief Imag part of longitudinal magnetization before the next pulse */
};

#endif

