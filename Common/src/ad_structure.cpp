/*!
 * \file ad_structure.cpp
 * \brief Main subroutines for the algorithmic differentiation (AD) structure.
 * \author T. Albring
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/datatype_structure.hpp"

namespace AD {
#ifdef CODI_REVERSE_TYPE
  /*--- Initialization of the global variables ---*/

  int adjointVectorPosition = 0;

  std::vector<su2double::GradientData> inputValues;
  std::vector<su2double::GradientData> localInputValues;
  std::vector<su2double*> localOutputValues;

  su2double::TapeType& globalTape = su2double::getGlobalTape();
  su2double::TapeType::Position StartPosition, EndPosition;

  bool Status = false;
  bool PreaccActive = false;
  bool PreaccEnabled = true;
  
  codi::PreaccumulationHelper<su2double> PreaccHelper;
  
  ExtFuncHelper* FuncHelper;

#endif
}


/*--- If we compile under OSX we have to overload some of the operators for
 *   complex numbers to avoid the use of the standard operators
 *  (they use a lot of functions that are only defined for doubles) ---*/

#ifdef __APPLE__

namespace std{
  template<>
  su2double abs(const complex<su2double>& x){
    return sqrt(x.real()*x.real() + x.imag()*x.imag());
  }

  template<>
  complex<su2double> operator/(const complex<su2double>& x, const complex<su2double>& y){

    su2double d    = (y.real()*y.real() + y.imag()*y.imag());
    su2double real = (x.real()*y.real() + x.imag()*y.imag())/d;
    su2double imag = (x.imag()*y.real() - x.real()*y.imag())/d;

    return complex<su2double>(real, imag);

  }

  template<>
  complex<su2double> operator*(const complex<su2double>& x, const complex<su2double>& y){

    su2double real = (x.real()*y.real() - x.imag()*y.imag());
    su2double imag = (x.imag()*y.real() + x.real()*y.imag());

    return complex<su2double>(real, imag);

  }
}
#endif
