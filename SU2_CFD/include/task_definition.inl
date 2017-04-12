/*!
 * \file task_definition.inl
 * \brief In-Line subroutines of the <i>task_definition.hpp</i> file.
 * \author E. van der Weide, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
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

#pragma once

inline CTaskDefinition::CTaskDefinition(void) {
  task                = NO_TASK;
  timeLevel           = 0;
  nIndMustBeCompleted = 0;
  for(int i=0; i<4; ++i) indMustBeCompleted[i] = -1;

  intPointADER = 0;
  secondPartTimeIntADER = false;
}

inline CTaskDefinition::CTaskDefinition(SOLVER_TASK    val_task,
                                        unsigned short val_timeLevel,
                                        int            val_ind0MustBeCompleted,
                                        int            val_ind1MustBeCompleted,
                                        int            val_ind2MustBeCompleted,
                                        int            val_ind3MustBeCompleted) {
  task                  = val_task;
  timeLevel             = val_timeLevel;
  indMustBeCompleted[0] = val_ind0MustBeCompleted;
  indMustBeCompleted[1] = val_ind1MustBeCompleted;
  indMustBeCompleted[2] = val_ind2MustBeCompleted;
  indMustBeCompleted[3] = val_ind3MustBeCompleted;

  for(nIndMustBeCompleted=0; nIndMustBeCompleted<4; ++nIndMustBeCompleted) {
    if(indMustBeCompleted[nIndMustBeCompleted] < 0) break;
  }

  intPointADER = 0;
  secondPartTimeIntADER = false;
}

inline CTaskDefinition::~CTaskDefinition(void) {}

inline CTaskDefinition::CTaskDefinition(const CTaskDefinition &other){Copy(other);}

inline CTaskDefinition& CTaskDefinition::operator=(const CTaskDefinition &other){Copy(other); return (*this);}

inline void CTaskDefinition::Copy(const CTaskDefinition &other) {
  task                  = other.task;
  timeLevel             = other.timeLevel;
  intPointADER          = other.intPointADER;
  secondPartTimeIntADER = other.secondPartTimeIntADER;
  nIndMustBeCompleted   = other.nIndMustBeCompleted;

  for(int i=0; i<4; ++i)
    indMustBeCompleted[i] = other.indMustBeCompleted[i];
}
