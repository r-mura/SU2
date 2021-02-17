/*!
 * \file CFEMStandardPrismVolumeSol.hpp
 * \brief Class for the FEM prism standard element for the grid.
 *        The functions are in the <i>CFEMStandardPrismVolumeSol.cpp</i> file.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "CFEMStandardPrismBase.hpp"

/*!
 * \class CFEMStandardPrismVolumeSol
 * \brief Class which defines the variables and methods for the
 *        prism standard element for the volume solution.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
 */
class CFEMStandardPrismVolumeSol final: public CFEMStandardPrismBase {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardPrismVolumeSol() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly       - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact  - Polynomial order that must be integrated exactly
   *                              by the integration rule.
   * \param[in] val_locGridDOFs - Location of the grid DOFs (LGL or Equidistant).
   * \param[in] val_nVar        - Number of variables in the jitted gemm calls.
   * \param[in] val_useLumpedMM - Whether or not the lumped mass matrix is used.
   */
  CFEMStandardPrismVolumeSol(const unsigned short val_nPoly,
                             const unsigned short val_orderExact,
                             const unsigned short val_locGridDOFs,
                             const unsigned short val_nVar,
                             const bool           val_useLumpedMM);

  /*!
   * \brief Destructor.
   */
  ~CFEMStandardPrismVolumeSol();

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, that makes available the correction factor for the inviscid
   *        spectral radius.
   * \return The correction factor for the inviscid spectral radius.
   */
  inline passivedouble GetFactorInviscidSpectralRadius(void) const override {return factInviscidRad;}

  /*!
   * \brief Function, that makes available the correction factor for the viscous
   *        spectral radius.
   * \return The correction factor for the viscous spectral radius.
   */
  inline passivedouble GetFactorViscousSpectralRadius(void) const override {return factViscousRad;}

  /*!
   * \brief Function that makes available the value of the first (constant)
   *        basis function of this element.
   * \return - The value of the first (constant) basis function.
   */
  inline passivedouble ValBasis0(void) override {return legBasisSolDOFs(0,0);}

  /*-----------------------------------------------------------------------------------*/
  /*---                     Public member functions.                                ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which determines the basis functions for the given parametric
   *        coordinates.
   * \param[in]  parCoor  - Double vector that contains the parametric coordinates
   *                        for which the basis functions must be determined.
   * \param[out] matBasis - Matrix that contains the values of the basis functions
   *                        in the given parametric coordinates.
   */
  void BasisFunctionsInPoints(const vector<vector<passivedouble> > &parCoor,
                              ColMajorMatrix<passivedouble>        &matBasis) override;

  /*!
   * \brief Function, that converts the modal form of the DOFs to the nodal form.
   * \param[in,out] solDOFs - On entry it contains the modal solution in the DOFs,
   *                          on exit it contains the nodal solution.
   */
  void ModalToNodal(ColMajorMatrix<su2double> &solDOFs) override;

  /*!
   * \brief Function, that converts the nodal form of the DOFs to the modal form.
   * \param[in,out] solDOFs - On entry it contains the nodal solution in the DOFs,
   *                          on exit it contains the modal solution.
   */
  void NodalToModal(ColMajorMatrix<su2double> &solDOFs) override;

  /*!
   * \brief Function that computes the gradients of the solution in integration points.
   * \param[in]  matSolDOF     - Matrix that contains the modal solution DOFs.
   * \param[out] matGradSolInt - Vector of matrices the contains the gradients of the
   *                             solution in the integration points.
   */
  void GradSolIntPoints(ColMajorMatrix<su2double>          &matSolDOF,
                        vector<ColMajorMatrix<su2double> > &matGradSolInt) override;

  /*!
   * \brief Function that computes the solution in integration points.
   * \param[in]  matSolDOF - Matrix that contains the modal solution DOFs.
   * \param[out] matSolInt - Matrix that contains the solution in the integration points.
   */
  void SolIntPoints(ColMajorMatrix<su2double> &matSolDOF,
                    ColMajorMatrix<su2double> &matSolInt) override;

private:
  passivedouble factInviscidRad = -1.0;  /*!< \brief Correction factor for the inviscid spectral radius
                                                     for the high order element. */
  passivedouble factViscousRad = -1.0;   /*!< \brief Correction factor for the viscous spectral radius
                                                     for the high order element. */

  vector<passivedouble> rLineSolDOFs;     /*!< \brief 1D parametric coordinates of the nodal solution DOFs. */
  vector<passivedouble> rTriangleSolDOFs; /*!< \brief Parametric r-coordinates of the triangle nodal solution DOFs. */
  vector<passivedouble> sTriangleSolDOFs; /*!< \brief Parametric s-coordinates of the triangle nodal solution DOFs. */

  ColMajorMatrix<passivedouble> legBasisInt;             /*!< \brief The values of the Legendre basis functions
                                                                     in the integration points. */
  vector<ColMajorMatrix<passivedouble> > derLegBasisInt; /*!< \brief The values of the derivatives of the Legendre
                                                                     basis functions in the integration points.
                                                                     It is a vector, because there are derivatives
                                                                     in three directions. */
  vector<ColMajorMatrix<passivedouble> > hesLegBasisInt; /*!< \brief The values of the 2nd derivatives of the Legendre
                                                                     basis functions in the integration points.
                                                                     It is a vector, because there are 6 2nd derivatives
                                                                     in three space dimensions. */

  ColMajorMatrix<passivedouble> legBasisSolDOFs;             /*!< \brief The values of the Legendre basis functions
                                                                         in the solution DOFs. */
  ColMajorMatrix<passivedouble> legBasisSolDOFsInv;          /*!< \brief The inverse of legBasisLineSolDOFs. */
  
  vector<ColMajorMatrix<passivedouble> > derLegBasisSolDOFs; /*!< \brief The values of the derivatives of the
                                                                         Legendre basis functions in the solution
                                                                         DOFs. It is a vector, because there
                                                                         are derivatives in three directions. */

  void *jitterDOFs2Int = nullptr;      /*!< \brief Pointer to the data for the jitted gemm function
                                                   to compute data in the integration points. */
  dgemm_jit_kernel_t gemmDOFs2Int;     /*!< \brief Pointer to the function to carry out jitterDOFs2Int. */

  void *jitterDOFs2SolDOFs = nullptr;  /*!< \brief Pointer to the data for the jitted gemm function
                                                   to compute data in the solution DOFs. */
  dgemm_jit_kernel_t gemmDOFs2SolDOFs; /*!< \brief Pointer to the function to carry out jitterDOFs2SolDOFs. */
};
