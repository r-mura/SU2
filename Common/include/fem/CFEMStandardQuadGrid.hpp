/*!
 * \file CFEMStandardQuadGrid.hpp
 * \brief Class for the FEM quadrilateral standard element for the grid.
 *        The functions are in the <i>CFEMStandardQuadGrid.cpp</i> file.
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

#include "CFEMStandardQuadBase.hpp"

/*!
 * \class CFEMStandardQuadGrid
 * \brief Class which defines the variables and methods for the
 *        quadrilateral standard element for the grid.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
 */
class CFEMStandardQuadGrid final: public CFEMStandardQuadBase {

public:
  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardQuadGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPolyGrid   - Polynomial degree of the grid for this element.
   * \param[in] val_nPolyGrid   - Polynomial degree of the solution for this element.
   * \param[in] val_orderExact  - Polynomial order that must be integrated exactly
   *                              by the integration rule.
   * \param[in] val_locGridDOFs - Location of the grid DOFS, either LGL or equidistant.
   */
  CFEMStandardQuadGrid(const unsigned short val_nPolyGrid,
                       const unsigned short val_nPolySol,
                       const unsigned short val_orderExact,
                       const unsigned short val_locGridDOFs);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardQuadGrid() = default;

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, that returns the number of solution DOFs.
   * \return The number of solution DOFs of the volume element.
   */
  inline unsigned short GetNSolDOFs(void) const override {
    const unsigned short nSol1D = rLineSolDOFs.size();
    return nSol1D*nSol1D;
  }

  /*!
   * \brief Function, that returns the padded number of solution DOFs.
   * \return The padded number of solution DOFs of the volume element.
   */
  inline unsigned short GetNSolDOFsPad(void) const override {
    const unsigned short nSol1D = rLineSolDOFs.size();
    const unsigned short nSol   = nSol1D*nSol1D;
    return PaddedValue(nSol);
  }

  /*-----------------------------------------------------------------------------------*/
  /*---                     Public member functions.                                ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which computes the coordinates in the integration points.
   * \param[in]  notUsed    - Argument present to be consistent with the base class
   *                          function, which is overwritten.
   * \param[in]  matCoorDOF - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matCoorInt - Matrix that contains the coordinates of the integration
   *                          points.
   */
  void CoorIntPoints(const bool                notUsed,
                     ColMajorMatrix<su2double> &matCoorDOF,
                     ColMajorMatrix<su2double> &matCoorInt) override;

  /*!
   * \brief Function, which computes the coordinates of the solution DOFs. These are
   *        only different from the grid DOFs when different polynomials degrees
   *        are used for the grid and solution.
   * \param[in]  matCoorDOF    - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matCoorSolDOF - Matrix that contains the coordinates of the solution DOFs.
   */
  void CoorSolDOFs(ColMajorMatrix<su2double> &matCoorDOF,
                   ColMajorMatrix<su2double> &matCoorSolDOF) override;

  /*!
   * \brief Function, which computes the derivatives of the coordinates in the
   *        integration points.
   * \param[in]  notUsed    - Argument present to be consistent with the base class
   *                          function, which is overwritten.
   * \param[in]  matCoor    - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matDerCoor - Vector of matrices to store the derivatives of the coordinates.
   */
  void DerivativesCoorIntPoints(const bool                         notUsed,
                                ColMajorMatrix<su2double>          &matCoor,
                                vector<ColMajorMatrix<su2double> > &matDerCoor) override;

  /*!
   * \brief Function, which computes the 2nd derivatives of the coordinates in the
   *        integration points.
   * \param[in]  matCoor       - Matrix that contains the coordinates of the grid DOFs
   * \param[out] matDer2ndCoor - Vector of matrices to store the 2nd derivatives of the coordinates.
   */
  void Derivatives2ndCoorIntPoints(ColMajorMatrix<su2double>          &matCoor,
                                   vector<ColMajorMatrix<su2double> > &matDer2ndCoor) override;

  /*!
   * \brief Function, which computes the derivatives of the coordinates in the
   *        solution DOFs.
   * \param[in]  matCoor    - Matrix that contains the coordinates of the grid DOFs.
   * \param[out] matDerCoor - Vector of matrices to store the derivatives of the coordinates.
   */
  void DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                              vector<ColMajorMatrix<su2double> > &matDerCoor) override;

  /*!
   * \brief Function, which evaluates the coordinates and gradients for the given
   *        parametric coordinates.
   * \param[in]  matCoor - Matrix that contains the coordinates of the grid DOFs.
   * \param[in]  par     - Parametric coordinates for which the physical coordinates
   *                       and derivatives must be determined.
   * \param[out] x       - Physical coordinates to be determined.
   * \param[out] dxdpar  - Derivatives of the physical coordinates w.r.t. the parametric
   *                       coordinates.
   */
  void EvalCoorAndGradCoor(ColMajorMatrix<su2double> &matCoor,
                           const su2double           *par,
                           su2double                 x[3],
                           su2double                 dxdpar[3][3]) override;

  /*!
   * \brief Function, which interpolates the parametric coordinate within the
   *        given linear sub-element.
   * \param[in]  subElem - Sub-element in which the coordinate must be interpolated.
   * \param[in]  weights - Interpolation weights in the linear sub-element.
   * \param[out] parCoor - Parametric coordinates to be interpolated.
   */
  void InterpolCoorSubElem(const unsigned short subElem,
                           const su2double      *weights,
                           su2double            *parCoor) override;

private:

  vector<passivedouble> rLineDOFs;    /*!< \brief 1D parametric coordinates of the grid DOFs. */
  vector<passivedouble> rLineSolDOFs; /*!< \brief 1D parametric coordinates of the nodal solution DOFs. */

  ColMajorMatrix<passivedouble> lagBasisLineInt;    /*!< \brief The values of the 1D Lagrangian basis functions
                                                                in the integration points. */
  ColMajorMatrix<passivedouble> derLagBasisLineInt; /*!< \brief The values of the derivatives of the 1D Lagrangian
                                                                basis functions in the integration points. */
  ColMajorMatrix<passivedouble> hesLagBasisLineInt; /*!< \brief The values of the 2nd derivatives, Hessian, of the
                                                                Lagrangian basis function in the integration points. */


  ColMajorMatrix<passivedouble> lagBasisLineSolDOFs;    /*!< \brief The values of the 1D Lagrangian basis functions
                                                                    in the nodal solution DOFs. */
  ColMajorMatrix<passivedouble> derLagBasisLineSolDOFs; /*!< \brief The values of the derivatives of the 1D Lagrangian
                                                                    basis functions in the nodal solution DOFs. */

  TPI2D TensorProductDataVolIntPoints = nullptr; /*!< \brief Function pointer to carry out the tensor product
                                                             to compute the data in the volume integration points. */
  TPI2D TensorProductDataVolSolDOFs   = nullptr; /*!< \brief Function pointer to carry out the tensor product
                                                             to compute the data in the volume nodal solution DOFs. */
};
