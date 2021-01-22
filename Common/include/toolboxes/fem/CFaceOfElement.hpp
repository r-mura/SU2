/*!
 * \file CFaceOfElement.hpp
 * \brief Header file for the class CFaceOfElement.
 *        The implementations are in the <i>CFaceOfElement.cpp</i> file.
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

#include "../../geometry/fem_grid/CVolumeElementFEM_DG.hpp"

#include <vector>
#include <algorithm>
#include <climits>

using namespace std;

/*!
 * \class CFaceOfElement
 * \brief Class used in the partitioning of the FEM grid as well as the building of
          the faces of DG. It stores a face of an element.
 */
class CFaceOfElement {
public:
  unsigned short nCornerPoints;          /*!< \brief Number of corner points of the face. */
  unsigned long  cornerPoints[4];        /*!< \brief Global ID's of ther corner points. */
  unsigned long  elemID0, elemID1;       /*!< \brief Element ID's to the left and right. */
  unsigned short nPolyGrid0, nPolyGrid1; /*!< \brief Polynomial degrees of the grid of the elements
                                                     to the left and right. */
  unsigned short nPolySol0,  nPolySol1;  /*!< \brief Polynomial degrees of the solution of the elements
                                                     to the left and right. */
  unsigned short nDOFsElem0, nDOFsElem1; /*!< \brief Number of DOFs of the elements to the left and right. */
  unsigned short elemType0,  elemType1;  /*!< \brief Type of the elements to the left and right. */
  unsigned short faceID0, faceID1;       /*!< \brief The local face ID in the corresponding elements
                                                     to the left and right of the face. */
  unsigned short periodicIndex;          /*!< \brief Periodic indicator of the face. A value of 0 means no
                                                     periodic face. A value larger than 0 gives the index of
                                                     the periodic boundary + 1. */
  unsigned short periodicIndexDonor;     /*!< \brief Periodic indicator of the donor face. A value of 0 means no
                                                     periodic donor face. A value larger than 0 gives the index of
                                                     the periodic donor boundary + 1. */
  short faceIndicator;                   /*!< \brief The corresponding boundary marker if this face is on a
                                                     boundary. For an internal face the value is -1,
                                                     while an invalidated face has the value -2. */
  bool JacFaceIsConsideredConstant;      /*!< \brief Whether or not the Jacobian of the transformation
                                                     to the standard element is considered constant. */
  bool elem0IsOwner;                     /*!< \brief Whether or not the neighboring element 0 is the owner
                                                     of the face. If false, element 1 is the owner. */

  /*!
   * \brief Default constructor of the class. Set the default values of the members.
   */
  CFaceOfElement();

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFaceOfElement() = default;

  /*!
   * \overload
   * \param[in] VTK_Type - Element type of the face following the VTK convention.
   * \param[in] nPoly    - Polynomial degree of the face.
   * \param[in] Nodes    - Pointer to the array containing the node numbers of the grid DOFs.
   */
  CFaceOfElement(const unsigned short VTK_Type,
                 const unsigned short nPoly,
                 const unsigned long  *Nodes);

  /*!
   * \brief Copy constructor.
   * \param[in] other - Object from which the data must be copied.
   */
  inline CFaceOfElement(const CFaceOfElement &other) { Copy(other); }

  /*!
   * \brief assignment operator.
   * \param[in] other - Object to which the current object must be assigned to.
   */
  inline CFaceOfElement& operator=(const CFaceOfElement &other) { Copy(other); return (*this); }

  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   * \return True if considered less and false otherwise.
   */
  bool operator<(const CFaceOfElement &other) const;

  /*!
   * \brief Equal operator. Needed for removing double entities.
   * \param[in] other - Object to compare against.
   * \return True if considered the same and false otherwise.
   */
  bool operator ==(const CFaceOfElement &other) const;

  /*!
   * \brief Function, which creates a unique numbering for the corner points.
   *        A sort in increasing order is OK for this purpose.
   */
  inline void CreateUniqueNumbering(void) { sort(cornerPoints, cornerPoints+nCornerPoints); }

  /*!
   * \brief Function, which creates a unique numbering for the corner points
   *        while the orientation is taken into account.
   */
  void CreateUniqueNumberingWithOrientation(void);

  /*!
   * \brief Function, which determines the orientation of the element on side 1
   *        compared to the orientation of the face.
   * \param[in] volElem - vector, which contains the local volume elements.
   * \return - Short hand for the orientation of the element on side 1.
   */
  unsigned short DetermineOrientationElemSide1(const vector<CVolumeElementFEM_DG> &volElem) const;

  /*!
   * \brief Function, which makes sure that the orientation of the face
   *        matches the orientation of the corresponding face of the element
   *        on side 0.
   * \param[in] volElem - vector, which contains the local volume elements.
   */
  void MatchOrientationElemSide0(const vector<CVolumeElementFEM_DG> &volElem);

  /*!
   * \brief Function, which swaps the sides of the face if the logic requires this.
   * \param[in] nVolElemOwned - Number of owned volume elements on this MPI rank.
   * \param[in] nVolElemTot   - Total number of volume elements on this MPI rank.
   */
  void SwapSidesIfNeeded(const unsigned long nVolElemOwned,
                         const unsigned long nVolElemTot);

private:
  /*!
   * \brief Copy function, which copies the data of the given object into the current object.
   * \param[in] other - Object to be copied.
   */
  void Copy(const CFaceOfElement &other);
};
