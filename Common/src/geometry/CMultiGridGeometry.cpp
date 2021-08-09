/*!
 * \file CMultiGridGeometry.cpp
 * \brief Implementation of the multigrid geometry class.
 * \author F. Palacios, T. Economon
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/geometry/CMultiGridGeometry.hpp"
#include "../../include/geometry/CMultiGridQueue.hpp"
#include "../../include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

CMultiGridGeometry::CMultiGridGeometry(CGeometry **geometry, CConfig *config_container, unsigned short iMesh) : CGeometry() {

  /*--- CGeometry & CConfig pointers to the fine grid level for clarity. We may
   need access to the other zones in the mesh for zone boundaries. ---*/

  CGeometry *fine_grid = geometry[iMesh-1];
  CConfig *config = config_container;

  /*--- Local variables ---*/

  unsigned long iPoint, Index_CoarseCV, iElem, iVertex, iteration, nVertexS, nVertexR,
                nBufferS_Vector, nBufferR_Vector, iParent, jVertex,Local_nPointCoarse, Local_nPointFine, Global_nPointCoarse, Global_nPointFine,
                *Buffer_Receive_Parent = nullptr, *Buffer_Send_Parent = nullptr, *Buffer_Receive_Children = nullptr, *Buffer_Send_Children = nullptr,
                *Parent_Remote = nullptr,         *Children_Remote = nullptr,    *Parent_Local = nullptr,            *Children_Local = nullptr;
  short marker_seed;
  bool agglomerate_seed = true;
  unsigned short nChildren, iNode, counter, iMarker, jMarker, priority, MarkerS, MarkerR, *nChildren_MPI;
  vector<unsigned long> Suitable_Indirect_Neighbors, Aux_Parent;
  vector<unsigned long>::iterator it;

  unsigned short nMarker_Max = config->GetnMarker_Max();

  unsigned short *copy_marker = new unsigned short [nMarker_Max];

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  nDim = fine_grid->GetnDim(); // Write the number of dimensions of the coarse grid.

  /*--- Create a queue system to deo the agglomeration
   1st) More than two markers ---> Vertices (never agglomerate)
   2nd) Two markers ---> Edges (agglomerate if same BC, never agglomerate if different BC)
   3rd) One marker ---> Surface (always agglomarate)
   4th) No marker ---> Internal Volume (always agglomarate) ---*/

  /*--- Set a marker to indicate indirect agglomeration ---*/

  if (iMesh == MESH_1) {

    for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++)
      fine_grid->nodes->SetAgglomerate_Indirect(iPoint, false);

    for (iElem = 0; iElem < fine_grid->GetnElem(); iElem++) {
      if ((fine_grid->elem[iElem]->GetVTK_Type() == HEXAHEDRON) ||
          (fine_grid->elem[iElem]->GetVTK_Type() == QUADRILATERAL)) {
        for (iNode = 0; iNode < fine_grid->elem[iElem]->GetnNodes(); iNode++) {
          iPoint = fine_grid->elem[iElem]->GetNode(iNode);
          fine_grid->nodes->SetAgglomerate_Indirect(iPoint, true);
        }
      }
    }

  }

  /*--- Create the coarse grid structure using as baseline the fine grid ---*/

  CMultiGridQueue MGQueue_InnerCV(fine_grid->GetnPoint());

  nPointNode = fine_grid->GetnPoint();

  nodes = new CPoint(fine_grid->GetnPoint(), nDim, iMesh, config);

  Index_CoarseCV = 0;

  /*--- The first step is the boundary agglomeration. ---*/

  for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++) {

    for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();

      /*--- If the element has not being previously agglomerated and it belongs
       to the physical domain, then the agglomeration is studied ---*/

      if ((fine_grid->nodes->GetAgglomerate(iPoint) == false) &&
          (fine_grid->nodes->GetDomain(iPoint)) &&
          (GeometricalCheck(iPoint, fine_grid, config))) {

        nChildren = 1;

        /*--- We set an index for the parent control volume ---*/

        fine_grid->nodes->SetParent_CV(iPoint, Index_CoarseCV);

        /*--- We add the seed point (child) to the parent control volume ---*/

        nodes->SetChildren_CV(Index_CoarseCV, 0, iPoint);
        agglomerate_seed = true; counter = 0; marker_seed = iMarker;

        /*--- For a particular point in the fine grid we save all the markers
         that are in that point ---*/

        for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
          if (fine_grid->nodes->GetVertex(iPoint, jMarker) != -1) {
            copy_marker[counter] = jMarker;
            counter++;
          }

        /*--- To aglomerate a vertex it must have only one physical bc!!
         This can be improved. If there is only a marker, it is a good
         candidate for agglomeration ---*/

        if (counter == 1) agglomerate_seed = true;

        /*--- If there are two markers, we will aglomerate if one of the
         marker is SEND_RECEIVE ---*/

        if (counter == 2) {
          if ((config->GetMarker_All_KindBC(copy_marker[0]) == SEND_RECEIVE) ||
              (config->GetMarker_All_KindBC(copy_marker[1]) == SEND_RECEIVE)) agglomerate_seed = true;
          else agglomerate_seed = false;
        }

        /*--- If there are more than 2 markers, the aglomeration will be discarted ---*/

        if (counter > 2) agglomerate_seed = false;

        /*--- If the seed can be agglomerated, we try to agglomerate more points ---*/

        if (agglomerate_seed) {

          /*--- Now we do a sweep over all the nodes that surround the seed point ---*/

          for (auto CVPoint : fine_grid->nodes->GetPoints(iPoint)) {

            /*--- The new point can be agglomerated ---*/

            if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config)) {

              /*--- We set the value of the parent ---*/

              fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

              /*--- We set the value of the child ---*/

              nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
              nChildren++;
            }

          }

          Suitable_Indirect_Neighbors.clear();

          if (fine_grid->nodes->GetAgglomerate_Indirect(iPoint))
            SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);

          /*--- Now we do a sweep over all the indirect nodes that can be added ---*/

          for (auto CVPoint : Suitable_Indirect_Neighbors) {

            /*--- The new point can be agglomerated ---*/

            if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config)) {

              /*--- We set the value of the parent ---*/

              fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

              /*--- We set the indirect agglomeration information ---*/

              if (fine_grid->nodes->GetAgglomerate_Indirect(CVPoint))
                nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);

              /*--- We set the value of the child ---*/

              nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
              nChildren++;
            }
          }


        }

        /*--- Update the number of child of the control volume ---*/

        nodes->SetnChildren_CV(Index_CoarseCV, nChildren);
        Index_CoarseCV++;
      }
    }
  }

  /*--- Agglomerate all the nodes that have more than one physical boundary condition,
   Maybe here we can add the posibility of merging the vertex that have the same number,
   and kind  of markers---*/

  for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++)
    for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      if ((fine_grid->nodes->GetAgglomerate(iPoint) == false) &&
          (fine_grid->nodes->GetDomain(iPoint))) {
        fine_grid->nodes->SetParent_CV(iPoint, Index_CoarseCV);
        nodes->SetChildren_CV(Index_CoarseCV, 0, iPoint);
        nodes->SetnChildren_CV(Index_CoarseCV, 1);
        Index_CoarseCV++;
      }
    }

  /*--- Update the queue with the results from the boundary agglomeration ---*/

  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {

    /*--- The CV has been agglomerated, remove form the list ---*/

    if (fine_grid->nodes->GetAgglomerate(iPoint) == true) {

      MGQueue_InnerCV.RemoveCV(iPoint);

    }

    else {

      /*--- Count the number of agglomerated neighbors, and modify the queue ---*/

      priority = 0;
      for (auto jPoint : fine_grid->nodes->GetPoints(iPoint)) {
        if (fine_grid->nodes->GetAgglomerate(jPoint) == true) priority++;
      }
      MGQueue_InnerCV.MoveCV(iPoint, priority);
    }
  }

  /*--- Agglomerate the domain nodes ---*/

  iteration = 0;
  while (!MGQueue_InnerCV.EmptyQueue() && (iteration < fine_grid->GetnPoint())) {

    iPoint = MGQueue_InnerCV.NextCV();
    iteration ++;

    /*--- If the element has not being previously agglomerated, belongs to the physical domain,
     and satisfies several geometrical criteria then the seed CV is acepted for agglomeration ---*/

    if ((fine_grid->nodes->GetAgglomerate(iPoint) == false) &&
        (fine_grid->nodes->GetDomain(iPoint)) &&
        (GeometricalCheck(iPoint, fine_grid, config))) {

      nChildren = 1;

      /*--- We set an index for the parent control volume ---*/

      fine_grid->nodes->SetParent_CV(iPoint, Index_CoarseCV);

      /*--- We add the seed point (child) to the parent control volume ---*/

      nodes->SetChildren_CV(Index_CoarseCV, 0, iPoint);

      /*--- Update the queue with the seed point (remove the seed and
       increase the priority of the neighbors) ---*/

      MGQueue_InnerCV.Update(iPoint, fine_grid);

      /*--- Now we do a sweep over all the nodes that surround the seed point ---*/

      for (auto CVPoint : fine_grid->nodes->GetPoints(iPoint)) {

        /*--- Determine if the CVPoint can be agglomerated ---*/

        if ((fine_grid->nodes->GetAgglomerate(CVPoint) == false) &&
            (fine_grid->nodes->GetDomain(CVPoint)) &&
            (GeometricalCheck(CVPoint, fine_grid, config))) {

          /*--- We set the value of the parent ---*/

          fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

          /*--- We set the value of the child ---*/

          nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
          nChildren++;

          /*--- Update the queue with the new control volume (remove the CV and
           increase the priority of the neighbors) ---*/

          MGQueue_InnerCV.Update(CVPoint, fine_grid);

        }

      }

      /*--- Subrotuine to identify the indirect neighbors ---*/

      Suitable_Indirect_Neighbors.clear();
      if (fine_grid->nodes->GetAgglomerate_Indirect(iPoint))
        SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);

      /*--- Now we do a sweep over all the indirect nodes that can be added ---*/

      for (auto CVPoint : Suitable_Indirect_Neighbors) {

        /*--- The new point can be agglomerated ---*/

        if ((fine_grid->nodes->GetAgglomerate(CVPoint) == false) &&
            (fine_grid->nodes->GetDomain(CVPoint))) {

          /*--- We set the value of the parent ---*/

          fine_grid->nodes->SetParent_CV(CVPoint, Index_CoarseCV);

          /*--- We set the indirect agglomeration information ---*/

          if (fine_grid->nodes->GetAgglomerate_Indirect(CVPoint))
            nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);

          /*--- We set the value of the child ---*/

          nodes->SetChildren_CV(Index_CoarseCV, nChildren, CVPoint);
          nChildren++;

          /*--- Update the queue with the new control volume (remove the CV and
           increase the priority of the neighbors) ---*/

          MGQueue_InnerCV.Update(CVPoint, fine_grid);

        }
      }

      /*--- Update the number of control of childrens ---*/

      nodes->SetnChildren_CV(Index_CoarseCV, nChildren);
      Index_CoarseCV++;
    }
    else {

      /*--- The seed point can not be agglomerated because of size, domain, streching, etc.
       move the point to the lowest priority ---*/

      MGQueue_InnerCV.MoveCV(iPoint, -1);
    }

  }

  /*--- Add all the elements that have not being agglomerated, in the previous stage ---*/

  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    if ((fine_grid->nodes->GetAgglomerate(iPoint) == false) && (fine_grid->nodes->GetDomain(iPoint))) {

      nChildren = 1;
      fine_grid->nodes->SetParent_CV(iPoint, Index_CoarseCV);
      if (fine_grid->nodes->GetAgglomerate_Indirect(iPoint))
        nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);
      nodes->SetChildren_CV(Index_CoarseCV, 0, iPoint);
      nodes->SetnChildren_CV(Index_CoarseCV, nChildren);
      Index_CoarseCV++;

    }
  }

  nPointDomain = Index_CoarseCV;

  /*--- Check that there are no hanging nodes ---*/

  unsigned long iFinePoint, iCoarsePoint, iCoarsePoint_Complete;
  unsigned short iChildren;

  /*--- Find the point surrounding a point ---*/
  {
    /*--- Temporary, CPoint (nodes) then compresses the information ---*/
    vector<vector<unsigned long> > points(fine_grid->GetnPoint());

    for (iCoarsePoint = 0; iCoarsePoint < nPointDomain; iCoarsePoint ++) {
      for (iChildren = 0; iChildren <  nodes->GetnChildren_CV(iCoarsePoint); iChildren ++) {
        iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        for (auto iFinePoint_Neighbor : fine_grid->nodes->GetPoints(iFinePoint)) {
          iParent = fine_grid->nodes->GetParent_CV(iFinePoint_Neighbor);
          if (iParent != iCoarsePoint) {
            auto End = points[iCoarsePoint].end();
            if (find(points[iCoarsePoint].begin(), End, iParent) == End)
              points[iCoarsePoint].push_back(iParent);
          }
        }
      }
    }
    nodes->SetPoints(points);
  }

  /*--- Detect isolated points and merge them with its correct neighbor ---*/

  for (iCoarsePoint = 0; iCoarsePoint < nPointDomain; iCoarsePoint ++) {

    if (nodes->GetnPoint(iCoarsePoint) == 1) {

      /*--- Find the neighbor of the isolated point. This neighbor is the right control volume ---*/

      iCoarsePoint_Complete = nodes->GetPoint(iCoarsePoint, 0);

      /*--- Add the children to the connected control volume (and modify it parent indexing).
       Identify the child CV from the finest grid and added to the correct control volume.
       Set the parent CV of iFinePoint. Instead of using the original
       (iCoarsePoint) one use the new one (iCoarsePoint_Complete) ---*/

      nChildren = nodes->GetnChildren_CV(iCoarsePoint_Complete);

      for (iChildren = 0; iChildren <  nodes->GetnChildren_CV(iCoarsePoint); iChildren ++) {
        iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        nodes->SetChildren_CV(iCoarsePoint_Complete, nChildren, iFinePoint);
        nChildren++;
        fine_grid->nodes->SetParent_CV(iFinePoint, iCoarsePoint_Complete);
      }

      /*--- Update the number of children control volumes ---*/

      nodes->SetnChildren_CV(iCoarsePoint_Complete, nChildren);
      nodes->SetnChildren_CV(iCoarsePoint, 0);

    }
  }

  //  unsigned long iPointFree = nPointDomain-1;
  //  iCoarsePoint = 0;
  //
  //  do {
  //
  //    if (nodes->GetnChildren_CV(iCoarsePoint) == 0) {
  //
  //      while (nodes->GetnChildren_CV(iPointFree) == 0) {
  //        Index_CoarseCV--;
  //        iPointFree--;
  //      }
  //
  //      nChildren = nodes->GetnChildren_CV(iPointFree);
  //      for (iChildren = 0; iChildren <  nChildren; iChildren ++) {
  //        iFinePoint = nodes->GetChildren_CV(iPointFree, iChildren);
  //        nodes->SetChildren_CV(iCoarsePoint, iChildren, iFinePoint);
  //        fine_grid->nodes->SetParent_CV(iFinePoint, iCoarsePoint);
  //      }
  //      nodes->SetnChildren_CV(iCoarsePoint, nChildren);
  //      nodes->SetnChildren_CV(iPointFree, 0);
  //
  //      Index_CoarseCV--;
  //      iPointFree--;
  //
  //    }
  //
  //    iCoarsePoint++;
  //
  //  } while ((iCoarsePoint-1) < Index_CoarseCV);
  //
  //  nPointDomain = Index_CoarseCV;

  /*--- Reset the point surrounding a point ---*/

  nodes->ResetPoints();

  /*--- Dealing with MPI parallelization, the objective is that the received nodes must be agglomerated
   in the same way as the donor nodes. Send the node agglomeration information of the donor
   (parent and children), Sending only occurs with MPI ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = fine_grid->nVertex[MarkerS];   nVertexR = fine_grid->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;               nBufferR_Vector = nVertexR;

      /*--- Allocate Receive and send buffers  ---*/

      Buffer_Receive_Children = new unsigned long [nBufferR_Vector];
      Buffer_Send_Children = new unsigned long [nBufferS_Vector];

      Buffer_Receive_Parent = new unsigned long [nBufferR_Vector];
      Buffer_Send_Parent = new unsigned long [nBufferS_Vector];

      /*--- Copy the information that should be sended ---*/

      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = fine_grid->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Children[iVertex] = iPoint;
        Buffer_Send_Parent[iVertex] = fine_grid->nodes->GetParent_CV(iPoint);
      }

#ifdef HAVE_MPI
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Children, nBufferS_Vector, MPI_UNSIGNED_LONG, send_to,0,
                   Buffer_Receive_Children, nBufferR_Vector, MPI_UNSIGNED_LONG, receive_from,0, SU2_MPI::GetComm(), &status);
      SU2_MPI::Sendrecv(Buffer_Send_Parent, nBufferS_Vector, MPI_UNSIGNED_LONG, send_to,1,
                   Buffer_Receive_Parent, nBufferR_Vector, MPI_UNSIGNED_LONG, receive_from,1, SU2_MPI::GetComm(), &status);
#else
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Children[iVertex] = Buffer_Send_Children[iVertex];
        Buffer_Receive_Parent[iVertex] = Buffer_Send_Parent[iVertex];
      }
#endif

      /*--- Deallocate send buffer ---*/

      delete [] Buffer_Send_Children;
      delete [] Buffer_Send_Parent;

      /*--- Create a list of the parent nodes without repeated parents ---*/

      Aux_Parent.clear();
      for (iVertex = 0; iVertex < nVertexR; iVertex++)
        Aux_Parent.push_back (Buffer_Receive_Parent[iVertex]);

      sort(Aux_Parent.begin(), Aux_Parent.end());
      it = unique(Aux_Parent.begin(), Aux_Parent.end());
      Aux_Parent.resize(it - Aux_Parent.begin());

      /*--- Allocate some structures ---*/

      Parent_Remote = new unsigned long[nVertexR];
      Children_Remote = new unsigned long[nVertexR];
      Parent_Local = new unsigned long[nVertexR];
      Children_Local = new unsigned long[nVertexR];

      /*--- Create the local vector and remote for the parents and the children ---*/

      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        Parent_Remote[iVertex] = Buffer_Receive_Parent[iVertex];

        /*--- We use the same sorting as in the donor domain ---*/

        for (jVertex = 0; jVertex < Aux_Parent.size(); jVertex++) {
          if (Parent_Remote[iVertex] == Aux_Parent[jVertex]) {
            Parent_Local[iVertex] = jVertex + Index_CoarseCV;
            break;
          }
        }

        Children_Remote[iVertex] = Buffer_Receive_Children[iVertex];
        Children_Local[iVertex] = fine_grid->vertex[MarkerR][iVertex]->GetNode();

      }

      Index_CoarseCV += Aux_Parent.size();

      nChildren_MPI = new unsigned short [Index_CoarseCV];
      for (iParent = 0; iParent < Index_CoarseCV; iParent++)
        nChildren_MPI[iParent] = 0;

      /*--- Create the final structure ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Be careful, it is possible that a node change the agglomeration configuration, the priority
         is always, when receive the information ---*/

        fine_grid->nodes->SetParent_CV(Children_Local[iVertex], Parent_Local[iVertex]);
        nodes->SetChildren_CV(Parent_Local[iVertex], nChildren_MPI[Parent_Local[iVertex]], Children_Local[iVertex]);
        nChildren_MPI[Parent_Local[iVertex]]++;
        nodes->SetnChildren_CV(Parent_Local[iVertex], nChildren_MPI[Parent_Local[iVertex]]);
        nodes->SetDomain(Parent_Local[iVertex], false);

      }

      /*--- Deallocate auxiliar structures ---*/

      delete[] nChildren_MPI;
      delete[] Parent_Remote;
      delete[] Children_Remote;
      delete[] Parent_Local;
      delete[] Children_Local;

      /*--- Deallocate receive buffer ---*/

      delete [] Buffer_Receive_Children;
      delete [] Buffer_Receive_Parent;

    }

  }

  /*--- Update the number of points after the MPI agglomeration ---*/

  nPoint = Index_CoarseCV;

  /*--- Console output with the summary of the agglomeration ---*/

  Local_nPointCoarse = nPoint;
  Local_nPointFine = fine_grid->GetnPoint();

  SU2_MPI::Allreduce(&Local_nPointCoarse, &Global_nPointCoarse, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&Local_nPointFine, &Global_nPointFine, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  su2double Coeff = 1.0, CFL = 0.0, factor = 1.5;

  if (iMesh != MESH_0) {
    if (nDim == 2) Coeff = pow(su2double(Global_nPointFine)/su2double(Global_nPointCoarse), 1./2.);
    if (nDim == 3) Coeff = pow(su2double(Global_nPointFine)/su2double(Global_nPointCoarse), 1./3.);
    CFL = factor*config->GetCFL(iMesh-1)/Coeff;
    config->SetCFL(iMesh, CFL);
  }

  su2double ratio = su2double(Global_nPointFine)/su2double(Global_nPointCoarse);

  if (((nDim == 2) && (ratio < 2.5)) ||
      ((nDim == 3) && (ratio < 2.5))) {
    config->SetMGLevels(iMesh-1);
  }
  else {
    if (rank == MASTER_NODE) {
      PrintingToolbox::CTablePrinter MGTable(&std::cout);
      MGTable.AddColumn("MG Level", 10);
      MGTable.AddColumn("CVs", 10);
      MGTable.AddColumn("Aggl. Rate", 10);
      MGTable.AddColumn("CFL", 10);
      MGTable.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);


      if (iMesh == 1){
        MGTable.PrintHeader();
        MGTable << iMesh - 1 << Global_nPointFine << "1/1.00" << config->GetCFL(iMesh -1);
      }
      stringstream ss;
      ss << "1/" << std::setprecision(3) << ratio;
      MGTable << iMesh << Global_nPointCoarse << ss.str() << CFL;
      if (iMesh == config->GetnMGLevels()){
        MGTable.PrintFooter();
      }
    }
  }

  edgeColorGroupSize = config->GetEdgeColoringGroupSize();

  delete [] copy_marker;

}

bool CMultiGridGeometry::SetBoundAgglomeration(unsigned long CVPoint, short marker_seed, CGeometry *fine_grid, CConfig *config) {

  bool agglomerate_CV = false;
  unsigned short counter, jMarker;

  unsigned short nMarker_Max = config->GetnMarker_Max();

  unsigned short *copy_marker = new unsigned short [nMarker_Max];

  /*--- Basic condition, the element has not being previously agglomerated, it belongs to the domain,
   and has passed some basic geometrical check ---*/

  if ((fine_grid->nodes->GetAgglomerate(CVPoint) == false) &&
      (fine_grid->nodes->GetDomain(CVPoint)) &&
      (GeometricalCheck(CVPoint, fine_grid, config))) {

    /*--- If the element belong to the boundary, we must be careful ---*/

    if (fine_grid->nodes->GetBoundary(CVPoint)) {

      /*--- Identify the markers of the vertex that we want to agglomerate ---*/

      counter = 0;
      for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
        if (fine_grid->nodes->GetVertex(CVPoint, jMarker) != -1) {
          copy_marker[counter] = jMarker;
          counter++;
        }

      /*--- The basic condition is that the aglomerated vertex must have the same physical marker,
       but eventually a send-receive condition ---*/

      /*--- Only one marker in the vertex that is going to be aglomerated ---*/

      if (counter == 1) {

        /*--- We agglomerate if there is only a marker and is the same marker as the seed marker ---*/

        if (copy_marker[0] == marker_seed)
          agglomerate_CV = true;

        /*--- If there is only a marker, but the marker is the SEND_RECEIVE ---*/

        if (config->GetMarker_All_KindBC(copy_marker[0]) == SEND_RECEIVE)
          agglomerate_CV = true;

      }

      /*--- If there are two markers in the vertex that is going to be aglomerated ---*/

      if (counter == 2) {

        /*--- First we verify that the seed is a physical boundary ---*/

        if (config->GetMarker_All_KindBC(marker_seed) != SEND_RECEIVE) {

          /*--- Then we check that one of the marker is equal to the seed marker, and the other is send/receive ---*/

          if (((copy_marker[0] == marker_seed) && (config->GetMarker_All_KindBC(copy_marker[1]) == SEND_RECEIVE)) ||
              ((config->GetMarker_All_KindBC(copy_marker[0]) == SEND_RECEIVE) && (copy_marker[1] == marker_seed)))
            agglomerate_CV = true;
        }

      }

    }

    /*--- If the element belong to the domain, it is allways aglomerated ---*/

    else { agglomerate_CV = true; }

  }

  delete [] copy_marker;

  return agglomerate_CV;

}


bool CMultiGridGeometry::GeometricalCheck(unsigned long iPoint, CGeometry *fine_grid, CConfig *config) {

  su2double max_dimension = 1.2;

  /*--- Evaluate the total size of the element ---*/

  bool Volume = true;
  su2double ratio = pow(fine_grid->nodes->GetVolume(iPoint), 1.0/su2double(nDim))*max_dimension;
  su2double limit = pow(config->GetDomainVolume(), 1.0/su2double(nDim));
  if ( ratio > limit ) Volume = false;

  /*--- Evaluate the stretching of the element ---*/

  bool Stretching = true;

  /* unsigned short iNode, iDim;
   unsigned long jPoint;
   su2double *Coord_i = fine_grid->nodes->GetCoord(iPoint);
   su2double max_dist = 0.0 ; su2double min_dist = 1E20;
   for (iNode = 0; iNode < fine_grid->nodes->GetnPoint(iPoint); iNode ++) {
   jPoint = fine_grid->nodes->GetPoint(iPoint, iNode);
   su2double *Coord_j = fine_grid->nodes->GetCoord(jPoint);
   su2double distance = 0.0;
   for (iDim = 0; iDim < nDim; iDim++)
   distance += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
   distance = sqrt(distance);
   max_dist = max(distance, max_dist);
   min_dist = min(distance, min_dist);
   }
   if ( max_dist/min_dist > 100.0 ) Stretching = false;*/

  return (Stretching && Volume);

}

void CMultiGridGeometry::SetSuitableNeighbors(vector<unsigned long> *Suitable_Indirect_Neighbors, unsigned long iPoint,
                                              unsigned long Index_CoarseCV, CGeometry *fine_grid) {

  unsigned short iNeighbor, jNeighbor;
  bool SecondNeighborSeed, ThirdNeighborSeed;
  vector<unsigned long>::iterator it;

  /*--- Create a list with the first neighbors, including the seed ---*/

  vector<unsigned long> First_Neighbor_Points;
  First_Neighbor_Points.push_back(iPoint);
  for (auto jPoint : fine_grid->nodes->GetPoints(iPoint))
    First_Neighbor_Points.push_back(jPoint);

  /*--- Create a list with the second neighbors, without first, and seed neighbors ---*/

  vector<unsigned long> Second_Neighbor_Points, Second_Origin_Points, Suitable_Second_Neighbors;

  for (auto jPoint : fine_grid->nodes->GetPoints(iPoint)) {

    for (auto kPoint : fine_grid->nodes->GetPoints(jPoint)) {

      /*--- Check that the second neighbor do not belong to the first neighbor or the seed ---*/

      SecondNeighborSeed = true;
      for (iNeighbor = 0; iNeighbor < First_Neighbor_Points.size(); iNeighbor ++)
        if (kPoint == First_Neighbor_Points[iNeighbor]) {
          SecondNeighborSeed = false; break;
        }

      if (SecondNeighborSeed) {
        Second_Neighbor_Points.push_back(kPoint);
        Second_Origin_Points.push_back(jPoint);
      }

    }
  }

  /*---  Identify those second neighbors that are repeated (candidate to be added) ---*/

  for (iNeighbor = 0; iNeighbor < Second_Neighbor_Points.size(); iNeighbor ++)

    for (jNeighbor = 0; jNeighbor < Second_Neighbor_Points.size(); jNeighbor ++)

    /*--- Repeated second neighbor with different origin ---*/

      if ((Second_Neighbor_Points[iNeighbor] == Second_Neighbor_Points[jNeighbor]) &&
          (Second_Origin_Points[iNeighbor] != Second_Origin_Points[jNeighbor]) &&
          (iNeighbor < jNeighbor)) {

        Suitable_Indirect_Neighbors->push_back(Second_Neighbor_Points[iNeighbor]);

        /*--- Create alist with the suitable second neighbor, that we will use
         to compute the third neighbors --*/

        Suitable_Second_Neighbors.push_back(Second_Neighbor_Points[iNeighbor]);

      }


  /*--- Remove repeated from the suitable second neighbors ---*/

  sort(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  it = unique(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  Suitable_Second_Neighbors.resize(it - Suitable_Second_Neighbors.begin());

  /*--- Remove repeated from first neighbors ---*/

  sort(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
  it = unique(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
  First_Neighbor_Points.resize(it - First_Neighbor_Points.begin());

  /*--- Create a list with the third neighbors, without first, second, and seed neighbors ---*/

  vector<unsigned long> Third_Neighbor_Points, Third_Origin_Points;

  for (auto kPoint : Suitable_Second_Neighbors) {
    for (auto lPoint : fine_grid->nodes->GetPoints(kPoint)) {

      /*--- Check that the third neighbor do not belong to the first neighbors or the seed ---*/

      ThirdNeighborSeed = true;

      for (iNeighbor = 0; iNeighbor < First_Neighbor_Points.size(); iNeighbor ++)
        if (lPoint == First_Neighbor_Points[iNeighbor]) {
          ThirdNeighborSeed = false;
          break;
        }

      /*--- Check that the third neighbor do not belong to the second neighbors ---*/

      for (iNeighbor = 0; iNeighbor < Suitable_Second_Neighbors.size(); iNeighbor ++)
        if (lPoint == Suitable_Second_Neighbors[iNeighbor]) {
          ThirdNeighborSeed = false;
          break;
        }

      if (ThirdNeighborSeed) {
        Third_Neighbor_Points.push_back(lPoint);
        Third_Origin_Points.push_back(kPoint);
      }

    }
  }

  /*---  Identify those third neighbors that are repeated (candidate to be added) ---*/

  for (iNeighbor = 0; iNeighbor < Third_Neighbor_Points.size(); iNeighbor ++)
    for (jNeighbor = 0; jNeighbor < Third_Neighbor_Points.size(); jNeighbor ++)

    /*--- Repeated second neighbor with different origin ---*/

      if ((Third_Neighbor_Points[iNeighbor] == Third_Neighbor_Points[jNeighbor]) &&
          (Third_Origin_Points[iNeighbor] != Third_Origin_Points[jNeighbor]) &&
          (iNeighbor < jNeighbor)) {

        Suitable_Indirect_Neighbors->push_back(Third_Neighbor_Points[iNeighbor]);

      }

  /*--- Remove repeated from Suitable Indirect Neighbors List ---*/

  sort(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
  it = unique(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
  Suitable_Indirect_Neighbors->resize(it - Suitable_Indirect_Neighbors->begin());

}

void CMultiGridGeometry::SetPoint_Connectivity(CGeometry *fine_grid) {

  unsigned long iFinePoint, iParent, iCoarsePoint;
  unsigned short iChildren;

  /*--- Set the point surrounding a point ---*/

  vector<vector<unsigned long> > points(nPoint);

  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    for (iChildren = 0; iChildren <  nodes->GetnChildren_CV(iCoarsePoint); iChildren ++) {
      iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
      for (auto iFinePoint_Neighbor : fine_grid->nodes->GetPoints(iFinePoint)) {
        iParent = fine_grid->nodes->GetParent_CV(iFinePoint_Neighbor);
        if (iParent != iCoarsePoint) {
          auto End = points[iCoarsePoint].end();
          if (find(points[iCoarsePoint].begin(), End, iParent) == End)
            points[iCoarsePoint].push_back(iParent);
        }
      }
    }
  }
  nodes->SetPoints(points);

  /*--- Set the number of neighbors variable, this is
   important for JST and multigrid in parallel ---*/

  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    nodes->SetnNeighbor(iCoarsePoint, nodes->GetnPoint(iCoarsePoint));

}

void CMultiGridGeometry::SetVertex(CGeometry *fine_grid, CConfig *config) {
  unsigned long  iVertex, iFinePoint, iCoarsePoint;
  unsigned short iMarker, iMarker_Tag, iChildren;

  nMarker = fine_grid->GetnMarker();
  unsigned short nMarker_Max = config->GetnMarker_Max();

  /*--- If any children node belong to the boundary then the entire control
   volume will belong to the boundary ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    for (iChildren = 0; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren ++) {
      iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
      if (fine_grid->nodes->GetBoundary(iFinePoint)) {
        nodes->SetBoundary(iCoarsePoint, nMarker);
        break;
      }
    }

  vertex = new CVertex**[nMarker];
  nVertex = new unsigned long [nMarker];

  Tag_to_Marker = new string [nMarker_Max];
  for (iMarker_Tag = 0; iMarker_Tag < nMarker_Max; iMarker_Tag++)
    Tag_to_Marker[iMarker_Tag] = fine_grid->GetMarker_Tag(iMarker_Tag);

  /*--- Compute the number of vertices to do the dimensionalization ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;


  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    if (nodes->GetBoundary(iCoarsePoint)) {
      for (iChildren = 0; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren ++) {
        iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        for (iMarker = 0; iMarker < nMarker; iMarker ++) {
          if ((fine_grid->nodes->GetVertex(iFinePoint, iMarker) != -1) && (nodes->GetVertex(iCoarsePoint, iMarker) == -1)) {
            iVertex = nVertex[iMarker];
            nodes->SetVertex(iCoarsePoint, iVertex, iMarker);
            nVertex[iMarker]++;
          }
        }
      }
    }
  }

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex* [fine_grid->GetnVertex(iMarker)+1];
    nVertex[iMarker] = 0;
  }

  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    if (nodes->GetBoundary(iCoarsePoint))
      for (iMarker = 0; iMarker < nMarker; iMarker ++)
        nodes->SetVertex(iCoarsePoint, -1, iMarker);

  for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;

  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    if (nodes->GetBoundary(iCoarsePoint)) {
      for (iChildren = 0; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren ++) {
        iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker ++) {
          if ((fine_grid->nodes->GetVertex(iFinePoint, iMarker) != -1) && (nodes->GetVertex(iCoarsePoint, iMarker) == -1)) {
            iVertex = nVertex[iMarker];
            vertex[iMarker][iVertex] = new CVertex(iCoarsePoint, nDim);
            nodes->SetVertex(iCoarsePoint, iVertex, iMarker);

            /*--- Set the transformation to apply ---*/
            unsigned long ChildVertex = fine_grid->nodes->GetVertex(iFinePoint, iMarker);
            unsigned short RotationKind = fine_grid->vertex[iMarker][ChildVertex]->GetRotation_Type();
            vertex[iMarker][iVertex]->SetRotation_Type(RotationKind);
            nVertex[iMarker]++;
          }
        }
      }
    }
  }
}

void CMultiGridGeometry::MatchNearField(CConfig *config) {

  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  int iProcessor = size;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (nodes->GetDomain(iPoint)) {
          vertex[iMarker][iVertex]->SetDonorPoint(iPoint, nodes->GetGlobalIndex(iPoint), iVertex, iMarker, iProcessor);
        }
      }
    }
  }

}

void CMultiGridGeometry::MatchActuator_Disk(CConfig *config) {

  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  int iProcessor = size;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
        (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (nodes->GetDomain(iPoint)) {
          vertex[iMarker][iVertex]->SetDonorPoint(iPoint, nodes->GetGlobalIndex(iPoint), iVertex, iMarker, iProcessor);
        }
      }
    }
  }

}

void CMultiGridGeometry::MatchPeriodic(CConfig *config, unsigned short val_periodic) {

  unsigned short iMarker, iPeriodic, nPeriodic;
  unsigned long iVertex, iPoint;
  int iProcessor = rank;

  /*--- Evaluate the number of periodic boundary conditions ---*/

  nPeriodic = config->GetnMarker_Periodic();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (nodes->GetDomain(iPoint)) {
            vertex[iMarker][iVertex]->SetDonorPoint(iPoint, nodes->GetGlobalIndex(iPoint), iVertex, iMarker, iProcessor);
          }
        }
      }
    }
  }

}

void CMultiGridGeometry::SetControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {

  SU2_OMP_MASTER {

  unsigned long iFinePoint, iCoarsePoint, iEdge, iParent;
  long FineEdge, CoarseEdge;
  unsigned short iChildren;
  bool change_face_orientation;
  su2double Coarse_Volume, Area;

  /*--- Compute the area of the coarse volume ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    nodes->SetVolume(iCoarsePoint, 0.0);
    Coarse_Volume = 0.0;
    for (iChildren = 0; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren ++) {
      iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
      Coarse_Volume += fine_grid->nodes->GetVolume(iFinePoint);
    }
    nodes->SetVolume(iCoarsePoint, Coarse_Volume);
  }

  /*--- Update or not the values of faces at the edge ---*/
  if (action != ALLOCATE) {
    edges->SetZeroValues();
  }

  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    for (iChildren = 0; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren ++) {
      iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);

      for (auto iFinePoint_Neighbor : fine_grid->nodes->GetPoints(iFinePoint)) {
        iParent = fine_grid->nodes->GetParent_CV(iFinePoint_Neighbor);
        if ((iParent != iCoarsePoint) && (iParent < iCoarsePoint)) {

          FineEdge = fine_grid->FindEdge(iFinePoint, iFinePoint_Neighbor);

          change_face_orientation = false;
          if (iFinePoint < iFinePoint_Neighbor) change_face_orientation = true;

          CoarseEdge = FindEdge(iParent, iCoarsePoint);

          const auto Normal = fine_grid->edges->GetNormal(FineEdge);

          if (change_face_orientation) {
            edges->SubNormal(CoarseEdge,Normal);
          }
          else {
            edges->AddNormal(CoarseEdge,Normal);
          }
        }
      }
    }

  /*--- Check if there is a normal with null area ---*/

  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    const auto NormalFace = edges->GetNormal(iEdge);
    Area = GeometryToolbox::Norm(nDim, NormalFace);
    if (Area == 0.0) {
      su2double DefaultNormal[3] = {EPS*EPS};
      edges->SetNormal(iEdge, DefaultNormal);
    }
  }

  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER
}

void CMultiGridGeometry::SetBoundControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {

  SU2_OMP_MASTER {

  unsigned long iCoarsePoint, iFinePoint, FineVertex, iVertex;
  unsigned short iMarker, iChildren, iDim;
  su2double *Normal, Area, *NormalFace = nullptr;

  Normal = new su2double [nDim];

  if (action != ALLOCATE) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        vertex[iMarker][iVertex]->SetZeroValues();
  }

  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      iCoarsePoint = vertex[iMarker][iVertex]->GetNode();
      for (iChildren = 0; iChildren < nodes->GetnChildren_CV(iCoarsePoint); iChildren ++) {
        iFinePoint = nodes->GetChildren_CV(iCoarsePoint, iChildren);
        if (fine_grid->nodes->GetVertex(iFinePoint, iMarker)!=-1) {
          FineVertex = fine_grid->nodes->GetVertex(iFinePoint, iMarker);
          fine_grid->vertex[iMarker][FineVertex]->GetNormal(Normal);
          vertex[iMarker][iVertex]->AddNormal(Normal);
        }
      }
    }

  delete[] Normal;

  /*--- Check if there is a normal with null area ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      NormalFace = vertex[iMarker][iVertex]->GetNormal();
      Area = GeometryToolbox::Norm(nDim, NormalFace);
      if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }

  }
  END_SU2_OMP_MASTER
  SU2_OMP_BARRIER
}

void CMultiGridGeometry::SetCoord(CGeometry *geometry) {

  SU2_OMP_FOR_STAT(roundUpDiv(nPoint, omp_get_max_threads()))
  for (auto Point_Coarse = 0ul; Point_Coarse < nPoint; Point_Coarse++) {
    auto Area_Parent = nodes->GetVolume(Point_Coarse);
    su2double Coordinates[3] = {0.0};
    for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
      auto Point_Fine = nodes->GetChildren_CV(Point_Coarse, iChildren);
      auto Area_Children = geometry->nodes->GetVolume(Point_Fine);
      auto Coordinates_Fine = geometry->nodes->GetCoord(Point_Fine);
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Coordinates[iDim] += Coordinates_Fine[iDim]*Area_Children/Area_Parent;
    }
    nodes->SetCoord(Point_Coarse, Coordinates);
  }
  END_SU2_OMP_FOR
}

void CMultiGridGeometry::SetMultiGridWallHeatFlux(CGeometry *geometry, unsigned short val_marker){

  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short iChildren;
  long Vertex_Fine;
  su2double Area_Parent, Area_Children;
  su2double WallHeatFlux_Fine, WallHeatFlux_Coarse;
  bool isVertex;
  int numberVertexChildren;

  for(iVertex=0; iVertex < nVertex[val_marker]; iVertex++){
    Point_Coarse = vertex[val_marker][iVertex]->GetNode();
    if (nodes->GetDomain(Point_Coarse)){
      Area_Parent = 0.0;
      WallHeatFlux_Coarse = 0.0;
      numberVertexChildren = 0;
      /*--- Compute area parent by taking into account only volumes that are on the marker ---*/
      for(iChildren=0; iChildren < nodes->GetnChildren_CV(Point_Coarse); iChildren++){
        Point_Fine = nodes->GetChildren_CV(Point_Coarse, iChildren);
        isVertex = (nodes->GetDomain(Point_Fine) && geometry->nodes->GetVertex(Point_Fine, val_marker) != -1);
        if (isVertex){
          numberVertexChildren += 1;
          Area_Parent += geometry->nodes->GetVolume(Point_Fine);
        }
      }

      /*--- Loop again and propagate values to the coarser level ---*/
      for(iChildren=0; iChildren < nodes->GetnChildren_CV(Point_Coarse); iChildren++){
        Point_Fine = nodes->GetChildren_CV(Point_Coarse, iChildren);
        Vertex_Fine = geometry->nodes->GetVertex(Point_Fine, val_marker);
        isVertex = (nodes->GetDomain(Point_Fine) && Vertex_Fine != -1);
        if(isVertex){
          Area_Children = geometry->nodes->GetVolume(Point_Fine);
          //Get the customized BC values on fine level and compute the values at coarse level
          WallHeatFlux_Fine = geometry->GetCustomBoundaryHeatFlux(val_marker, Vertex_Fine);
          WallHeatFlux_Coarse += WallHeatFlux_Fine*Area_Children/Area_Parent;
        }

      }
      //Set the customized BC values at coarse level
      CustomBoundaryHeatFlux[val_marker][iVertex] = WallHeatFlux_Coarse;
    }
  }

}

void CMultiGridGeometry::SetMultiGridWallTemperature(CGeometry *geometry, unsigned short val_marker){

  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short iChildren;
  long Vertex_Fine;
  su2double Area_Parent, Area_Children;
  su2double WallTemperature_Fine, WallTemperature_Coarse;
  bool isVertex;
  int numberVertexChildren;

  for(iVertex=0; iVertex < nVertex[val_marker]; iVertex++){
    Point_Coarse = vertex[val_marker][iVertex]->GetNode();
    if (nodes->GetDomain(Point_Coarse)){
      Area_Parent = 0.0;
      WallTemperature_Coarse = 0.0;
      numberVertexChildren = 0;
      /*--- Compute area parent by taking into account only volumes that are on the marker ---*/
      for(iChildren=0; iChildren < nodes->GetnChildren_CV(Point_Coarse); iChildren++){
        Point_Fine = nodes->GetChildren_CV(Point_Coarse, iChildren);
        isVertex = (nodes->GetDomain(Point_Fine) && geometry->nodes->GetVertex(Point_Fine, val_marker) != -1);
        if (isVertex){
          numberVertexChildren += 1;
          Area_Parent += geometry->nodes->GetVolume(Point_Fine);
        }
      }

      /*--- Loop again and propagate values to the coarser level ---*/
      for(iChildren=0; iChildren < nodes->GetnChildren_CV(Point_Coarse); iChildren++){
        Point_Fine = nodes->GetChildren_CV(Point_Coarse, iChildren);
        Vertex_Fine = geometry->nodes->GetVertex(Point_Fine, val_marker);
        isVertex = (nodes->GetDomain(Point_Fine) && Vertex_Fine != -1);
        if(isVertex){
          Area_Children = geometry->nodes->GetVolume(Point_Fine);
          //Get the customized BC values on fine level and compute the values at coarse level
          WallTemperature_Fine = geometry->GetCustomBoundaryTemperature(val_marker, Vertex_Fine);
          WallTemperature_Coarse += WallTemperature_Fine*Area_Children/Area_Parent;
        }

      }
      //Set the customized BC values at coarse level
      CustomBoundaryTemperature[val_marker][iVertex] = WallTemperature_Coarse;
    }
  }

}

void CMultiGridGeometry::SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config) {

  /*--- Loop over all coarse mesh points. ---*/
  SU2_OMP_FOR_STAT(roundUpDiv(nPoint,omp_get_max_threads()))
  for (unsigned long Point_Coarse = 0; Point_Coarse < nPoint; Point_Coarse++) {
    su2double Area_Parent = nodes->GetVolume(Point_Coarse);

    /*--- Initialize coarse grid velocity to zero. ---*/
    su2double Grid_Vel[3] = {0.0, 0.0, 0.0};

    /*--- Loop over all of the children for this coarse CV and compute
     a grid velocity based on the values in the child CVs (fine mesh). ---*/
    for (unsigned short iChild = 0; iChild < nodes->GetnChildren_CV(Point_Coarse); iChild++) {
      unsigned long Point_Fine       = nodes->GetChildren_CV(Point_Coarse, iChild);
      su2double Area_Child           = fine_mesh->nodes->GetVolume(Point_Fine);
      const su2double* Grid_Vel_Fine = fine_mesh->nodes->GetGridVel(Point_Fine);
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Grid_Vel[iDim] += Grid_Vel_Fine[iDim]*Area_Child/Area_Parent;
    }

    /*--- Set the grid velocity for this coarse node. ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      nodes->SetGridVel(Point_Coarse, iDim, Grid_Vel[iDim]);
  }
  END_SU2_OMP_FOR
}


void CMultiGridGeometry::FindNormal_Neighbor(CConfig *config) {

  unsigned short iMarker, iDim;
  unsigned long iPoint, iVertex;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY ) {

      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {

        iPoint = vertex[iMarker][iVertex]->GetNode();

        /*--- If the node belong to the domain ---*/
        if (nodes->GetDomain(iPoint)) {

          /*--- Compute closest normal neighbor ---*/
          su2double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord;
          unsigned long Point_Normal = 0;
          su2double *Normal = vertex[iMarker][iVertex]->GetNormal();
          cos_max = -1.0;
          for (auto jPoint : nodes->GetPoints(iPoint)) {
            scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              diff_coord = nodes->GetCoord(jPoint, iDim)-nodes->GetCoord(iPoint, iDim);
              scalar_prod += diff_coord*Normal[iDim];
              norm_vect += diff_coord*diff_coord;
              norm_Normal += Normal[iDim]*Normal[iDim];
            }
            norm_vect = sqrt(norm_vect);
            norm_Normal = sqrt(norm_Normal);
            cos_alpha = scalar_prod/(norm_vect*norm_Normal);

            /*--- Get maximum cosine (not minimum because normals are oriented inwards) ---*/
            if (cos_alpha >= cos_max) {
              Point_Normal = jPoint;
              cos_max = cos_alpha;
            }
          }
          vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
        }
      }
    }
  }
}

void CMultiGridGeometry::ComputeNSpan(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate) {
  unsigned short iMarker, jMarker, iMarkerTP, iSpan, jSpan;
  unsigned long iPoint, iVertex;
  long jVertex;
  int nSpan, nSpan_loc;
  su2double *coord, *valueSpan, min, max, radius, delta;
  short PeriodicBoundary;
  unsigned short SpanWise_Kind = config->GetKind_SpanWise();

  unsigned short iSize;
  int nSpan_max;
  su2double MyMin, MyMax;

  nSpan = 0;
  nSpan_loc = 0;
  if (nDim == 2){
    nSpanWiseSections[marker_flag-1] = 1;
    //TODO (turbo) make it more genral
    if(marker_flag == OUTFLOW) config->SetnSpanWiseSections(1);

    /*---Initilize the vector of span-wise values that will be ordered ---*/
    SpanWiseValue[marker_flag -1] = new su2double[1];
    for (iSpan = 0; iSpan < 1; iSpan++){
      SpanWiseValue[marker_flag -1][iSpan] = 0;
    }
  }
  else{
    if(SpanWise_Kind == AUTOMATIC){
      /*--- loop to find inflow of outflow marker---*/
      for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){

          if (config->GetMarker_All_Turbomachinery(iMarker) != iMarkerTP) continue;
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) != marker_flag) continue;

          /*--- loop to find the vertex that ar both of inflow or outflow marker and on the periodic
           * in order to caount the number of Span ---*/
          for (jMarker = 0; jMarker < nMarker; jMarker++){
            if (config->GetMarker_All_KindBC(jMarker) != PERIODIC_BOUNDARY) continue;

            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {

              iPoint = vertex[iMarker][iVertex]->GetNode();
              if (!nodes->GetDomain(iPoint)) continue;

              PeriodicBoundary = config->GetMarker_All_PerBound(jMarker);
              jVertex = nodes->GetVertex(iPoint, jMarker);

              if ((jVertex != -1) && (PeriodicBoundary == (val_iZone + 1))){
                nSpan++;
              }
            }
          }
        }
      }

      /*--- storing the local number of span---*/
      nSpan_loc = nSpan;
      SU2_MPI::Allreduce(&nSpan_loc, &nSpan, 1, MPI_INT, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&nSpan_loc, &nSpan_max, 1, MPI_INT, MPI_MAX, SU2_MPI::GetComm());

      /*--- initialize the vector that will contain the disordered values span-wise ---*/
      nSpanWiseSections[marker_flag -1] = nSpan;
      valueSpan = new su2double[nSpan];

      for (iSpan = 0; iSpan < nSpan; iSpan ++ ){
        valueSpan[iSpan] = -1001.0;
      }

      /*--- store the local span-wise value for each processor ---*/
      nSpan_loc = 0;
      for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){

          if (config->GetMarker_All_Turbomachinery(iMarker) != iMarkerTP) continue;
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) != marker_flag) continue;

          for (jMarker = 0; jMarker < nMarker; jMarker++){
            if (config->GetMarker_All_KindBC(jMarker) != PERIODIC_BOUNDARY) continue;

            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {

              iPoint = vertex[iMarker][iVertex]->GetNode();
              if (!nodes->GetDomain(iPoint)) continue;

              PeriodicBoundary = config->GetMarker_All_PerBound(jMarker);
              jVertex = nodes->GetVertex(iPoint, jMarker);

              if ((jVertex != -1) && (PeriodicBoundary == (val_iZone + 1))){
                coord = nodes->GetCoord(iPoint);
                radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                switch (config->GetKind_TurboMachinery(val_iZone)){
                case CENTRIFUGAL:
                  valueSpan[nSpan_loc] = coord[2];
                  break;
                case CENTRIPETAL:
                  valueSpan[nSpan_loc] = coord[2];
                  break;
                case AXIAL:
                  valueSpan[nSpan_loc] = radius;
                  break;
                case CENTRIPETAL_AXIAL:
                  if (marker_flag == OUTFLOW){
                    valueSpan[nSpan_loc] = radius;
                  }
                  else{
                    valueSpan[nSpan_loc] = coord[2];
                  }
                  break;
                case AXIAL_CENTRIFUGAL:
                  if (marker_flag == INFLOW){
                    valueSpan[nSpan_loc] = radius;
                  }
                  else{
                    valueSpan[nSpan_loc] = coord[2];
                  }
                  break;
                }
                nSpan_loc++;
              }
            }
          }
        }
      }

      /*--- Gather the span-wise values on all the processor ---*/

      vector<su2double> MyTotValueSpan(nSpan_max*size, -1001.0);
      vector<su2double> MyValueSpan(nSpan_max, -1001.0);
      vector<int> My_nSpan_loc(size);

      for(iSpan = 0; iSpan<nSpan_loc; iSpan++){
        MyValueSpan[iSpan] = valueSpan[iSpan];
      }

      SU2_MPI::Allgather(MyValueSpan.data(), nSpan_max , MPI_DOUBLE, MyTotValueSpan.data(), nSpan_max, MPI_DOUBLE, SU2_MPI::GetComm());
      SU2_MPI::Allgather(&nSpan_loc, 1 , MPI_INT, My_nSpan_loc.data(), 1, MPI_INT, SU2_MPI::GetComm());

      jSpan = 0;
      for (iSize = 0; iSize< size; iSize++){
        for(iSpan = 0; iSpan < My_nSpan_loc[iSize]; iSpan++){
          valueSpan[jSpan] = MyTotValueSpan[iSize*nSpan_max + iSpan];
          jSpan++;
        }
      }
      if (jSpan != nSpan) SU2_MPI::Error("Panic!",CURRENT_FUNCTION);

      /*--- Terrible stuff to do but so is this entire bloody function goodness me... ---*/
      SpanWiseValue[marker_flag -1] = valueSpan;

      sort(SpanWiseValue[marker_flag -1], SpanWiseValue[marker_flag -1]+nSpan);

      /*--- Find the minimum value among the span-wise values  ---*/
      min = SpanWiseValue[marker_flag -1][0];
    }
    /*--- Compute equispaced Span-wise sections using number of section specified by the User---*/
    else{
      /*--- Initialize number of span---*/
      nSpanWiseSections[marker_flag-1] = config->Get_nSpanWiseSections_User();
      SpanWiseValue[marker_flag -1] = new su2double[config->Get_nSpanWiseSections_User()];
      for (iSpan = 0; iSpan < config->Get_nSpanWiseSections_User(); iSpan++){
        SpanWiseValue[marker_flag -1][iSpan] = 0;
      }
      /*--- Compute maximum and minimum value span-wise---*/
      min = 1E+07;
      max =-1E+07;
      for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){

          if (config->GetMarker_All_Turbomachinery(iMarker) != iMarkerTP) continue;
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) != marker_flag) continue;

          for (jMarker = 0; jMarker < nMarker; jMarker++){
            if (config->GetMarker_All_KindBC(jMarker) != PERIODIC_BOUNDARY) continue;

            for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {

              iPoint = vertex[iMarker][iVertex]->GetNode();
              if (!nodes->GetDomain(iPoint)) continue;

              PeriodicBoundary = config->GetMarker_All_PerBound(jMarker);
              jVertex = nodes->GetVertex(iPoint, jMarker);

              if ((jVertex != -1) && (PeriodicBoundary == (val_iZone + 1))){

                coord = nodes->GetCoord(iPoint);
                radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                switch (config->GetKind_TurboMachinery(val_iZone)){
                case CENTRIFUGAL:
                case CENTRIPETAL:
                  if (coord[2] < min) min = coord[2];
                  if (coord[2] > max) max = coord[2];
                  break;
                case AXIAL:
                  if (radius < min) min = radius;
                  if (radius > max) max = radius;
                  break;
                case CENTRIPETAL_AXIAL:
                  if (marker_flag == OUTFLOW){
                    if (radius < min) min = radius;
                    if (radius > max) max = radius;
                  }
                  else{
                    if (coord[2] < min) min = coord[2];
                    if (coord[2] > max) max = coord[2];
                  }
                  break;
                case AXIAL_CENTRIFUGAL:
                  if (marker_flag == INFLOW){
                    if (radius < min) min = radius;
                    if (radius > max) max = radius;
                  }
                  else{
                    if (coord[2] < min) min = coord[2];
                    if (coord[2] > max) max = coord[2];
                  }
                  break;
                }
              }
            }
          }
        }
      }
      /*--- compute global minimum and maximum value on span-wise ---*/
      MyMin= min;  min = 0;
      MyMax= max;  max = 0;
      SU2_MPI::Allreduce(&MyMin, &min, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&MyMax, &max, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

      /*--- compute height value for each spanwise section---*/
      delta = (max - min)/(nSpanWiseSections[marker_flag-1] -1);
      for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
        SpanWiseValue[marker_flag - 1][iSpan]= min + delta*iSpan;
      }
    }

    if(marker_flag == OUTFLOW){
      if(nSpanWiseSections[INFLOW -1] != nSpanWiseSections[OUTFLOW - 1]){
        char buf[100];
        SPRINTF(buf, "nSpan inflow %u, nSpan outflow %u", nSpanWiseSections[INFLOW-1], nSpanWiseSections[OUTFLOW-1]);
        SU2_MPI::Error(string(" At the moment only turbomachinery with the same amount of span-wise section can be simulated\n") + buf, CURRENT_FUNCTION);
      }
      else{
        config->SetnSpanWiseSections(nSpanWiseSections[OUTFLOW -1]);
      }
    }

  }

}

void CMultiGridGeometry::SetTurboVertex(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate) {
  unsigned long  iPoint, **ordered, **disordered, **oldVertex3D, iInternalVertex;
  unsigned long nVert, nVertMax;
  unsigned short iMarker, iMarkerTP, iSpan, jSpan, iDim;
  su2double min, minInt, max, *coord, dist, Normal2, *TurboNormal, *NormalArea, target = 0.0, **area, ***unitnormal, Area = 0.0;
  bool **checkAssign;
  min    =  10.0E+06;
  minInt =  10.0E+06;
  max    = -10.0E+06;

  su2double radius;
  long iVertex, iSpanVertex, jSpanVertex, kSpanVertex = 0;
  int *nTotVertex_gb, *nVertexSpanHalo;
  su2double **x_loc, **y_loc, **z_loc, **angCoord_loc, **deltaAngCoord_loc, **angPitch, **deltaAngPitch, *minIntAngPitch,
  *minAngPitch, *maxAngPitch;
  int       **rank_loc;

#ifdef HAVE_MPI
  unsigned short iSize, kSize = 0, jSize;
  su2double MyMin,MyIntMin, MyMax;
  su2double *x_gb = NULL, *y_gb = NULL, *z_gb = NULL, *angCoord_gb = NULL, *deltaAngCoord_gb = NULL;
  bool *checkAssign_gb =NULL;
  unsigned long My_nVert;

#endif
  string multizone_filename;

  x_loc              = new su2double*[nSpanWiseSections[marker_flag-1]];
  y_loc              = new su2double*[nSpanWiseSections[marker_flag-1]];
  z_loc              = new su2double*[nSpanWiseSections[marker_flag-1]];
  angCoord_loc       = new su2double*[nSpanWiseSections[marker_flag-1]];
  deltaAngCoord_loc  = new su2double*[nSpanWiseSections[marker_flag-1]];
  angPitch           = new su2double*[nSpanWiseSections[marker_flag-1]];
  deltaAngPitch      = new su2double*[nSpanWiseSections[marker_flag-1]];
  rank_loc           = new int*[nSpanWiseSections[marker_flag-1]];
  minAngPitch        = new su2double[nSpanWiseSections[marker_flag-1]];
  minIntAngPitch     = new su2double[nSpanWiseSections[marker_flag-1]];
  maxAngPitch        = new su2double[nSpanWiseSections[marker_flag-1]];

  nTotVertex_gb      = new int[nSpanWiseSections[marker_flag-1]];
  nVertexSpanHalo    = new int[nSpanWiseSections[marker_flag-1]];
  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    nTotVertex_gb[iSpan]   = -1;
    nVertexSpanHalo[iSpan] = 0;
    minAngPitch[iSpan]     = 10.0E+06;
    minIntAngPitch[iSpan]  = 10.0E+06;
    maxAngPitch[iSpan]     = -10.0E+06;
  }

  /*--- Initialize auxiliary pointers ---*/
  TurboNormal        = new su2double[3];
  NormalArea         = new su2double[3];
  ordered            = new unsigned long* [nSpanWiseSections[marker_flag-1]];
  disordered         = new unsigned long* [nSpanWiseSections[marker_flag-1]];
  oldVertex3D        = new unsigned long* [nSpanWiseSections[marker_flag-1]];
  area               = new su2double* [nSpanWiseSections[marker_flag-1]];
  unitnormal         = new su2double** [nSpanWiseSections[marker_flag-1]];
  checkAssign        = new bool* [nSpanWiseSections[marker_flag-1]];

  /*--- Initialize the new Vertex structure. The if statement ensures that these vectors are initialized
   * only once even if the routine is called more than once.---*/

  if (allocate){
      for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
          if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
            if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
            nSpanSectionsByMarker[iMarker]       = nSpanWiseSections[marker_flag-1];
            nVertexSpan[iMarker]                 = new long[nSpanWiseSections[marker_flag-1]];
            turbovertex[iMarker]                 = new CTurboVertex** [nSpanWiseSections[marker_flag-1]];
            nTotVertexSpan[iMarker]              = new unsigned long [nSpanWiseSections[marker_flag-1] +1];
            MaxAngularCoord[iMarker]             = new su2double [nSpanWiseSections[marker_flag-1]];
            MinAngularCoord[iMarker]             = new su2double [nSpanWiseSections[marker_flag-1]];
            MinRelAngularCoord[iMarker]          = new su2double [nSpanWiseSections[marker_flag-1]];
            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
              nVertexSpan[iMarker][iSpan]        = 0;
              turbovertex[iMarker][iSpan]        = nullptr;
              MinAngularCoord[iMarker][iSpan]    = 10.0E+06;
              MaxAngularCoord[iMarker][iSpan]    = -10.0E+06;
              MinRelAngularCoord[iMarker][iSpan] = 10.0E+06;
            }
            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1] +1; iSpan++){
              nTotVertexSpan[iMarker][iSpan]     = 0;
            }
          }
        }
      }
    }
  }

  //this works only for turbomachinery rotating around the Z-Axes.
  // the reordering algorithm pitch-wise assumes that X-coordinate of each boundary vertex is positive so that reordering can be based on the Y-coordinate.
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){

            /*--- compute the amount of vertexes for each span-wise section to initialize the CTurboVertex pointers and auxiliary pointers  ---*/
            for (iVertex = 0; (unsigned long)iVertex  < nVertex[iMarker]; iVertex++) {
              iPoint = vertex[iMarker][iVertex]->GetNode();
              if (nDim == 3){
                dist = 10E+06;
                jSpan = std::numeric_limits<unsigned short>::max();
                coord = nodes->GetCoord(iPoint);

                switch (config->GetKind_TurboMachinery(val_iZone)){
                case CENTRIFUGAL: case CENTRIPETAL:
                  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                    if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                      dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                      jSpan=iSpan;
                    }
                  }
                  break;
                case AXIAL:
                  radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                    if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                      dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                      jSpan=iSpan;
                    }
                  }
                  break;
                case CENTRIPETAL_AXIAL:
                  if (marker_flag == OUTFLOW){
                    radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  else{
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  break;

                case AXIAL_CENTRIFUGAL:
                  if (marker_flag == INFLOW){
                    radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  else{
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  break;
                }
              }

              /*--- 2D problem do not need span-wise separation---*/
              else{
                jSpan = 0;
              }

              if(nodes->GetDomain(iPoint)){
                nVertexSpan[iMarker][jSpan]++;
              }
              nVertexSpanHalo[jSpan]++;
            }

            /*--- initialize the CTurboVertex pointers and auxiliary pointers  ---*/
            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
              if (allocate){
                turbovertex[iMarker][iSpan] = new CTurboVertex* [nVertexSpan[iMarker][iSpan]];
                for (iVertex = 0; iVertex < nVertexSpan[iMarker][iSpan]; iVertex++){
                  turbovertex[iMarker][iSpan][iVertex] = nullptr;
                }
              }
              ordered[iSpan]                           = new unsigned long [nVertexSpanHalo[iSpan]];
              disordered[iSpan]                        = new unsigned long [nVertexSpanHalo[iSpan]];
              oldVertex3D[iSpan]                       = new unsigned long [nVertexSpanHalo[iSpan]];
              checkAssign[iSpan]                       = new bool [nVertexSpanHalo[iSpan]];
              area[iSpan]                              = new su2double [nVertexSpanHalo[iSpan]];
              unitnormal[iSpan]                        = new su2double* [nVertexSpanHalo[iSpan]];
              for (iVertex = 0; iVertex < nVertexSpanHalo[iSpan]; iVertex++){
                unitnormal[iSpan][iVertex]             = new su2double [nDim];
              }
              angPitch[iSpan]                          = new su2double [nVertexSpanHalo[iSpan]];
              deltaAngPitch[iSpan]                     = new su2double [nVertexSpanHalo[iSpan]];
              nVertexSpanHalo[iSpan]                   = 0;
            }

            /*--- store the vertexes in a ordered manner in span-wise directions but not yet ordered pitch-wise ---*/
            for (iVertex = 0; (unsigned long)iVertex < nVertex[iMarker]; iVertex++) {
              iPoint = vertex[iMarker][iVertex]->GetNode();
              if(nDim == 3){
                dist  = 10E+06;
                jSpan = std::numeric_limits<unsigned short>::max();

                coord = nodes->GetCoord(iPoint);
                switch (config->GetKind_TurboMachinery(val_iZone)){
                case CENTRIFUGAL: case CENTRIPETAL:
                  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                    if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                      dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                      jSpan=iSpan;
                    }
                  }
                  break;
                case AXIAL:
                  radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                    if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                      dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                      jSpan=iSpan;
                    }
                  }
                  break;
                case CENTRIPETAL_AXIAL:
                  if(marker_flag == OUTFLOW){
                    radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }else{
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  break;

                case AXIAL_CENTRIFUGAL:
                  if(marker_flag == INFLOW){
                    radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }else{
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  break;
                }
              }
              /*--- 2D problem do not need span-wise separation---*/
              else{
                jSpan = 0;
              }
              /*--- compute the face area associated with the vertex ---*/
              vertex[iMarker][iVertex]->GetNormal(NormalArea);
              for (iDim = 0; iDim < nDim; iDim++) NormalArea[iDim] = -NormalArea[iDim];
              Area = GeometryToolbox::Norm(nDim, NormalArea);

              for (iDim = 0; iDim < nDim; iDim++) NormalArea[iDim] /= Area;
              /*--- store all the all the info into the auxiliary containers ---*/
              disordered[jSpan][nVertexSpanHalo[jSpan]]  = iPoint;
              oldVertex3D[jSpan][nVertexSpanHalo[jSpan]] = iVertex;
              area[jSpan][nVertexSpanHalo[jSpan]]        = Area;
              for (iDim = 0; iDim < nDim; iDim++){
                unitnormal[jSpan][nVertexSpanHalo[jSpan]][iDim] = NormalArea[iDim];
              }
              checkAssign[jSpan][nVertexSpanHalo[jSpan]] = false;
              nVertexSpanHalo[jSpan]++;
            }

            /*--- using the auxiliary container reordered the vertexes pitch-wise direction at each span ---*/
            // the reordering algorithm can be based on the Y-coordinate.
            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){

              /*--- find the local minimum and maximum pitch-wise for each processor---*/
              min    = 10E+06;
              minInt = 10E+06;
              max    = -10E+06;
              for(iSpanVertex = 0; iSpanVertex < nVertexSpanHalo[iSpan]; iSpanVertex++){
                iPoint = disordered[iSpan][iSpanVertex];
                coord = nodes->GetCoord(iPoint);
                /*--- find nodes at minimum pitch among all nodes---*/
                if (coord[1]<min){
                  min = coord[1];
                  if (nDim == 2 && config->GetKind_TurboMachinery(val_iZone) == AXIAL){
                    MinAngularCoord[iMarker][iSpan] = coord[1];
                  }
                  else{
                  MinAngularCoord[iMarker][iSpan] = atan(coord[1]/coord[0]);
                  }
                  minAngPitch[iSpan]= MinAngularCoord[iMarker][iSpan];
                  kSpanVertex =iSpanVertex;
                }

                /*--- find nodes at minimum pitch among the internal nodes---*/
                if (coord[1]<minInt){
                  if(nodes->GetDomain(iPoint)){
                    minInt = coord[1];
                    if (nDim == 2 && config->GetKind_TurboMachinery(val_iZone) == AXIAL){
                      minIntAngPitch[iSpan] = coord[1];
                    }
                    else{
                      minIntAngPitch[iSpan] = atan(coord[1]/coord[0]);
                    }
                  }
                }

                /*--- find nodes at maximum pitch among the internal nodes---*/
                if (coord[1]>max){
                  if(nodes->GetDomain(iPoint)){
                    max =coord[1];
                    if (nDim == 2 && config->GetKind_TurboMachinery(val_iZone) == AXIAL){
                      MaxAngularCoord[iMarker][iSpan] = coord[1];
                    }
                    else{
                      MaxAngularCoord[iMarker][iSpan] = atan(coord[1]/coord[0]);
                    }
                    maxAngPitch[iSpan]= MaxAngularCoord[iMarker][iSpan];
                  }
                }
              }

              iInternalVertex = 0;

              /*--- reordering the vertex pitch-wise, store the ordered vertexes span-wise and pitch-wise---*/
              for(iSpanVertex = 0; iSpanVertex<nVertexSpanHalo[iSpan]; iSpanVertex++){
                dist = 10E+06;
                ordered[iSpan][iSpanVertex] = disordered[iSpan][kSpanVertex];
                checkAssign[iSpan][kSpanVertex] = true;
                coord = nodes->GetCoord(ordered[iSpan][iSpanVertex]);
                target = coord[1];
                if (nDim == 2 && config->GetKind_TurboMachinery(val_iZone) == AXIAL){
                   angPitch[iSpan][iSpanVertex]=coord[1];
                }
                else{
                  angPitch[iSpan][iSpanVertex]=atan(coord[1]/coord[0]);
                }
                if(iSpanVertex == 0){
                  deltaAngPitch[iSpan][iSpanVertex]=0.0;
                }
                else{
                  deltaAngPitch[iSpan][iSpanVertex]= angPitch[iSpan][iSpanVertex] - angPitch[iSpan][iSpanVertex - 1];
                }
                /*---create turbovertex structure only for the internal nodes---*/
                if(nodes->GetDomain(ordered[iSpan][iSpanVertex])){
                  if (allocate){
                    turbovertex[iMarker][iSpan][iInternalVertex] = new CTurboVertex(ordered[iSpan][iSpanVertex], nDim);
                  }
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetArea(area[iSpan][kSpanVertex]);
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetNormal(unitnormal[iSpan][kSpanVertex]);
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetOldVertex(oldVertex3D[iSpan][kSpanVertex]);
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetAngularCoord(angPitch[iSpan][iSpanVertex]);
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetDeltaAngularCoord(deltaAngPitch[iSpan][iSpanVertex]);
                  switch (config->GetKind_TurboMachinery(val_iZone)){
                  case CENTRIFUGAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if (marker_flag == INFLOW){
                      TurboNormal[0] = -coord[0]/sqrt(Normal2);
                      TurboNormal[1] = -coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }else{
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }
                    break;
                  case CENTRIPETAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if (marker_flag == OUTFLOW){
                      TurboNormal[0] = -coord[0]/sqrt(Normal2);
                      TurboNormal[1] = -coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }else{
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }
                    break;
                  case AXIAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if(nDim == 3){
                      if (marker_flag == INFLOW){
                        TurboNormal[0] = coord[0]/sqrt(Normal2);
                        TurboNormal[1] = coord[1]/sqrt(Normal2);
                        TurboNormal[2] = 0.0;
                      }else{
                        TurboNormal[0] = coord[0]/sqrt(Normal2);
                        TurboNormal[1] = coord[1]/sqrt(Normal2);
                        TurboNormal[2] = 0.0;
                      }
                    }
                    else{
                      if (marker_flag == INFLOW){
                        TurboNormal[0] = -1.0;
                        TurboNormal[1] = 0.0;
                        TurboNormal[2] = 0.0;
                      }else{
                        TurboNormal[0] = 1.0;
                        TurboNormal[1] = 0.0;
                        TurboNormal[2] = 0.0;
                      }
                    }

                    break;
                  case CENTRIPETAL_AXIAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if (marker_flag == INFLOW){
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }else{
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }
                    break;

                  case AXIAL_CENTRIFUGAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if (marker_flag == INFLOW){
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }else{
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }
                    break;

                  }
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetTurboNormal(TurboNormal);
                  iInternalVertex++;
                }


                for(jSpanVertex = 0; jSpanVertex<nVertexSpanHalo[iSpan]; jSpanVertex++){
                  coord = nodes->GetCoord(disordered[iSpan][jSpanVertex]);
                  if(dist >= (coord[1] - target) && !checkAssign[iSpan][jSpanVertex] && (coord[1] - target) >= 0.0){
                    dist= coord[1] - target;
                    kSpanVertex =jSpanVertex;
                  }
                }
              }
            }

            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){

              delete [] ordered[iSpan];
              delete [] disordered[iSpan];
              delete [] oldVertex3D[iSpan];
              delete [] checkAssign[iSpan];
              delete [] area[iSpan];
              delete [] angPitch[iSpan];
              delete [] deltaAngPitch[iSpan];

              for(iVertex=0; iVertex < nVertexSpanHalo[iSpan]; iVertex++){
                delete [] unitnormal[iSpan][iVertex];
              }
              delete [] unitnormal[iSpan];
            }
          }
        }
      }
    }

  /*--- to be set for all the processor to initialize an appropriate number of frequency for the NR BC ---*/
  nVertMax = 0;

  /*--- compute global max and min pitch per span ---*/
  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    nVert    = 0;

#ifdef HAVE_MPI
    MyMin     = minAngPitch[iSpan];      minAngPitch[iSpan]    = 10.0E+6;
    MyIntMin  = minIntAngPitch[iSpan];   minIntAngPitch[iSpan] = 10.0E+6;
    MyMax     = maxAngPitch[iSpan];      maxAngPitch[iSpan]    = -10.0E+6;

    SU2_MPI::Allreduce(&MyMin, &minAngPitch[iSpan], 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&MyIntMin, &minIntAngPitch[iSpan], 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&MyMax, &maxAngPitch[iSpan], 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
#endif


    /*--- compute the relative angular pitch with respect to the minimum value ---*/

    for (iMarker = 0; iMarker < nMarker; iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
            nVert = nVertexSpan[iMarker][iSpan];
            MinAngularCoord[iMarker][iSpan]    = minAngPitch[iSpan];
            MaxAngularCoord[iMarker][iSpan]    = maxAngPitch[iSpan];
            MinRelAngularCoord[iMarker][iSpan] = minIntAngPitch[iSpan] - minAngPitch[iSpan];
            for(iSpanVertex = 0; iSpanVertex< nVertexSpan[iMarker][iSpan]; iSpanVertex++){
             turbovertex[iMarker][iSpan][iSpanVertex]->SetRelAngularCoord(MinAngularCoord[iMarker][iSpan]);
            }
          }
        }
      }
    }


#ifdef HAVE_MPI
    My_nVert = nVert;nVert = 0;
    SU2_MPI::Allreduce(&My_nVert, &nVert, 1, MPI_INT, MPI_SUM, SU2_MPI::GetComm());
#endif

    /*--- to be set for all the processor to initialize an appropriate number of frequency for the NR BC ---*/
    if(nVert > nVertMax){
      SetnVertexSpanMax(marker_flag,nVert);
    }
    /*--- for all the processor should be known the amount of total turbovertex per span  ---*/
    nTotVertex_gb[iSpan]= (int)nVert;

    for (iMarker = 0; iMarker < nMarker; iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
            nTotVertexSpan[iMarker][iSpan]= nVert;
            nTotVertexSpan[iMarker][nSpanWiseSections[marker_flag-1]]+= nVert;
          }
        }
      }
    }
  }


  /*--- Printing Tec file to check the global ordering of the turbovertex pitch-wise ---*/
  /*--- Send all the info to the MASTERNODE ---*/

  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    x_loc[iSpan]             = new su2double[nTotVertex_gb[iSpan]];
    y_loc[iSpan]             = new su2double[nTotVertex_gb[iSpan]];
    z_loc[iSpan]             = new su2double[nTotVertex_gb[iSpan]];
    angCoord_loc[iSpan]      = new su2double[nTotVertex_gb[iSpan]];
    deltaAngCoord_loc[iSpan] = new su2double[nTotVertex_gb[iSpan]];
    rank_loc[iSpan]          = new int[nTotVertex_gb[iSpan]];
    for(iSpanVertex = 0; iSpanVertex<nTotVertex_gb[iSpan]; iSpanVertex++){
      x_loc[iSpan][iSpanVertex]             = -1.0;
      y_loc[iSpan][iSpanVertex]             = -1.0;
      z_loc[iSpan][iSpanVertex]             = -1.0;
      angCoord_loc[iSpan][iSpanVertex]      = -1.0;
      deltaAngCoord_loc[iSpan][iSpanVertex] = -1.0;
      rank_loc[iSpan][iSpanVertex]          = -1;
    }
  }

  for (iMarker = 0; iMarker < nMarker; iMarker++){
    for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
      if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
          for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
            for(iSpanVertex = 0; iSpanVertex<nVertexSpan[iMarker][iSpan]; iSpanVertex++){
              iPoint = turbovertex[iMarker][iSpan][iSpanVertex]->GetNode();
              coord  = nodes->GetCoord(iPoint);
              x_loc[iSpan][iSpanVertex]   = coord[0];
              y_loc[iSpan][iSpanVertex]   = coord[1];
              if (nDim == 3){
                z_loc[iSpan][iSpanVertex] = coord[2];
              }
              else{
                z_loc[iSpan][iSpanVertex] = 0.0;
              }
              angCoord_loc[iSpan][iSpanVertex]      = turbovertex[iMarker][iSpan][iSpanVertex]->GetRelAngularCoord();
              deltaAngCoord_loc[iSpan][iSpanVertex] = turbovertex[iMarker][iSpan][iSpanVertex]->GetDeltaAngularCoord();
            }
          }
        }
      }
    }
  }

#ifdef HAVE_MPI

  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    if (rank == MASTER_NODE){
      x_gb                = new su2double[nTotVertex_gb[iSpan]*size];
      y_gb                = new su2double[nTotVertex_gb[iSpan]*size];
      z_gb                = new su2double[nTotVertex_gb[iSpan]*size];
      angCoord_gb         = new su2double[nTotVertex_gb[iSpan]*size];
      deltaAngCoord_gb    = new su2double[nTotVertex_gb[iSpan]*size];
      checkAssign_gb      = new bool[nTotVertex_gb[iSpan]*size];

     for(iSize= 0; iSize < size; iSize++){
       for(iSpanVertex = 0; iSpanVertex < nTotVertex_gb[iSpan]; iSpanVertex++){
         checkAssign_gb[iSize*nTotVertex_gb[iSpan] + iSpanVertex] = false;
       }
     }
    }
    SU2_MPI::Gather(y_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, y_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(x_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, x_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(z_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, z_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(angCoord_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, angCoord_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(deltaAngCoord_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, deltaAngCoord_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());

    if (rank == MASTER_NODE){
      for(iSpanVertex = 0; iSpanVertex<nTotVertex_gb[iSpan]; iSpanVertex++){
        x_loc[iSpan][iSpanVertex]             = -1.0;
        y_loc[iSpan][iSpanVertex]             = -1.0;
        z_loc[iSpan][iSpanVertex]             = -1.0;
        angCoord_loc[iSpan][iSpanVertex]      = -1.0;
        deltaAngCoord_loc[iSpan][iSpanVertex] = -1.0;
      }



      min = 10.0E+06;
      for(iSize= 0; iSize < size; iSize++){
        if (angCoord_gb[iSize*nTotVertex_gb[iSpan]] < min && angCoord_gb[iSize*nTotVertex_gb[iSpan]] >= 0.0){
          kSize = iSize;
          min = angCoord_gb[iSize*nTotVertex_gb[iSpan]];
        }
      }

      kSpanVertex = 0;
      for(iSpanVertex = 0; iSpanVertex < nTotVertex_gb[iSpan]; iSpanVertex++){
        x_loc[iSpan][iSpanVertex]              = x_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        y_loc[iSpan][iSpanVertex]              = y_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        z_loc[iSpan][iSpanVertex]              = z_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        angCoord_loc[iSpan][iSpanVertex]       = angCoord_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        deltaAngCoord_loc[iSpan][iSpanVertex]  = deltaAngCoord_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        rank_loc[iSpan][iSpanVertex]           = kSize;
        target = angCoord_loc[iSpan][iSpanVertex];
        checkAssign_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex] = true;
        min = 10.0E+06;
        for(jSize= 0; jSize < size; jSize++){
          for(jSpanVertex = 0; jSpanVertex < nTotVertex_gb[iSpan]; jSpanVertex++){
            if ((angCoord_gb[jSize*nTotVertex_gb[iSpan] + jSpanVertex] < min) &&
                (angCoord_gb[jSize*nTotVertex_gb[iSpan] + jSpanVertex] >= target) &&
                !checkAssign_gb[jSize*nTotVertex_gb[iSpan] + jSpanVertex]) {
              kSize = jSize;
              kSpanVertex = jSpanVertex;
              min = angCoord_gb[jSize*nTotVertex_gb[iSpan] + jSpanVertex];
            }
          }
        }
      }

      delete [] x_gb; delete [] y_gb; delete [] z_gb;  delete [] angCoord_gb; delete [] deltaAngCoord_gb; delete[] checkAssign_gb;

    }
  }

#endif

  if (rank == MASTER_NODE){
    if (marker_flag == INFLOW && val_iZone ==0){
      std::string sPath = "TURBOMACHINERY";
      int nError = 0;
#if defined(_WIN32)
#ifdef __MINGW32__
      nError = mkdir(sPath.c_str());  // MINGW on Windows
#else
      nError = _mkdir(sPath.c_str()); // can be used on Windows
#endif
#else
      mode_t nMode = 0733; // UNIX style permissions
      //nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows // LOOK AT THIS
#endif
      if (nError != 0) {
        cout << "TURBOMACHINERY folder creation failed." <<endl;
      }
    }
    if (marker_flag == INFLOW){
      multizone_filename = "TURBOMACHINERY/spanwise_division_inflow.dat";
    }
    else{
      multizone_filename = "TURBOMACHINERY/spanwise_division_outflow.dat";
    }
    char buffer[50];

    if (GetnZone() > 1){
      unsigned short lastindex = multizone_filename.find_last_of(".");
      multizone_filename = multizone_filename.substr(0, lastindex);
      SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(val_iZone));
      multizone_filename.append(string(buffer));
    }

    // File to print the vector x_loc, y_loc, z_loc, globIdx_loc to check vertex ordering
    ofstream myfile;
    myfile.open (multizone_filename.data(), ios::out | ios::trunc);
    myfile.setf(ios::uppercase | ios::scientific);
    myfile.precision(8);

    myfile << "TITLE = \"Global index visualization file\"" << endl;
    myfile << "VARIABLES =" << endl;
    myfile.width(10); myfile << "\"iSpan\"";
    myfile.width(20); myfile << "\"x_coord\"" ;
    myfile.width(20); myfile << "\"y_coord\"" ;
    myfile.width(20); myfile << "\"z_coord\"" ;
    myfile.width(20); myfile << "\"radius\"" ;
    myfile.width(20); myfile << "\"Relative Angular Coord \"" ;
    myfile.width(20); myfile << "\"Delta Angular Coord \"" ;
    myfile.width(20); myfile << "\"processor\"" <<endl;
    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
      for(iSpanVertex = 0; iSpanVertex < nTotVertex_gb[iSpan]; iSpanVertex++){
        radius = sqrt(x_loc[iSpan][iSpanVertex]*x_loc[iSpan][iSpanVertex] + y_loc[iSpan][iSpanVertex]*y_loc[iSpan][iSpanVertex]);
        myfile.width(10); myfile << iSpan;
        myfile.width(20); myfile << x_loc[iSpan][iSpanVertex];
        myfile.width(20); myfile << y_loc[iSpan][iSpanVertex];
        myfile.width(20); myfile << z_loc[iSpan][iSpanVertex];
        myfile.width(20); myfile << radius;
        if (nDim ==2 && config->GetKind_TurboMachinery(val_iZone)){
          myfile.width(20); myfile << angCoord_loc[iSpan][iSpanVertex];
          myfile.width(20); myfile << deltaAngCoord_loc[iSpan][iSpanVertex];
        }
        else{
          myfile.width(20); myfile << angCoord_loc[iSpan][iSpanVertex]*180.0/PI_NUMBER;
          myfile.width(20); myfile << deltaAngCoord_loc[iSpan][iSpanVertex]*180.0/PI_NUMBER;
        }
        myfile.width(20); myfile << rank_loc[iSpan][iSpanVertex]<<endl;
      }
      myfile << endl;
    }
  }


  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    delete [] x_loc[iSpan];
    delete [] y_loc[iSpan];
    delete [] z_loc[iSpan];
    delete [] angCoord_loc[iSpan];
    delete [] deltaAngCoord_loc[iSpan];
    delete [] rank_loc[iSpan];

  }


  delete [] area;
  delete [] ordered;
  delete [] disordered;
  delete [] oldVertex3D;
  delete [] checkAssign;
  delete [] TurboNormal;
  delete [] unitnormal;
  delete [] NormalArea;
  delete [] x_loc;
  delete [] y_loc;
  delete [] z_loc;
  delete [] angCoord_loc;
  delete [] nTotVertex_gb;
  delete [] nVertexSpanHalo;
  delete [] angPitch;
  delete [] deltaAngPitch;
  delete [] deltaAngCoord_loc;
  delete [] rank_loc;
  delete [] minAngPitch;
  delete [] maxAngPitch;
  delete [] minIntAngPitch;

}

void CMultiGridGeometry::SetAvgTurboValue(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate) {

  unsigned short iMarker, iMarkerTP, iSpan, iDim;
  unsigned long iPoint;
  su2double *TurboNormal,*coord, *Normal, turboNormal2, Normal2, *gridVel, TotalArea, TotalRadius, radius;
  su2double *TotalTurboNormal,*TotalNormal, *TotalGridVel, Area;
  long iVertex;
  /*-- Variables declaration and allocation ---*/
  TotalTurboNormal = new su2double[nDim];
  TotalNormal      = new su2double[nDim];
  TurboNormal      = new su2double[nDim];
  TotalGridVel     = new su2double[nDim];
  Normal           = new su2double[nDim];

  bool grid_movement        = config->GetGrid_Movement();
#ifdef HAVE_MPI
  su2double MyTotalArea, MyTotalRadius, *MyTotalTurboNormal= NULL, *MyTotalNormal= NULL, *MyTotalGridVel= NULL;
#endif

  /*--- Intialization of the vector for the interested boundary ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
      if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
          if(allocate){
            AverageTurboNormal[iMarker]               = new su2double *[nSpanWiseSections[marker_flag-1] + 1];
            AverageNormal[iMarker]                    = new su2double *[nSpanWiseSections[marker_flag-1] + 1];
            AverageGridVel[iMarker]                   = new su2double *[nSpanWiseSections[marker_flag-1] + 1];
            AverageTangGridVel[iMarker]               = new su2double [nSpanWiseSections[marker_flag-1] + 1];
            SpanArea[iMarker]                         = new su2double [nSpanWiseSections[marker_flag-1] + 1];
            TurboRadius[iMarker]                      = new su2double [nSpanWiseSections[marker_flag-1] + 1];
            for (iSpan= 0; iSpan < nSpanWiseSections[marker_flag-1] + 1; iSpan++){
              AverageTurboNormal[iMarker][iSpan]      = new su2double [nDim];
              AverageNormal[iMarker][iSpan]           = new su2double [nDim];
              AverageGridVel[iMarker][iSpan]          = new su2double [nDim];
            }
          }
          for (iSpan= 0; iSpan < nSpanWiseSections[marker_flag-1] + 1; iSpan++){
            AverageTangGridVel[iMarker][iSpan]          = 0.0;
            SpanArea[iMarker][iSpan]                    = 0.0;
            TurboRadius[iMarker][iSpan]                 = 0.0;
            for(iDim=0; iDim < nDim; iDim++){
              AverageTurboNormal[iMarker][iSpan][iDim]  = 0.0;
              AverageNormal[iMarker][iSpan][iDim]       = 0.0;
              AverageGridVel[iMarker][iSpan][iDim]      = 0.0;
            }
          }
        }
      }
    }
  }



  /*--- start computing the average quantities span wise --- */
  for (iSpan= 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){

    /*--- Forces initialization for contenitors to zero ---*/
    for (iDim=0; iDim<nDim; iDim++) {
      TotalTurboNormal[iDim]  =0.0;
      TotalNormal[iDim]         =0.0;
      TotalGridVel[iDim]        =0.0;
    }
    TotalArea = 0.0;
    TotalRadius = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
            for(iVertex = 0; iVertex < nVertexSpan[iMarker][iSpan]; iVertex++){
              iPoint = turbovertex[iMarker][iSpan][iVertex]->GetNode();
              turbovertex[iMarker][iSpan][iVertex]->GetTurboNormal(TurboNormal);
              turbovertex[iMarker][iSpan][iVertex]->GetNormal(Normal);
              coord  = nodes->GetCoord(iPoint);

              if (nDim == 3){
                radius = sqrt(coord[0]*coord[0] + coord[1]*coord[1]);
              }
              else{
                radius = 0.0;
              }
              Area = turbovertex[iMarker][iSpan][iVertex]->GetArea();
              TotalArea   += Area;
              TotalRadius += radius;
              for (iDim = 0; iDim < nDim; iDim++) {
                TotalTurboNormal[iDim]  +=TurboNormal[iDim];
                TotalNormal[iDim]       +=Normal[iDim];
              }
              if (grid_movement){
                gridVel = nodes->GetGridVel(iPoint);
                for (iDim = 0; iDim < nDim; iDim++) TotalGridVel[iDim] +=gridVel[iDim];
              }
            }
          }
        }
      }
    }

#ifdef HAVE_MPI

    MyTotalArea            = TotalArea;                 TotalArea            = 0;
    MyTotalRadius          = TotalRadius;               TotalRadius          = 0;
    SU2_MPI::Allreduce(&MyTotalArea, &TotalArea, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&MyTotalRadius, &TotalRadius, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

    MyTotalTurboNormal     = new su2double[nDim];
    MyTotalNormal          = new su2double[nDim];
    MyTotalGridVel         = new su2double[nDim];

    for (iDim = 0; iDim < nDim; iDim++) {
      MyTotalTurboNormal[iDim]  = TotalTurboNormal[iDim];
      TotalTurboNormal[iDim]    = 0.0;
      MyTotalNormal[iDim]       = TotalNormal[iDim];
      TotalNormal[iDim]         = 0.0;
      MyTotalGridVel[iDim]      = TotalGridVel[iDim];
      TotalGridVel[iDim]        = 0.0;
    }

    SU2_MPI::Allreduce(MyTotalTurboNormal, TotalTurboNormal, nDim, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(MyTotalNormal, TotalNormal, nDim, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(MyTotalGridVel, TotalGridVel, nDim, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

    delete [] MyTotalTurboNormal;delete [] MyTotalNormal; delete [] MyTotalGridVel;

#endif

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){


            SpanArea[iMarker][iSpan]           = TotalArea;
            TurboRadius[iMarker][iSpan]        = TotalRadius/nTotVertexSpan[iMarker][iSpan];

            turboNormal2    = 0.0;
            Normal2         = 0.0;
            for (iDim = 0; iDim < nDim; iDim++){
              turboNormal2 += TotalTurboNormal[iDim]*TotalTurboNormal[iDim];
              Normal2      += TotalNormal[iDim]*TotalNormal[iDim];
            }
            for (iDim = 0; iDim < nDim; iDim++){
              AverageTurboNormal[iMarker][iSpan][iDim] = TotalTurboNormal[iDim]/sqrt(turboNormal2);
              AverageNormal[iMarker][iSpan][iDim]      = TotalNormal[iDim]/sqrt(Normal2);
            }
            if (grid_movement){
              for (iDim = 0; iDim < nDim; iDim++){
                AverageGridVel[iMarker][iSpan][iDim]   =TotalGridVel[iDim]/nTotVertexSpan[iMarker][iSpan];
              }
              switch (config->GetKind_TurboMachinery(val_iZone)){
              case CENTRIFUGAL:case CENTRIPETAL:
                if (marker_flag == INFLOW ){
                  AverageTangGridVel[iMarker][iSpan]= -(AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0]);
                }
                else{
                  AverageTangGridVel[iMarker][iSpan]= AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }
                break;
              case AXIAL:
                if (marker_flag == INFLOW && nDim == 2){
                  AverageTangGridVel[iMarker][iSpan]= -AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1] + AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }
                else{
                  AverageTangGridVel[iMarker][iSpan]= AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }
                  break;
              case CENTRIPETAL_AXIAL:
                if (marker_flag == OUTFLOW){
                  AverageTangGridVel[iMarker][iSpan]= (AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0]);
                }
                else{
                  AverageTangGridVel[iMarker][iSpan]= -(AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0]);
                }
                break;
              case AXIAL_CENTRIFUGAL:
                if (marker_flag == INFLOW)
                {
                  AverageTangGridVel[iMarker][iSpan]= AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }else
                {
                  AverageTangGridVel[iMarker][iSpan]= AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }
                break;

              default:
                  SU2_MPI::Error("Tang grid velocity NOT IMPLEMENTED YET for this configuration", CURRENT_FUNCTION);
                break;
              }
            }

            /*--- Compute the 1D average values ---*/
            AverageTangGridVel[iMarker][nSpanWiseSections[marker_flag-1]]             += AverageTangGridVel[iMarker][iSpan]/nSpanWiseSections[marker_flag-1];
            SpanArea[iMarker][nSpanWiseSections[marker_flag-1]]                       += SpanArea[iMarker][iSpan];
            for(iDim=0; iDim < nDim; iDim++){
              AverageTurboNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim]     += AverageTurboNormal[iMarker][iSpan][iDim];
              AverageNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim]          += AverageNormal[iMarker][iSpan][iDim];
              AverageGridVel[iMarker][nSpanWiseSections[marker_flag-1]][iDim]         += AverageGridVel[iMarker][iSpan][iDim]/nSpanWiseSections[marker_flag-1];

            }
          }
        }
      }
    }
  }

  /*--- Normalize 1D normals---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
      if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
          turboNormal2 = 0.0;
          Normal2     = 0.0;

          for (iDim = 0; iDim < nDim; iDim++){
            turboNormal2 += AverageTurboNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim]*AverageTurboNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim];
            Normal2      += AverageNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim]*AverageNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim];
          }
          for (iDim = 0; iDim < nDim; iDim++){
            AverageTurboNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim] /=sqrt(turboNormal2);
            AverageNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim] /=sqrt(Normal2);
          }
        }
      }
    }
  }


  delete [] TotalTurboNormal;
  delete [] TotalNormal;
  delete [] TotalGridVel;
  delete [] TurboNormal;
  delete [] Normal;

}
