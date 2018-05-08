//! @file CellularAutomaton.cpp
//! @brief CellularAutomaton functions

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "CellularAutomaton.h"

//! To keep track of simulation TIME globally.
extern int GLOBALTIMESTEP;
extern bool VesselRemoval;
double stresstau(double P);

//-------------------------------------------------------------------
//! CA constructor.
//! This function creates an empty CA.
CellularAutomaton::CellularAutomaton()
{
	TimeUpdateCells = 0;
	TimeCellMovement = 0;
	TimeRunStep = 0;
}

//-------------------------------------------------------------------
//! CA destructor.
//! This function destroys a CA, writing timing info and freeing all Cells.
CellularAutomaton::~CellularAutomaton()
{
	printf("@@@@@@ FINAL CellularAutomaton TIMING @@@@@@:\n");
	printf("\t Time in RunStep = %g seconds\n",TimeRunStep);
	printf("\t Time in CellMovement = %g seconds\n",TimeCellMovement);
	printf("\t Time in UpdateCells = %g seconds\n",TimeUpdateCells);

	DeleteMemory();
}

//-------------------------------------------------------------------
//! Adds the pressure at a Node to its CA location
//! @todo NOT USED!
void CellularAutomaton::AddVesselPressureToCA()
{
	//Set all pressures to 0.0, to ensure that there are no pressures at points where the vessel was deleted. (Not sure if this is needed)!
	for(int i=0; i<param_GridSizeX; ++i)
	{
		for(int j=0; j<param_GridSizeY; ++j)
		{
			for(int k=0; k<param_GridSizeZ; ++k)
			{
				Cells[i][j][k][0].Pressure=0.0;
			}
		}
	}

	//Add the pressures at each node of the vascular system to their associated CA element
	for(int k=0;k<Vasc.ActiveVessels.size();++k)
	{
		int j = Vasc.ActiveVessels[k];
		int startnode = Vasc.Vessels[j].Start;
		int endnode = Vasc.Vessels[j].End;
		Cells[Vasc.Nodes[startnode].x][Vasc.Nodes[startnode].y][Vasc.Nodes[startnode].z][0].Pressure=Vasc.Nodes[startnode].Pressure;
		Cells[Vasc.Nodes[endnode].x][Vasc.Nodes[endnode].y][Vasc.Nodes[endnode].z][0].Pressure=Vasc.Nodes[endnode].Pressure;
	}
}

//-------------------------------------------------------------------
void CellularAutomaton::AddVesselSiteDiffusibleToCA(int v, int curx, int cury, int curz)
{
	double twopiRL_dx3 = twopi*Vasc.Vessels[v].Radius*Vasc.Vessels[v].Length/pow(param_gridsize,3);
	// if there is no flow in the vessel, it does not contribute to diffusibles
	if (fabs(Vasc.Vessels[v].Flow) < param_NoFlowThreshold) twopiRL_dx3 = 0.0;

	if (param_isOxygenDiffusible) Oxygen.SetVascConstSource(curx,cury,curz,twopiRL_dx3*Oxygen.param_OxygenHaemScaling*Vasc.Vessels[v].Haematocrit);
	if (param_isDrugDiffusible) Drug.SetVascConstSource(curx,cury,curz,Vasc.DrugBloodConcentration*twopiRL_dx3*(1.0-Vasc.Vessels[v].Haematocrit));
	if (param_isProDrugDiffusible) ProDrug.SetVascConstSource(curx,cury,curz,Vasc.ProdrugBloodConcentration*twopiRL_dx3*(1.0-Vasc.Vessels[v].Haematocrit));

	SetVascPropSource(curx,cury,curz,twopiRL_dx3);
}

//-------------------------------------------------------------------
void CellularAutomaton::AddVesselSiteToCA(int v, int curx, int cury, int curz)
{
	int k;
	Location pos = Location(curx,cury,curz,0);

	// when location 0 is occupied by non-vessel:
	// move its cell up the stack, or kill the cell
	if (Cells[curx][cury][curz][0].Type != NONE && Cells[curx][cury][curz][0].Type != VESSEL)
	{
		bool foundspace = false;
		k = 1;
		while (k<StackCells && foundspace==false)
		{
			if (Cells[curx][cury][curz][k].Type == NONE)
				foundspace = true;
			else
				k++;
		}
		// if the cell being moved is part of a protovessel, we need
		// to make sure its xyzlist is updated.
		// If  there is NO room (foundspace==false) we must replace some
		// other cell type with the ENDOCELL or TIPCELL
		if ((foundspace == false) && ((Cells[curx][cury][curz][0].Type == ENDOCELL) || (Cells[curx][cury][curz][0].Type == TIPCELL)))
		{
			k = 1;
			while (k<StackCells && foundspace==false)
			{
				if (Cells[curx][cury][curz][k].Type != ENDOCELL && Cells[curx][cury][curz][k].Type != TIPCELL)
					foundspace = true;
				else
					k++;
			}
			if (foundspace == false)
			{
				printf("ERROR: No space to move ENDOCELL - stack is full of them already!\n");
				exit(1);
			}
		}

		if (foundspace==true)
		{
			Cells[curx][cury][curz][k] = Cells[curx][cury][curz][0];
			// if the cell being moved is part of a protovessel, we need
			// to make sure its xyzlist is updated: this should fix it
			// unless there is NO room (foundspace==false). In this case
			// we should probably replace some other cell type with the
			// ENDOCELL or TIPCELL
			if (Cells[curx][cury][curz][0].Type == ENDOCELL || Cells[curx][cury][curz][0].Type == TIPCELL)
			{
				int PVID = Cells[curx][cury][curz][0].VesselID;
				std::vector<Location>::iterator it;
				it = find (Vasc.ProtoVessels[PVID].xyzstacklist.begin(), Vasc.ProtoVessels[PVID].xyzstacklist.end(), pos);
				Location pvpos = *it;
				pvpos.stack = k;
				*it = pvpos;
			}
			// ENDOCELLs are not active, otherwise:
			if (Cells[curx][cury][curz][0].Type != ENDOCELL)
			{
				std::vector<Location>::iterator it = std::find(ActiveCells.begin(), ActiveCells.end(), pos);
				int curActiveCell = it - ActiveCells.begin();
				pos.stack = k;
				if (!(it == ActiveCells.end()))
				ActiveCells[curActiveCell] = pos;
			}
			if (VERBOSE) printf("$$$$$$$$$ Cell type %d moved from vessel location at (%d,%d,%d) to stack=%d.\n",Cells[curx][cury][curz][0].Type,curx,cury,curz,k);
		}
		else
		{
			KillCell(curx,cury,curz,0);
			printf("$$$$$$$$$ WARNING: No space to move cell from vessel location at (%d,%d,%d).\n",curx,cury,curz);
		}
	}

	// overwrite whatever was at stack=0:
	Cells[curx][cury][curz][0].Type = VESSEL;
	Cells[curx][cury][curz][0].VesselID = v;
	AddVesselSiteDiffusibleToCA(v, curx, cury, curz);
}

//-------------------------------------------------------------------
//! Adds the vasculature to the CellularAutomaton CA.  Calls AddVesselToCA
//! for every vessel in the Vasculature.  Note that this function writes the
//! vasculature over the current Cells array in the CA, meaning that any
//! elements of Cells containing part of the vasculature have their contents
//! overwritten.
void CellularAutomaton::AddVesselsToCA()
{
	int i, k;
	for(k=0;k<Vasc.ActiveVessels.size();k++)
	{
		i = Vasc.ActiveVessels[k];
		AddVesselToCA(i);
	}
}

//-------------------------------------------------------------------
//! Adds a single vessel, vessel (given as an index into the Vessels array),
//! to the CellularAutomaton CA. This is done by adding cells of type VESSEL,
//! at the Start and End of the vessel, to CA's Cells grid.
//!
//! This function writes vessels to location 0 in the stack and moves
//! any other cell type at 0 to the first available stack position if there is
//! one, otherwise destroying it and posting a WARNING.
//!
//! The function also calls AddVesselSiteToCA which updates the vessel Oxygen,
//! Drug, VEGF source terms (VascPropSource, VascConstSource) according to the
//! haematocrit of the vessel, and hence should be called at each time step.
//!
//! TO DO: this should be split so that the add vessel part is separate to
//! passing this oxygen info. Currently AddVesselToCA has to be called twice
//! in each Update for any new vessel, since the flow calculation needs the
//! local info on haematocrit, etc. In addition, AddToCA at the end of Update
//! adds existing Vessels each time as well as updating their properties.
//!
void CellularAutomaton::AddVesselToCA(int v)
{

	int curx, cury, curz;

	// the startnode is a Node
	int startnode = Vasc.Vessels[v].Start;
	curx = Vasc.Nodes[startnode].x;
	cury = Vasc.Nodes[startnode].y;
	curz = Vasc.Nodes[startnode].z;
	AddVesselSiteToCA(v, curx, cury, curz);
	// recall that the NumBonds property records the number of vessels connected to a Node
	// and is hence zero if the CA location is NOT a Node
	std::vector<int>::iterator itIn = find (Vasc.InflowNode.begin(), Vasc.InflowNode.end(), startnode);
	std::vector<int>::iterator itOut = find (Vasc.OutflowNode.begin(), Vasc.OutflowNode.end(), startnode);
	if (itIn != Vasc.InflowNode.end() | itOut != Vasc.OutflowNode.end())
		Cells[curx][cury][curz][0].NumBonds = 3;
	else
		{
		Cells[curx][cury][curz][0].NumBonds = Vasc.Nodes[startnode].NumVesselsIn + Vasc.Nodes[startnode].NumVesselsOut;
		if (taillocs[curx][cury][curz] > 0) Cells[curx][cury][curz][0].NumBonds++;
		}

	if (Cells[curx][cury][curz][0].NumBonds > 3)
		printf("AddVesselToCA: Node %d, start of vessel %d at (%d,%d,%d), NumBonds = %d > 3\n",startnode, v, curx, cury, curz, Cells[curx][cury][curz][0].NumBonds);

	// the endnode is a Node
	int endnode = Vasc.Vessels[v].End;
	curx = Vasc.Nodes[endnode].x;
	cury = Vasc.Nodes[endnode].y;
	curz = Vasc.Nodes[endnode].z;
	AddVesselSiteToCA(v, curx, cury, curz);
	itIn = find (Vasc.InflowNode.begin(), Vasc.InflowNode.end(), endnode);
	itOut = find (Vasc.OutflowNode.begin(), Vasc.OutflowNode.end(), endnode);
	if (itIn != Vasc.InflowNode.end() | itOut != Vasc.OutflowNode.end())
		Cells[curx][cury][curz][0].NumBonds = 3;
	else
		{
		Cells[curx][cury][curz][0].NumBonds = Vasc.Nodes[endnode].NumVesselsIn + Vasc.Nodes[endnode].NumVesselsOut;
		if (taillocs[curx][cury][curz] > 0) Cells[curx][cury][curz][0].NumBonds++;
		}

	if (Cells[curx][cury][curz][0].NumBonds > 3)
		printf("AddVesselToCA: Node %d, end of vessel %d at (%d,%d,%d), NumBonds = %d > 3\n",endnode, v, curx, cury, curz, Cells[curx][cury][curz][0].NumBonds);
}

//-------------------------------------------------------------------
//! This function ages all cells by TimeStep.
//! It iterates through all dividing cells and increments Age.
//! If a cell is a MACROPHAGE that is too old, it is killed.
void CellularAutomaton::AgeCells(double TimeStep)
{
	int NumActiveCells = ActiveCells.size();
	int i, count=0;
	Location pos;

	for(i=0;i<ActiveCells.size();i++)
    {
		pos = ActiveCells[i];
		Cells[pos.x][pos.y][pos.z][pos.stack].Age += TimeStep;
		if (Cells[pos.x][pos.y][pos.z][pos.stack].Type==MACROPHAGE || Cells[pos.x][pos.y][pos.z][pos.stack].Type==MONOCYTE)
		{
			if (Cells[pos.x][pos.y][pos.z][pos.stack].Age > param_MacLifespan)
			{
				KillCell(pos.x,pos.y,pos.z,pos.stack);
				printf("Killed Macrophage/Monocyte at (%d,%d,%d). Active cells %zu, i=%d\n",pos.x,pos.y,pos.z,ActiveCells.size(),i);
				i--;
			}
		}
		count++;
    }

	if (count != NumActiveCells)
		printf("Active cells at end %zu, i=%d, cells checked=%d\n",ActiveCells.size(),i,count);

}

//-------------------------------------------------------------------
//! Find new vessel connections.
//! Loops over all ProtoVessels and checks if their tips coincide with other (Proto)Vessels.
//! If so, create new Nodes and new Vessels.
bool CellularAutomaton::Anastomosis()
{
	// need these when setting up new vessels
	double radius = 2.3456e-4;
	bool NewVessels = false;
	int NewNode, NewVesselID;

	//////////////////////////////
	// loop through all ProtoVessels
	for(int i=0;i<Vasc.NumProtoVessels;i++)
	{
		// don't do anything if the ProtoVessel is already connected (i.e. inactive)
		if (Vasc.ProtoVessels[i].Connected == false)
		{
			// age protovessels
			Vasc.ProtoVessels[i].Age += param_timestep;
			if (Vasc.ProtoVessels[i].Age > param_ProtoVesselDeathTime)
			{
				// remove the current ProtoVessel
				printf("Sprout %d too old - removing it ....\n",i);
				Vasc.ProtoVessels[i].Connected = true;
				RemoveProtoVesselFromCA(i);
				Vasc.AvailableProtoVessels.push_back(i);
			}
		}
	}

	//////////////////////////////
	// loop through all ProtoVessels
	for(int i=0;i<Vasc.NumProtoVessels;i++)
	{
		// don't do anything if the ProtoVessel is already connected!
		if (Vasc.ProtoVessels[i].Connected == false)
		{
			// number of entries in the list for this ProtoVessel
			int nxyzstacklist = Vasc.ProtoVessels[i].xyzstacklist.size();

			// only check for connections if there is at least param_minProtoVesselLengthForAnastomosis
			// members of the ProtoVessel. i.e. check that the ProtoVessel has actually gone somewhere
			if (nxyzstacklist>=param_minProtoVesselLengthForAnastomosis)
			{
				// location of the tip of the current ProtoVessel
				Location tiploc = Vasc.ProtoVessels[i].xyzstacklist[nxyzstacklist-1];
				// loop through the local cell stack checking for ENDOCELLs or VESSELs
				// always find VESSELs first since they are always at stack=0
				for(int j=0;j<StackCells;j++)
				{
					// This is needed because the ProtoVessel may get connected for an earlier j in the stack
					if (Vasc.ProtoVessels[i].Connected == false)
					{
						if ((Cells[tiploc.x][tiploc.y][tiploc.z][j].Type==ENDOCELL) || (Cells[tiploc.x][tiploc.y][tiploc.z][j].Type==TIPCELL))
						{
							// ID of ProtoVessel connected to
							int PVID = Cells[tiploc.x][tiploc.y][tiploc.z][j].VesselID;
							if (PVID == i && Cells[tiploc.x][tiploc.y][tiploc.z][j].Type==ENDOCELL)
								{
								// if the connection is to the same ProtoVessel
								// connect the current ProtoVessel to itself
								printf("Sprout %d connected to itself ....\n",i);
								Vasc.ProtoVessels[i].Connected = true;
								RemoveProtoVesselFromCA(i);
								Vasc.AvailableProtoVessels.push_back(i);
								////////////////////////////////////////////////////////////
								// testing sprout-sprout connection
								////////////////////////////////////////////////////////////
								//            Y<<<<<<<<
								// T>>>>>>>>>X        ^
								//           >>>>>>>>>>
								//
								// T marks the tail of the sprout (node nn1)
								// X marks the tip AND the point of anastomosis, Y the last but one step (node nn2)
								int OldNumVessels = Vasc.NumVessels;
								int nn1 = Vasc.ProtoVessels[i].TailNode;

								// assign TWSSsubthreshold=Age of the sprout
								// make vessels for each CA size segment
								// only go as far as the last but one location
								int nn2, nnx;
								for (int k=1; k<nxyzstacklist-1; k++)
								{
									nn2 = Vasc.MakeNewNode(Vasc.ProtoVessels[i].xyzstacklist[k].x,Vasc.ProtoVessels[i].xyzstacklist[k].y,Vasc.ProtoVessels[i].xyzstacklist[k].z);
									NewVesselID = Vasc.MakeNewVessel(nn1, nn2, param_InflowHaematocrit, radius, 0.0, 0.0, Vasc.ProtoVessels[i].Age);
									AddVesselToCA(NewVesselID);
									// the new Node must inherit its connectivity from the parent vessel
									Vasc.Nodes[nn2].Connected = Vasc.Nodes[nn1].Connected;
									nn1 = nn2;
									if (Vasc.ProtoVessels[i].xyzstacklist[k].x==tiploc.x && Vasc.ProtoVessels[i].xyzstacklist[k].y==tiploc.y && Vasc.ProtoVessels[i].xyzstacklist[k].z==tiploc.z)
									{
										// if the current point in the sprout is the same as the tip,
										// this is the node to connect to later (node X above)
										nnx = Vasc.Vessels[NewVesselID].End;
									}
								}

								// now make one last vessel from the last node created (node Y above, nn2) to the anastomosis node (X, now nn1)
								// insert a new connection from nn1 to nn2 (X to Y above)
								// assign TWSSsubthreshold=Age of the sprout
								NewVesselID = Vasc.MakeNewVessel(nnx, nn2, param_InflowHaematocrit, radius, 0.0, 0.0, Vasc.ProtoVessels[i].Age);
								AddVesselToCA(NewVesselID);

								printf("New vessels created, NumVessels = %d\n",Vasc.NumVessels);
								NewVessels = true;
								////////////////////////////////////////////////////////////
								////////////////////////////////////////////////////////////
								}
							// if the connection is to a different ProtoVessel
							else if (PVID != i)
								{
								printf("Tip-sprout connection %d from (%d,%d,%d) ->%d at (%d,%d,%d)..............\n",i,Vasc.ProtoVessels[i].xyzstacklist[0].x,Vasc.ProtoVessels[i].xyzstacklist[0].y,Vasc.ProtoVessels[i].xyzstacklist[0].z,PVID,tiploc.x,tiploc.y,tiploc.z);
								// make the tip cell inactive
								ActiveCells.erase(std::remove(ActiveCells.begin(), ActiveCells.end(), tiploc), ActiveCells.end());
								// (tipx,tipy,tipz,tail node) for sprout
								int n1x = Vasc.ProtoVessels[i].xyzstacklist[0].x;
								int n1y = Vasc.ProtoVessels[i].xyzstacklist[0].y;
								int n1z = Vasc.ProtoVessels[i].xyzstacklist[0].z;
								int nn1 = Vasc.ProtoVessels[i].TailNode;
								// (tipx,tipy,tipz,tail node) for ProtoVessel connected to
								int n2x = Vasc.ProtoVessels[PVID].xyzstacklist[0].x;
								int n2y = Vasc.ProtoVessels[PVID].xyzstacklist[0].y;
								int n2z = Vasc.ProtoVessels[PVID].xyzstacklist[0].z;
								int nn2 = Vasc.ProtoVessels[PVID].TailNode;
								// (x,y,z) for anastomosis location
								int nax = tiploc.x;
								int nay = tiploc.y;
								int naz = tiploc.z;
								printf("(x1,y1,z1,nn1)=(%d,%d,%d,%d), (x2,y2,z2,nn2)=(%d,%d,%d,%d), (xa,ya,za)=(%d,%d,%d)\n", n1x, n1y, n1z, nn1, n2x, n2y, n2z, nn2, nax, nay, naz);
								if (Vasc.ProtoVessels[PVID].Connected == false)
									{
									int OldNumVessels = Vasc.NumVessels;

									// remove the current ProtoVessel
									// do this prior to adding any new vessels, since that potentially
									// moves ENDOCELLS up in the cell stack in order to place VESSELs
									// at stack=0.
									Vasc.ProtoVessels[i].Connected = true;
									RemoveProtoVesselFromCA(i);

									Location aloc = tiploc;
									aloc.stack = j;
									// starting index for the remainder of the protovessel
									std::vector<Location>::iterator it, itPVIDstart;
									itPVIDstart = find (Vasc.ProtoVessels[PVID].xyzstacklist.begin(), Vasc.ProtoVessels[PVID].xyzstacklist.end(), aloc);

									// decide what to do with the ProtoVessel connected to:
									if (Cells[tiploc.x][tiploc.y][tiploc.z][j].Type==TIPCELL)
										{
										// both ProtoVessels removed if a Tip-Tip connection
										printf("TIP-TIP...\n");
										Vasc.ProtoVessels[PVID].Connected = true;
										RemoveProtoVesselFromCA(PVID);
										}
									else
										{
										// reset the tail to zero
										taillocs[Vasc.ProtoVessels[PVID].xyzstacklist[0].x][Vasc.ProtoVessels[PVID].xyzstacklist[0].y][Vasc.ProtoVessels[PVID].xyzstacklist[0].z]=0;
										// remove elements from CA
										for(it = Vasc.ProtoVessels[PVID].xyzstacklist.begin(); !(*it == aloc); it++)
											{
											Location pos = *it;
											// check that it really is an ENDOCELL
											if ((Cells[pos.x][pos.y][pos.z][pos.stack].Type == ENDOCELL) || (Cells[pos.x][pos.y][pos.z][pos.stack].Type == TIPCELL))
												{
												Cells[pos.x][pos.y][pos.z][pos.stack].Type = NONE;
												ActiveCells.erase(std::remove(ActiveCells.begin(), ActiveCells.end(), pos), ActiveCells.end());
												}
											}
										}

									// insert a new connection from nn1 to nn2 (between the tails of the two sprouts)
									// assign TWSSsubthreshold=Age of the primary sprout

									// make new nodes and vessels along each step of ProtoVessel i
									for (int k=1; k<nxyzstacklist; k++)
									{
										nn2 = Vasc.MakeNewNode(Vasc.ProtoVessels[i].xyzstacklist[k].x,Vasc.ProtoVessels[i].xyzstacklist[k].y,Vasc.ProtoVessels[i].xyzstacklist[k].z);
										//int test; cin>>test;
										NewVesselID = Vasc.MakeNewVessel(nn1, nn2, param_InflowHaematocrit, radius, 0.0, 0.0, Vasc.ProtoVessels[i].Age);
										AddVesselToCA(NewVesselID);
										// the new Node must inherit its connectivity from the parent vessel
										Vasc.Nodes[nn2].Connected = Vasc.Nodes[nn1].Connected;
										nn1 = nn2;
									}
									int AnastomosisNode = Vasc.Vessels[NewVesselID].End;

									// make new nodes and vessels along each step of ProtoVessel PVID
									// starting at aloc, the location of anastomosis
									// assign TWSSsubthreshold=Age of the secondary sprout
									for(it = itPVIDstart-1; it != Vasc.ProtoVessels[PVID].xyzstacklist.begin(); it--)
									{
										Location pos = *it;
										nn2 = Vasc.MakeNewNode(pos.x,pos.y,pos.z);
										//int test; cin>>test;
										NewVesselID = Vasc.MakeNewVessel(nn1, nn2, param_InflowHaematocrit, radius, 0.0, 0.0, Vasc.ProtoVessels[PVID].Age);
										AddVesselToCA(NewVesselID);
										// the new Node must inherit its connectivity from the parent vessel
										Vasc.Nodes[nn2].Connected = Vasc.Nodes[nn1].Connected;
										nn1 = nn2;
									}
									// the last step connects to the existing node at the tail of PVID
									NewVesselID = Vasc.MakeNewVessel(nn1, Vasc.ProtoVessels[PVID].TailNode, param_InflowHaematocrit, radius, 0.0, 0.0, Vasc.ProtoVessels[PVID].Age);
									AddVesselToCA(NewVesselID);

									// decide what to do with the ProtoVessel connected to:
									if (Vasc.ProtoVessels[PVID].Connected == false)
										{
										// remove elements from xyzlist up to the connection point
										Vasc.ProtoVessels[PVID].xyzstacklist.erase(Vasc.ProtoVessels[PVID].xyzstacklist.begin(),itPVIDstart);
										Vasc.ProtoVessels[PVID].TailNode = AnastomosisNode;
										// move the tail
										taillocs[Vasc.ProtoVessels[PVID].xyzstacklist[0].x][Vasc.ProtoVessels[PVID].xyzstacklist[0].y][Vasc.ProtoVessels[PVID].xyzstacklist[0].z]=PVID+1;
										}

									printf("New vessels created, NumVessels = %d\n",Vasc.NumVessels);
									NewVessels = true;
									}

									// FINALLY: check that the original tail nodes are connected
									// they may have been disconnected by a prior removal in
									// RemoveLowFlowVessels.
									// FOR NOW DO THIS GLOBALLY AT END OF THIS FUNCTION
								}
						}
						else if (Cells[tiploc.x][tiploc.y][tiploc.z][j].Type==VESSEL)
						{
							// if the tip coincides with a VESSEL
							// we want to split that vessel at the point of anastomosis
							// and make a new vessel connecting that new node with the
							// start point of the ProtoVessel

							// make the tip cell inactive
							ActiveCells.erase(std::remove(ActiveCells.begin(), ActiveCells.end(), tiploc), ActiveCells.end());
							// ID of Vessel connected to
							int VID = Cells[tiploc.x][tiploc.y][tiploc.z][j].VesselID;
							printf("Tip-vessel connection %d from (%d,%d,%d) ->vessel %d at (%d,%d,%d)..............\n",i,Vasc.ProtoVessels[i].xyzstacklist[0].x,Vasc.ProtoVessels[i].xyzstacklist[0].y,Vasc.ProtoVessels[i].xyzstacklist[0].z,VID,tiploc.x,tiploc.y,tiploc.z);
							// Node at the start of the ProtoVessel
							int nn1 = Vasc.ProtoVessels[i].TailNode;
							int nn2;

							// remove the current ProtoVessel
							// do this prior to adding any new vessels, since that potentially
							// moves ENDOCELLS up in the cell stack in order to place VESSELs
							// at z=0.
							Vasc.ProtoVessels[i].Connected = true;
							RemoveProtoVesselFromCA(i);

							////////////////////////////////////////
							// The vessel connected to should already have a node at tiploc
							////////////////////////////////////////
							int nax = Vasc.Nodes[Vasc.Vessels[VID].Start].x;
							int nay = Vasc.Nodes[Vasc.Vessels[VID].Start].y;
							int naz = Vasc.Nodes[Vasc.Vessels[VID].Start].z;
							// either anastomosis is at the start of vessel VID or at the end
							if (nax == tiploc.x && nay == tiploc.y && naz == tiploc.z)
								nn2 = Vasc.Vessels[VID].Start;
							else
								nn2 = Vasc.Vessels[VID].End;

							int OldNumVessels = Vasc.NumVessels;
							// insert a new connection from nn1 to nn2
							// assign TWSSsubthreshold=Age of the primary sprout
							// make new nodes and vessels from nn1 along each step of ProtoVessel i
							int nni;
							for (int k=1; k<nxyzstacklist-1; k++)
							{
								nni = Vasc.MakeNewNode(Vasc.ProtoVessels[i].xyzstacklist[k].x,Vasc.ProtoVessels[i].xyzstacklist[k].y,Vasc.ProtoVessels[i].xyzstacklist[k].z);
								NewVesselID = Vasc.MakeNewVessel(nn1, nni, param_InflowHaematocrit, radius, 0.0, 0.0, Vasc.ProtoVessels[i].Age);
								AddVesselToCA(NewVesselID);
								nn1 = nni;
							}
							// final connection is to nn2
							NewVesselID = Vasc.MakeNewVessel(nn1, nn2, param_InflowHaematocrit, radius, 0.0, 0.0, Vasc.ProtoVessels[i].Age);
							AddVesselToCA(NewVesselID);

							printf("New vessels created, NumVessels = %d\n",Vasc.NumVessels);
							NewVessels = true;
						}
					}
				}
			}
		}
	}
	if (NewVessels)
	{
		Vasc.UpdateVesselLengths();
	}
	return NewVessels;
}

//-------------------------------------------------------------------
//! Calculate summary statistics during or after simulation.
//! This will update the CellularAutomaton variables TotalCancer,
//! TotalVessel, MeanVEGF, etc.
void CellularAutomaton::CalculateStatistics()
{
	int x, y, z, stack;
	int tempEmpty;
	TotalEmpty = 0;
	TotalNormal = 0;
	TotalCancer = 0;
	TotalQuiescent = 0;
	TotalMacs = 0;
	TotalMonos = 0;
	TotalVessel = 0;
	MeanDrug = 0;
	MeanOxygen = 0;
	MeanProDrug = 0;
	MeanVEGF = 0;

	for(x=0;x<XCells;x++)
	{
		for(y=0;y<YCells;y++)
		{
			for(z=0;z<ZCells;z++)
			{
				tempEmpty = 0;
				for(stack=0;stack<StackCells;stack++)
				{
					if (Cells[x][y][z][stack].Type==CANCER)
					{
						if (Cells[x][y][z][stack].Quiescent)
							TotalQuiescent++;
						else
							TotalCancer++;
					}
					else if (Cells[x][y][z][stack].Type==VESSEL)
					{
						TotalVessel++;
					}
					else if (Cells[x][y][z][stack].Type==NORMALCELL && Cells[x][y][z][stack].Phenotype==AEROBIC)	TotalNormal++;
					else if (Cells[x][y][z][stack].Type==MACROPHAGE)	TotalMacs++;
					else if (Cells[x][y][z][stack].Type==MONOCYTE)		TotalMonos++;
					else if (Cells[x][y][z][stack].Type==NONE)			tempEmpty++;
				}
				if (tempEmpty == StackCells) TotalEmpty++;
				MeanDrug += Drug.GetLevelDirect(x,y,z);
				drugsum[x][y][z] += Drug.GetLevelDirect(x,y,z);
				MeanOxygen += Oxygen.GetLevelDirect(x,y,z);
				MeanProDrug += ProDrug.GetLevelDirect(x,y,z);
				MeanVEGF += VEGF.GetLevelDirect(x,y,z);
			}
		}
	}

	MeanDrug = MeanDrug/(double)(XCells*YCells*ZCells);
	MeanOxygen = MeanOxygen/(double)(XCells*YCells*ZCells);
	MeanProDrug = MeanProDrug/(double)(XCells*YCells*ZCells);
	MeanVEGF = MeanVEGF/(double)(XCells*YCells*ZCells);
	if (VERBOSE) printf("TotalCancer = %d, TotalVessel = %d, MeanVEGF = %f\n",TotalCancer,TotalVessel,MeanVEGF);
}

//-------------------------------------------------------------------
//! This function examines CA's Diffusible grids along the paths of each of the
//! blood vessels (using CalculateVesselDiffusible) and determines the
//! average level of each relevant Diffusible along the path of the
//! blood vessel.  This is then stored in the Vessel's MeanDiffusible variable.
void CellularAutomaton::CalculateVascularDiffusibles()
{
	int i,k;

	for(k=0;k<Vasc.ActiveVessels.size();k++)
	{
		i = Vasc.ActiveVessels[k];
		if (param_isVEGFDiffusible) Vasc.Vessels[i].MeanVEGF = Vasc.CalculateVesselDiffusible(VEGF, i);
		else Vasc.Vessels[i].MeanVEGF = 0;
	}
}

//-------------------------------------------------------------------
//! Calculate the vessel radii, using information about cells and diffusibles
double CellularAutomaton::CalculateVascularRadii()
{
	Vasc.mystart=clock();

	// This version includes the CA as input as cells affect Vessel behaviour
	// e.g. via co-option
	// Calculates the vessel radii and returns the maximum (absolute) change in radius since the last update
	double MaxProportionalChange=0;
	int i,j,k;
	double DownStreamStimulus;
	double UpStreamStimulus;
	double bigexpression;
	int newstate;

	for(k=0;k<Vasc.ActiveVessels.size();k++)
	{
		i = Vasc.ActiveVessels[k];

		if (fabs(Vasc.Vessels[i].Flow) > param_NoFlowThreshold)
		{
			//Calculate the average pressure in the vessel
			double P = (Vasc.Nodes[Vasc.Vessels[i].Start].Pressure + Vasc.Nodes[Vasc.Vessels[i].End].Pressure)/2.0;
			double H = Vasc.Vessels[i].Haematocrit;
			double VEGFLevel = Vasc.Vessels[i].MeanVEGF;
			double mult;

 			int NumTrueIn=0;
			double Jm=0.;
			double Jmd=0.;

			double affinityratio = 10.;
			// scaling mult due to anti-angiogenic drug may be equivalent to
			mult = 1+param_VEGFMaxMetabolic*VEGFLevel/(VEGFLevel+param_V0);

			double km = mult*param_km;
			double km2 = param_km2;
			double kc = param_upstreamkc;
			double J0 = param_upstreamJ0;

			for(j=0;j<Vasc.Nodes[Vasc.Vessels[i].Start].NumVesselsIn;j++)
			{
				if(Vasc.Vessels[Vasc.Nodes[Vasc.Vessels[i].Start].VesselsIn[j]].Flow>0.0)
				{
					if(Vasc.Vessels[Vasc.Nodes[Vasc.Vessels[i].Start].VesselsIn[j]].MeanVEGF>0.)
					{
						Jm+=0.01*Vasc.Vessels[Vasc.Nodes[Vasc.Vessels[i].Start].VesselsIn[j]].Flow;
						//Vasc.Vessels[Vasc.Nodes[Vasc.Vessels[i].Start].VesselsIn[j]].Length;
					}
				}
			}

			for(j=0;j<Vasc.Nodes[Vasc.Vessels[i].End].NumVesselsOut;j++)
			{
				if(Vasc.Vessels[Vasc.Nodes[Vasc.Vessels[i].End].VesselsOut[j]].MeanVEGF>0.)
					Jmd+=Vasc.Vessels[Vasc.Nodes[Vasc.Vessels[i].End].VesselsOut[j]].Length/pow(2.7182818,Vasc.Vessels[Vasc.Nodes[Vasc.Vessels[i].End].VesselsOut[j]].Length);
			}

			UpStreamStimulus=km2*kc*Jmd/(Jmd+J0);
			DownStreamStimulus=km2*log(1+(Jm/(fabs(Vasc.Vessels[i].Flow)+param_Qdotref2)));

			if(Vasc.Vessels[i].WallStress+param_Tauref <= 0) printf("Wallstress+Tauref=%f\n",Vasc.Vessels[i].WallStress+param_Tauref);
			if(stresstau(P) <= 0) printf("stresstau=%f\n",stresstau(P));
			double Stimulus_WallShearStress = log10(Vasc.Vessels[i].WallStress+param_Tauref);
			double Stimulus_Pressure = param_kp*log10(stresstau(P));
			double Stimulus_Metabolic = km*log10(1.+param_Qdotref/(fabs(Vasc.Vessels[i].Flow)*H));

			newstate = Vasc.Vessels[i].State;
			if (newstate == COLLAPSEVESSEL)
			{
				double LCOLLAPSE=param_CollapseScaling*param_ks;
				if(VesselRemoval == false) bigexpression = Stimulus_WallShearStress - Stimulus_Pressure - param_ks - LCOLLAPSE;
			}
			else if (newstate == ANGIOVESSEL)
			{
				// ``Angiogenesis''
				bigexpression = Stimulus_WallShearStress - Stimulus_Pressure + Stimulus_Metabolic - param_ks;
			}
			else
			{
				bigexpression = Stimulus_WallShearStress - Stimulus_Pressure
					+ Stimulus_Metabolic + DownStreamStimulus + UpStreamStimulus - param_ks;
			}

			if(isnan(bigexpression))
			{
				printf("ERROR: Bigexpression = NaN\n");
				exit(1);
			}
		}
		else
		{
			bigexpression = (param_RMIN-Vasc.Vessels[i].Radius)/(Vasc.Vessels[i].Radius*param_VesselRadiusTimestep);
		}

		double NewRadius = Vasc.Vessels[i].Radius + param_VesselRadiusTimestep*Vasc.Vessels[i].Radius*bigexpression;
		if (NewRadius > param_RMAX) {
			NewRadius = param_RMAX;
			bigexpression = (NewRadius - Vasc.Vessels[i].Radius)/(param_VesselRadiusTimestep*Vasc.Vessels[i].Radius);
		}
		if (NewRadius < param_RMIN) {
			NewRadius = param_RMIN;
			bigexpression = (NewRadius - Vasc.Vessels[i].Radius)/(param_VesselRadiusTimestep*Vasc.Vessels[i].Radius);
		}
		Vasc.Vessels[i].Radius = NewRadius;

		if(fabs(bigexpression)>MaxProportionalChange) MaxProportionalChange = fabs(bigexpression);

		if(isnan(bigexpression))
		{
			printf("ERROR: Bigexpression = NaN\n");
			exit(1);
		}

	}
	Vasc.mystop=clock();
	Vasc.CVRtime += Vasc.mystop-Vasc.mystart;
	return MaxProportionalChange;
}

//-------------------------------------------------------------------
//! Loop over a Vessel CooptionList (CA locations adjacent to a vessel)
//! accumulating the number of CANCER and NORMALCELL. If the proportion of
//! CANCER is > param_cooptionThreshold, the vessel is coopted.
//!
//! @todo Check if the proportion should take account of other non-cancer cell types.
//!
bool CellularAutomaton::CalculateVesselCooption(int v)
{
	int numcancer=0;
	int numnormal=0;
	bool COOPTED=false;
	Location pos;
	int curx, cury, curz, curstack, k;

	for(int i=0;i<Vasc.Vessels[v].cooptionlist.size();i++)
	{
		pos = Vasc.Vessels[v].cooptionlist[i];
		curx = pos.x;
		cury = pos.y;
		curz = pos.z;
		curstack = pos.stack;
		if(Cells[curx][cury][curz][curstack].Type==CANCER) numcancer += 1;
		else if(Cells[curx][cury][curz][curstack].Type==NORMALCELL) numnormal += 1;
	}
	if ((double)numcancer/((double)numcancer+(double)numnormal) > param_cooptionThreshold) COOPTED=true;
	return COOPTED;
}

//-------------------------------------------------------------------
bool CellularAutomaton::CanDivide(int x, int y, int z, int stack, Location& NewPos)
{
  if(STENCIL == DIRECT4)
    {
      return CanDivideDirect4(x,y,z,stack,NewPos);
    }
  else if(STENCIL == TRANSVESSEL)
    {
      return CanDivideTransVessel(x,y,z,stack,NewPos);
    }
  else if(STENCIL == LOCAL)
    {
      bool found = false;
      found = CanDivideLocal(x,y,z,stack,NewPos);
      if (!found)
	  {
		if (param_RandomDaughterCellPlacementXYZ) {
			found = CanDivideRandom(x,y,z,stack,NewPos);
			if (found && VERBOSE) printf("Division into random neighbour after carrycap reached\n");
		}
		else {
			found = CanDivideTransVessel(x,y,z,stack,NewPos);
			if (found && VERBOSE) printf("Division into neighbour after carrycap reached\n");
		}
	  }
      return found;
    }
  return false;

}

//-------------------------------------------------------------------
//! only consider direct neighbours (right, left, up, down, front, back)
//!
//! to get the old DIRECT4 model just set StackCells=1
//! so no need to check CONTACT_INHIBITION
bool CellularAutomaton::CanDivideDirect4(int x, int y, int z, int stack, Location& NewPos)
{
  int cellcount;
      if(x<XCells-1)
		{
		  cellcount=StackCount(x+1, y, z, NewPos);
		  if (cellcount == 0)
			// only divide into x+1 if it is completely EMPTY!
			{
			  NewPos.x=x+1;
			  NewPos.y=y;
			  NewPos.z=z;
			  NewPos.stack=0;
			  return true;
			}
		}
      if(x>0) // x-1
		{
		  cellcount=StackCount(x-1, y, z, NewPos);
		  if (cellcount == 0)
			{
			  NewPos.x=x-1;
			  NewPos.y=y;
			  NewPos.z=z;
			  NewPos.stack=0;
			  return true;
			}
		}
      if(y<YCells-1) // y+1
		{
		  cellcount=StackCount(x, y+1, z, NewPos);
		  if (cellcount == 0)
			{
			  NewPos.x=x;
			  NewPos.y=y+1;
			  NewPos.z=z;
			  NewPos.stack=0;
			  return true;
			}
		}
      if(y>0) // y-1
		{
		  cellcount=StackCount(x, y-1, z, NewPos);
		  if (cellcount == 0)
			{
			  NewPos.x=x;
			  NewPos.y=y-1;
			  NewPos.z=z;
			  NewPos.stack=0;
			  return true;
			}
		}
	  if(z<ZCells-1) // z+1
        {
		  cellcount=StackCount(x, y, z+1, NewPos);
          if (cellcount == 0)
            {
              NewPos.x=x;
              NewPos.y=y;
              NewPos.z=z+1;
              NewPos.stack=0;
              return true;
            }
        }
      if(z>0) // z-1
        {
		  cellcount=StackCount(x, y+1, z, NewPos);
          if (cellcount == 0)
            {
              NewPos.x=x;
              NewPos.y=y;
              NewPos.z=z-1;
              NewPos.stack=0;
              return true;
            }
        }
      else if(Cells[x][y][z][stack].Type==CANCER)
	// if you didn't find any empty locations
	// and the cell is a cancer cell
	// allow division if number of Cancer cells is less than StackCells
	{
	  if(Cells[x][y][z][stack].Quiescent)
	    return false;
	  else
	    {
		  cellcount = StackCount(x, y, z, NewPos);
	      if (cellcount<StackCells)
		{
		  NewPos.x=x;
		  NewPos.y=y;
		  NewPos.z=z;
		  printf("CANCER CELLS < CAP...\n");
		  return true;
		}
	    }
	}
      //******************************
      // end DIRECT4
      //******************************
}

//-------------------------------------------------------------------
bool CellularAutomaton::CanDivideLocal(int x, int y, int z, int stack, Location& NewPos)
// ************************************************************
// New local stencil, only divide if room in current location
// 160506 Markus Owen
// ************************************************************
{

	int filled=0;
	double mass=0;

	filled = StackCount(x, y, z, NewPos);
	if (filled<Cells[x][y][z][stack].DivCarryCap && filled<StackCells)
	//if (mass<double(Cells[x][y][z][stack].DivCarryCap))
	{
		NewPos.x=x;
		NewPos.y=y;
		NewPos.z=z;
		return true;
	}
	return false;
}

//-------------------------------------------------------------------
//! Check if a dividing cell can place its daughter in a neighbouring site
//! @todo Amalgamate with CanDivideRandom by randomising the start of the while loops if the flags indicate.
bool CellularAutomaton::CanDivideTransVessel(int x, int y, int z, int stack, Location& NewPos)
{
	int cellcount;
	// Making sure cell can't divide further than the size of grid
	int minx = MAX(x-param_DivisionRadius, 0);
	int maxx = MIN(x+param_DivisionRadius, XCells-1);
	int miny = MAX(y-param_DivisionRadius, 0);
	int maxy = MIN(y+param_DivisionRadius, YCells-1);
	int minz = MAX(z-param_DivisionRadius,0);
	int maxz = MIN(z+param_DivisionRadius, ZCells-1);

	if(PERIODIC_X)
	{
		minx = x-param_DivisionRadius;
		maxx = x+param_DivisionRadius;
	}
	if(PERIODIC_Y)
	{
		miny = y-param_DivisionRadius;
		maxy = y+param_DivisionRadius;
	}
	if(PERIODIC_Z)
	{
		minz = z-param_DivisionRadius;
		maxz = z+param_DivisionRadius;
	}

	// To store the locations of the surrounding elements to check
	Location test[3];
	int NumNeighbours = ((maxx-minx+1)*(maxy-miny+1)*(maxz-minz+1));
	int Numtests;
	int count = 0;
	bool found = false;
	double maxoxygen = -1.0;

	int tx = minx;
	int ty = miny;
	int tz = minz;
	// randomise the start of the search if necessary
	if (param_RandomDaughterCellPlacement) {
		if (maxx>minx) tx += mtrand.randInt(maxx-minx);
		if (maxy>miny) ty += mtrand.randInt(maxy-miny);
		if (maxz>minz) tz += mtrand.randInt(maxz-minz); // ensures randInt only called twice for 2D
	}

	// loop round the test elements starting at (tx,ty,tz)
	while(count < NumNeighbours)
	{
		if(PERIODIC_X) { tx=MOD(tx,(XCells));}
		if(PERIODIC_Y) { ty=MOD(ty,(YCells));}
		if(PERIODIC_Z) { tz=MOD(tz,(ZCells));}

		// if test location is a vessel, up to three alternative test
		// locations are used, hence Numtests is set to 1, 2 or 3
		Numtests = TransVesselLocations(x, y, z, tx, ty, tz, test);

		// loop through the locations to test (1 or 2 or 3)
		// choose a random number for random start to loop through the test locations
		int count_test = 0;
		int t = 0;
		// randomise the start of the search if necessary
		if (param_RandomDaughterCellPlacement) t = mtrand.randInt(Numtests-1);
	    while(count_test < Numtests)
		{
			if(PERIODIC_X) { test[t].x=MOD(test[t].x,(XCells));}
			if(PERIODIC_Y) { test[t].y=MOD(test[t].y,(YCells));}
			if(PERIODIC_Z) { test[t].z=MOD(test[t].z,(ZCells));}
			// if test element within grid
			if(test[t].x>=0 && test[t].x<XCells && test[t].y>=0 && test[t].y<YCells && test[t].z>=0 && test[t].z<ZCells)
			{
				// count cells at each location
				cellcount = StackCount(test[t].x, test[t].y, test[t].z, test[t]);

				// if test element has space
				if (cellcount < Cells[x][y][z][stack].DivCarryCap && cellcount < StackCells)
				{
					if (!param_OxygenDependentDaughterCellPlacement)
					{
						  NewPos = test[t];
						  found=true;
						  return found;
					}
					else
					{
						// if the level of oxygen at test element > maxoxygen
						// Then can divide
						if(Oxygen.GetLevelDirect(test[t].x, test[t].y, test[t].z) > maxoxygen)
						{
						  NewPos = test[t];
						  // and update maxoxygen for subsequent comparison
						  maxoxygen = Oxygen.GetLevelDirect(test[t].x, test[t].y, test[t].z);
						  found=true;
						}
					}
				}
			}
			// next transvessel test location
			t++;
			if(t == Numtests) t = 0;
			// increment count_test
			count_test++;
		}
		// move to the next neighbourhood test location
		if (tz == maxz) {
			tz = minz;
			if (ty == maxy) {
				ty = miny;
				if (tx == maxx) tx = minx;
				else tx++;
				}
			else ty++;
			}
		else tz++;
		count++;
    }
	return found;
}

//-------------------------------------------------------------------
//! Alternative to CanDivideTransvessel.
//! Search for daughter placement starts at a random neighbour.
//! If there is an oxygen gradient it should not matter where the search starts.
//! param_RandomDaughterCellPlacement==true means the daughter location is purely random
//! param_RandomDaughterCellPlacement==false means the search starts at random but the
//! daughter is placed at the location with maximum oxygen
//! Because the search is preceded by three calls to randInt, archived simulations
//! will not be exactly reproduced.
bool CellularAutomaton::CanDivideRandom(int x, int y, int z, int stack, Location& NewPos)
{
  // ************************************************************
  // RANDOM
  // ************************************************************
  int cellcount;
  //printf("DOING CanDivide with RANDOM\n");
  // Making sure cell can't divide further than the size of grid
  int minx = MAX(x-param_DivisionRadius, 0);
  int maxx = MIN(x+param_DivisionRadius, XCells-1);
  int miny = MAX(y-param_DivisionRadius, 0);
  int maxy = MIN(y+param_DivisionRadius, YCells-1);
  int minz = MAX(z-param_DivisionRadius,0);
  int maxz = MIN(z+param_DivisionRadius, ZCells-1);

  if(PERIODIC_X)
  {
      minx = x-param_DivisionRadius;
      maxx = x+param_DivisionRadius;
  }
  if(PERIODIC_Y)
  {
      miny = y-param_DivisionRadius;
      maxy = y+param_DivisionRadius;
  }
  if(PERIODIC_Z)
  {
      minz = z-param_DivisionRadius;
      maxz = z+param_DivisionRadius;
  }

  // To store the locations of the surrounding elements to check
  Location test[3];
  int NumNeighbours = ((maxx-minx+1)*(maxy-miny+1)*(maxz-minz+1));
  int Numtests;
  int count = 0;
  bool found = false;
  double maxoxygen = -1.0;

  // To find a random numbers, within range of the neighbouring node indices, for a random start
  int tx = minx;
  int ty = miny;
  int tz = minz;
  if (maxx>minx) tx += mtrand.randInt(maxx-minx);
  if (maxy>miny) ty += mtrand.randInt(maxy-miny);
  if (maxz>minz) tz += mtrand.randInt(maxz-minz); // ensures randInt not called for 2D

  // loop round the test elements starting at (tx,ty,tz)
  while(count < NumNeighbours)
      {
		if(PERIODIC_X) { tx=MOD(tx,(XCells));}
		if(PERIODIC_Y) { ty=MOD(ty,(YCells));}
		if(PERIODIC_Z) { tz=MOD(tz,(ZCells));}
      // if test location is a vessel, up to two alternative test
      // locations are used, hence Numtests is set to 1 or 2
      Numtests = TransVesselLocations(x, y, z, tx, ty, tz, test);

      // chose a random number for random start to loop through the test locations
      int count_test = 0;
      int t = mtrand.randInt(Numtests-1);

      // loop through the locations to test (1 or 2)
      while(count_test < Numtests)
        {
				if(PERIODIC_X) { test[t].x=MOD(test[t].x,(XCells));}
				if(PERIODIC_Y) { test[t].y=MOD(test[t].y,(YCells));}
				if(PERIODIC_Z) { test[t].z=MOD(test[t].z,(ZCells));}
          // check if test element is within grid
          if (test[t].x>=0 && test[t].x<XCells && test[t].y>=0 && test[t].y<YCells && test[t].z>=0 && test[t].z<ZCells)
            {
			  // count cells at each location
			  cellcount = StackCount(test[t].x, test[t].y, test[t].z, test[t]);
              // if test element has space
              if (cellcount < Cells[x][y][z][stack].DivCarryCap && cellcount < StackCells)
                {
                    if (!param_OxygenDependentDaughterCellPlacement)
                    {
                          NewPos = test[t];
                          found=true;
                          return found;
                    }
                    else
                    {
                        // if the level of oxygen at test element > maxoxygen
                        // Then can divide
                        if(Oxygen.GetLevelDirect(test[t].x, test[t].y, test[t].z) > maxoxygen)
                        {
                          NewPos = test[t];
                          // and update maxoxygen for subsequent comparison
                          maxoxygen = Oxygen.GetLevelDirect(test[t].x, test[t].y, test[t].z);
                          found=true;
                        }
                    }
                }
            }
           // next test location
           t++;
           if(t == Numtests) t = 0;
           // increment count_test
           count_test++;
        }
      // move to the next location to test
		if(tx == maxx) {
			tx = minx;
			if(ty == maxy) {
				ty = miny;
				if(tz == maxz)
					tz = minz;
				else
					tz++;
				}
			else
				ty++;
			}
		else {
			tx++;
			}
      count++;
    }
   return found;
}

//-------------------------------------------------------------------
//! Function to determine if the ActiveCell at (x,y,z) can move (Local stencil).
//! Returns the new position NewPos.
//!
//! This implements a (biased) random walk, e.g. depending on VEGF gradients.
//! Uses a Moore neighbourhood. Hence DCell_dt_dx2 is divided by 2.
bool CellularAutomaton::CanMove(int x, int y, int z, int stack, Location& NewPos)
{
	int cellcount;

	// Do this whether or not there is > 1 cell (all active cells can move)
	// Making sure cell can't divide further than the size of grid
	int minxc = MAX(x-1, 0);
	int maxxc = MIN(x+1, XCells-1);
	int minyc = MAX(y-1, 0);
	int maxyc = MIN(y+1, YCells-1);
	int minzc = MAX(z-1, 0);
	int maxzc = MIN(z+1, ZCells-1);

      if(PERIODIC_X)
      {
          minxc = x-param_DivisionRadius;
          maxxc = x+param_DivisionRadius;
      }
      if(PERIODIC_Y)
      {
          minyc = y-param_DivisionRadius;
          maxyc = y+param_DivisionRadius;
      }
      if(PERIODIC_Z)
      {
          minzc = z-param_DivisionRadius;
          maxzc = z+param_DivisionRadius;
      }

      int Numtests = (maxxc-minxc+1)*(maxyc-minyc+1)*(maxzc-minzc+1);

	double Pij[30];
	double myrand;
	// an empty stack-location at each test site
	Location Test[30];
	// the denominator of the movement probability
	double move_denom=Cells[x][y][z][stack].MoveCarryCap;
	double DCell_dt_dx2 = param_DCell[Cells[x][y][z][stack].Type]*param_timestep/(2*param_gridsize*param_gridsize);
	double Chemo_dx = param_VEGFChemo[Cells[x][y][z][stack].Type]/(2*param_DCell[Cells[x][y][z][stack].Type]);
	bool found = false;

	// don't allow TIPCELL movement FROM a location that has another TIPCELL already
	// this must have occurred earlier in a CellMovement loop and hence we should
	// leave the two TIPCELLS together for anastomosis
	// No need to check any other locations since we don't let the TIPCELL move whatever.
	// Hence exit this function if another TIPCELL is found and return false.
	if (Cells[x][y][z][stack].Type == TIPCELL)
	{
		for (int k=0; k<StackCells; k++)
		{
			if(Cells[x][y][z][k].Type == TIPCELL && k != stack)
				return false;
		}
	}

	// loop round the test elements assigning locations
	int locint;
	int loccount=0;
	double distance = 1.0;
	int tx,ty,tz;
	for(int txc=minxc;txc<=maxxc;txc++)
	{
		for(int tyc=minyc;tyc<=maxyc;tyc++)
		{
			for(int tzc=minzc;tzc<=maxzc;tzc++)
			{
			    tx=txc;
                ty=tyc;
                tz=tzc;

                if(PERIODIC_X) { tx=MOD(txc,(XCells));}
                if(PERIODIC_Y) { ty=MOD(tyc,(YCells));}
                if(PERIODIC_Z) { tz=MOD(tzc,(ZCells));}

				if (tx == x && ty == y && tz == z)
				{
				}
				else
				{
					// if not at current site, increment location index
					locint = loccount;
					loccount++;
					distance = sqrt((double)((txc-x)*(txc-x)+(tyc-y)*(tyc-y)+(tzc-z)*(tzc-z)));
					Test[locint].x = tx;
					Test[locint].y = ty;
					Test[locint].z = tz;
					// count cells at each location
					cellcount = StackCount(tx, ty, tz, Test[locint]);
					// numerator = empty spaces
					// +ve VEGF gradient increases movement probability
					Pij[locint] = (Cells[x][y][z][stack].MoveCarryCap - cellcount)*(1.0+Chemo_dx*(VEGF.GetLevelDirect(tx,ty,tz)-VEGF.GetLevelDirect(x,y,z)));
					// scale probability by distance to target site
					Pij[locint] = DCell_dt_dx2*Pij[locint]/(distance*distance);
					if (Cells[x][y][z][stack].Type == TIPCELL)
					{
					// sprouts can move diagonally but not if a VESSEL or PROTOVESSEL is CROSSED...
					if ((tx==x+1 && ty==y+1 && tz==z) | (tx==x+1 && ty==y-1 && tz==z) | (tx==x-1 && ty==y+1 && tz==z) | (tx==x-1 && ty==y-1 && tz==z))
						{
						if (Cells[tx][y][z][0].Type == VESSEL && Cells[x][ty][z][0].Type == VESSEL) Pij[locint] = 0.0;
						for (int k=0; k<StackCells; k++)
							{
							if((Cells[tx][y][z][k].Type == ENDOCELL) || (Cells[tx][y][z][k].Type == TIPCELL))
								{
								int p1endocellID=Cells[tx][y][z][k].VesselID;
								for (int j=0; j<StackCells; j++)
									{
									if (((Cells[x][ty][z][j].Type == ENDOCELL) || (Cells[x][ty][z][j].Type == TIPCELL)) && (Cells[x][ty][z][j].VesselID==p1endocellID)) Pij[locint] = 0.0;
									}
								}
							}
						}
                            if ((tx==x && ty==y+1 && tz==z+1) | (tx==x && ty==y+1 && tz==z-1) | (tx==x && ty==y-1 && tz==z+1) | (tx==x && ty==y-1 && tz==z-1))
                                {
                                if (Cells[x][ty][z][0].Type == VESSEL && Cells[x][y][tz][0].Type == VESSEL) Pij[locint] = 0.0;
                                for (int k=0; k<StackCells; k++)
                                    {
                                    if((Cells[x][ty][z][k].Type == ENDOCELL) || (Cells[x][ty][z][k].Type == TIPCELL))
                                        {
                                        int p1endocellID=Cells[x][ty][z][k].VesselID;
                                        for (int j=0; j<StackCells; j++)
                                            {
                                            if (((Cells[x][y][tz][j].Type == ENDOCELL) || (Cells[x][y][tz][j].Type == TIPCELL)) && (Cells[x][y][tz][j].VesselID==p1endocellID)) Pij[locint] = 0.0;
                                            }
                                        }
                                    }
                                }
                            if ((tx==x+1 && ty==y && tz==z+1) | (tx==x+1 && ty==y && tz==z-1) | (tx==x-1 && ty==y && tz==z+1) | (tx==x-1 && ty==y && tz==z-1))
                                {
                                if (Cells[tx][y][z][0].Type == VESSEL && Cells[x][y][tz][0].Type == VESSEL) Pij[locint] = 0.0;
                                for (int k=0; k<StackCells; k++)
                                    {
                                    if((Cells[tx][y][z][k].Type == ENDOCELL) || (Cells[tx][y][z][k].Type == TIPCELL))//Holger do I have to change tx to tz??
                                        {
                                        int p1endocellID=Cells[tx][y][z][k].VesselID;
                                        for (int j=0; j<StackCells; j++)
                                            {
                                            if (((Cells[x][y][tz][j].Type == ENDOCELL) || (Cells[x][y][tz][j].Type == TIPCELL)) && (Cells[x][y][tz][j].VesselID==p1endocellID)) Pij[locint] = 0.0;
                                            }
                                        }
                                    }
                                }
                            // Now test the diagonal elements
                            if ((tx==x+1 && ty==y+1 && tz==z+1)|(tx==x-1 && ty==y+1 && tz==z+1)|(tx==x+1 && ty==y-1 && tz==z+1)|(tx==x-1 && ty==y-1 && tz==z+1)|(tx==x+1 && ty==y+1 && tz==z-1)|(tx==x-1 && ty==y+1 && tz==z-1)|(tx==x+1 && ty==y-1 && tz==z-1)|(tx==x-1 && ty==y-1 && tz==z-1))
                            {
                                if ((Cells[tx][y][z][0].Type == VESSEL && Cells[x][ty][tz][0].Type == VESSEL)|(Cells[x][ty][z][0].Type == VESSEL && Cells[tx][y][tz][0].Type == VESSEL)) Pij[locint] = 0.0;
                                for (int k=0; k<StackCells; k++)
                                    {
                                    if((Cells[tx][y][z][k].Type == ENDOCELL) || (Cells[tx][y][z][k].Type == TIPCELL))
                                        {
                                        int p1endocellID=Cells[tx][y][z][k].VesselID;
                                        for (int j=0; j<StackCells; j++)
                                            {
                                            if (((Cells[x][ty][tz][j].Type == ENDOCELL) || (Cells[x][ty][tz][j].Type == TIPCELL)) && (Cells[x][ty][tz][j].VesselID==p1endocellID)) Pij[locint] = 0.0;
                                            }
                                        }
                                    }
                                for (int k=0; k<StackCells; k++)
                                    {
                                    if((Cells[x][ty][z][k].Type == ENDOCELL) || (Cells[x][ty][z][k].Type == TIPCELL))
                                        {
                                        int p1endocellID=Cells[x][ty][z][k].VesselID;
                                        for (int j=0; j<StackCells; j++)
                                            {
                                            if (((Cells[tx][y][tz][j].Type == ENDOCELL) || (Cells[tx][y][tz][j].Type == TIPCELL)) && (Cells[tx][y][tz][j].VesselID==p1endocellID)) Pij[locint] = 0.0;
                                            }
                                        }
                                    }
                            }

					// don't allow TIPCELL movement to a Node with 3 or more connections or Inflow/OutflowNode
					if (Cells[tx][ty][tz][0].NumBonds > 2) Pij[locint] = 0.0;
					// don't allow 1st step to a vessel as this could allow unwanted crossing
					// easy way is to use the fact that if the current location is a node then
					// this must be the 1st movement of the tipcell
					if (Cells[x][y][z][0].NumBonds > 0 && Cells[tx][ty][tz][0].Type == VESSEL) Pij[locint] = 0.0;
					// can't move to a location with the same protovessel?
					int isInTissue = 0;
					for (int k=0; k<StackCells; k++)
					{
						if(Cells[tx][ty][tz][k].Type == ENDOCELL && Cells[tx][ty][tz][k].VesselID == Cells[x][y][z][stack].VesselID)
							Pij[locint] = param_TipSelf*Pij[locint];
						if(Cells[tx][ty][tz][k].Type != NONE) isInTissue = 1;
					}
					if(isInTissue == 0 && param_isVesselMovementInTissueOnly) Pij[locint] = 0.0;
					// can't move if a protovessel tail is already there
					// 1) this could give a fourfold connection if the protovessel also connects later
					// 2) this could mess up when a sprout connects to another sprout
					if (taillocs[tx][ty][tz]>0) Pij[locint] = 0.0;
					// can't move to (tx,ty) if a VESSEL is there
					if (Cells[tx][ty][tz][0].Type == VESSEL) Pij[locint] = param_TipVessel*Pij[locint];
					}
					if (Pij[locint]<0)
					{
						Pij[locint] = 0.0;
					}
				}
			}
		}
	}

	//pick a random number in [0,1)
	myrand = mtrand.rand();


#if SHUFFLELOCATIONS
	shuffleLocationsAndProbabilities(Test, Pij, Numtests-1);
#endif
	// loop round the test locations: no need to check if myrand is LARGE enough,
	// since if so the function CanMove would have returned already.
	// initialise the test probability
	double PijTest = 0.0;

	for(locint=0;locint<Numtests-1;locint++)
	{
		PijTest += Pij[locint]/move_denom;
		if (myrand < PijTest)
		{
			NewPos = Test[locint];
			if (Test[locint].stack==-1)
			{
				return false;
			}
			else
				return true;
			}
		}
	return false;

}

//-------------------------------------------------------------------
//! Function to determine if the ActiveCell at (x,y,z) can move (Local stencil).
//! Returns the new position NewPos.
//!
//! This implements a pure random walk, using a Moore neighbourhood.
bool CellularAutomaton::CanMoveRandomMoore(int x, int y, int z, int stack, Location& NewPos)
{
	// Making sure cell can't divide further than the size of grid
	int minx = MAX(x-1, 0);
	int maxx = MIN(x+1, XCells-1);
	int miny = MAX(y-1, 0);
	int maxy = MIN(y+1, YCells-1);
	int minz = MAX(z-1, 0);
	int maxz = MIN(z+1, ZCells-1);

	double Pij[30];
	double myrand;
	// an empty stack-location at each test site
	Location Test[30];
	double DCell_dt_dx2 = param_DCell[Cells[x][y][z][stack].Type]*param_timestep/(2*param_gridsize*param_gridsize);
	bool found = false;

	// loop round the test elements assigning locations
	int locint;
	int loccount=0;
	double distance = 1.0;
	for(int tx=minx;tx<=maxx;tx++)
	{
		for(int ty=miny;ty<=maxy;ty++)
		{
			for(int tz=minz;tz<=maxz;tz++)
			{
				if (tx == x && ty == y && tz == z)
				{
				}
				else
				{
					// if not at current site, increment location index
					locint = loccount;
					loccount++;
					distance = sqrt((double)((tx-x)*(tx-x)+(ty-y)*(ty-y)+(tz-z)*(tz-z)));
					Test[locint].x = tx;
					Test[locint].y = ty;
					Test[locint].z = tz;
					Test[locint].stack = 0;
					Pij[locint] = DCell_dt_dx2/(distance*distance);
				}
			}
		}
	}

	//pick a random number in [0,1)
	myrand = mtrand.rand();
	// initialise the test probability
	double PijTest = 0.0;
#if SHUFFLELOCATIONS
	shuffleLocationsAndProbabilities(Test, Pij, loccount-1);
#endif
	// loop round the test locations: no need to check if myrand is LARGE enough,
	// since if so the function CanMove would have returned already.
	for(locint=0;locint<loccount;locint++)
	{
		PijTest += Pij[locint];
		if (myrand < PijTest)
		{
			NewPos = Test[locint];
			return true;
		}
	}
	return false;
}

//-------------------------------------------------------------------
//! Function to determine if the ActiveCell at (x,y,z,stack) can move (Local stencil).
//! Returns the new position NewPos.
//!
//! This implements a pure random walk, using a von Neumann neighbourhood.
bool CellularAutomaton::CanMoveRandomNeumann(int x, int y, int z, int stack, Location& NewPos)
{
	// Making sure cell can't divide further than the size of grid
	int minx = MAX(x-1, 0);
	int maxx = MIN(x+1, XCells-1);
	int miny = MAX(y-1, 0);
	int maxy = MIN(y+1, YCells-1);
	int minz = MAX(z-1, 0);
	int maxz = MIN(z+1, ZCells-1);

	double Pij[6]; // only 6 target locations
	double myrand;
	// an empty stack-location at each test site
	Location Test[6];
	double DCell_dt_dx2 = param_DCell[Cells[x][y][z][stack].Type]*param_timestep/(param_gridsize*param_gridsize);
	bool found = false;

	// loop round the test elements assigning locations
	int locint;
	int loccount=0;
	double distance;
	for(int tx=minx;tx<=maxx;tx++)
	{
		for(int ty=miny;ty<=maxy;ty++)
		{
			for(int tz=minz;tz<=maxz;tz++)
			{
				if (tx-x == ty-y || tx-x == y-ty || tx-x == tz-z || tx-x == z-tz || ty-y == tz-z || ty-y == z-tz)
				{
				}
				else
				{
					// if not at current site or a diagonal move, increment location index
					locint = loccount;
					loccount++;
					Test[locint].x = tx;
					Test[locint].y = ty;
					Test[locint].z = tz;
					Test[locint].stack = 0;
					Pij[locint] = DCell_dt_dx2;
				}
			}
		}
	}

	//pick a random number in [0,1)
	myrand = mtrand.rand();
	// initialise the test probability
	double PijTest = 0.0;
#if SHUFFLELOCATIONS
	shuffleLocationsAndProbabilities(Test, Pij, loccount-1);
#endif
	// loop round the test locations: no need to check if myrand is LARGE enough,
	// since if so the function CanMove would have returned already.
	for(locint=0;locint<loccount;locint++)
	{
		PijTest += Pij[locint];
		if (myrand < PijTest)
		{
			NewPos = Test[locint];
			return true;
		}
	}
	return false;
}

//-------------------------------------------------------------------
//! Cell cycle model test
void CellularAutomaton::CellCycleTest(int nsteps, string outstem)
{
	int i = 1;
	int j;
	bool appendsubcellular = true;
	char filename[200];
	sprintf(filename, "./%s/subcellular-temporal.sub", outstem.c_str());
	SubCellularToFile(filename, appendsubcellular);
	while(i<=param_nsumwrites)
	{
		for(j=1;j<=nsteps;j++)
        {
			RunCellCycleStep();
			GLOBALTIMESTEP+=(int)param_timestep;
			sprintf(filename, "./%s/subcellular-temporal.sub", outstem.c_str());
			SubCellularToFile(filename, appendsubcellular);
		}
		printf("***** End of step %d *****\n",i);
		i++;
	}
	printf("End of CellCycleTest\n");
}

//-------------------------------------------------------------------
void CellularAutomaton::CellDeath()
{
	Location pos;
	int count=0, NumActiveCells=ActiveCells.size();

	for(int i=0;i<ActiveCells.size();i++)
	{
		pos = ActiveCells[i];
		if(Cells[pos.x][pos.y][pos.z][pos.stack].InternalCellDeath())
		{
			if (VERBOSE) printf("Celldeath = true at (%d,%d,%d,%d), p53=%f\n",pos.x,pos.y,pos.z,pos.stack,Cells[pos.x][pos.y][pos.z][pos.stack].p53VEGF.p53);
			KillCell(pos.x, pos.y, pos.z, pos.stack);
			aposum[pos.x][pos.y][pos.z]++;
			// Killing a cell decreases the number of members
			// decrement i by one to check the member that would otherwise be skipped
			i--;
		}
		count++;
	}
}

//-------------------------------------------------------------------
//! @todo Implement resetPhaseAfterFailedDivision as a parameter?
void CellularAutomaton::CellDivision()
{
	clock_t mystart = clock();
	Location pos;
	int i, imax=ActiveCells.size(), count=0, killcount=0, divcount=0;
	int NumActiveCells=ActiveCells.size();

	if (VERBOSE) printf("Starting CellularAutomaton::CellDivision\n");

//	for(i=0;i<ActiveCells.size();i++)
	for(i=0;i<imax;i++)
    {
		pos = ActiveCells[i];
		if(Cells[pos.x][pos.y][pos.z][pos.stack].Type == NONE)
        {
			printf("Checking for division of type NONE, (%d,%d,%d,%d)\n",pos.x,pos.y,pos.z,pos.stack);
			exit(1);
		}
		if(Cells[pos.x][pos.y][pos.z][pos.stack].InternalCellDivision())
		{
			//Internal factors are right for division - what about external ones?
			Location NewPos;
			if(CanDivide(pos.x, pos.y, pos.z, pos.stack, NewPos))
				{
				if (Cells[pos.x][pos.y][pos.z][pos.stack].Intercalated == false)
					{
					//NewPos now contains a location for the cell to divide into
					MakeCellDaughterFromParent(pos,NewPos);
					//Reset current cell
					Cells[pos.x][pos.y][pos.z][pos.stack].DivisionReset();
					//this is a local counter to check how the number of active cells changes with divisions
					divcount++;
					// increment the tally of total divisions at this CA location
					divsum[pos.x][pos.y][pos.z]++;
				}
				else
				{
					if (VERBOSE) printf("Killing dividing cell due to drug intercalation at (%d,%d,%d,%d)\n",pos.x,pos.y,pos.z,pos.stack);
					KillCell(pos.x, pos.y,pos.z,pos.stack);
					killsum[pos.x][pos.y][pos.z]++;
					//this is a local counter to check how the number of active cells changes with killing
					killcount++;
					// Killing a cell decreases the number of members
					// decrement i by one to check the member that would otherwise be skipped
					i--;
					imax--;
				}
			}
			else
			{
				if (VERBOSE) printf("Failed cell division at (%d,%d,%d,%d) due to lack of space\n",pos.x,pos.y,pos.z,pos.stack);
				int resetPhaseAfterFailedDivision=1;
				if (resetPhaseAfterFailedDivision)
				{
					//Reset current cell
					Cells[pos.x][pos.y][pos.z][pos.stack].DivisionReset();
				}
				else
				{
					KillCell(pos.x, pos.y, pos.z, pos.stack);
					diesum[pos.x][pos.y][pos.z]++;
					//this is a local counter to check how the number of active cells changes with killing
					killcount++;
					// Killing a cell decreases the number of members
					// decrement i by one to check the member that would otherwise be skipped
					i--;
					// maybe also decrement imax, and don't increase imax when a cell divides
					// as there should be no way a fresh daughter could divide anyway
					imax--;
				}
					if(Cells[pos.x][pos.y][pos.z][pos.stack].Type == CANCER && Cells[pos.x][pos.y][pos.z][0].Type == VESSEL)
						Vasc.CompressVesselDueToProliferation(pos.x, pos.y, pos.z);
            }
        }
		count++;
    }

	// print a warning message if we haven't checked all initially active cells
	if (count != NumActiveCells)
		printf("WARNING: CellDivision: Cells at start=%d, end=%zu, checked=%d, killed=%d, divided=%d\n",NumActiveCells,ActiveCells.size(),count,killcount,divcount);

	clock_t mystop = clock();
	if (VERBOSE) printf("FINISHED CellularAutomaton::CellDivision, elapsed time = %g seconds\n",(mystop-mystart)/(double)CLOCKS_PER_SEC);
}

//-------------------------------------------------------------------
//! Function to calculate movement of all Active Cells.
//! Calls CanMove for each Active Cell and if this returns TRUE,
//! call MoveCell to move the cell to NewPos.
void CellularAutomaton::CellMovement()
{
	clock_t mystart = clock();
	if (VERBOSE) printf("Starting CellularAutomaton::CellMovement\n");

	Location pos;
	Location NewPos;
	for(int i=0;i<ActiveCells.size();i++)
		{
		pos = ActiveCells[i];

		if (Cells[pos.x][pos.y][pos.z][pos.stack].isMoving)
			{

			if(Cells[pos.x][pos.y][pos.z][pos.stack].Type != VESSEL)
				{
				Cells[pos.x][pos.y][pos.z][pos.stack].MovedThisStep = false;

				if(CanMove(pos.x,pos.y,pos.z,pos.stack,NewPos))
					{
					if(Cells[pos.x][pos.y][pos.z][pos.stack].Type == NONE)
						printf("Trying to move cell type NONE\n");
					Cells[pos.x][pos.y][pos.z][pos.stack].MovedThisStep = true;
					MoveCell(pos,NewPos);
					ActiveCells[i]=NewPos;
					if(Cells[NewPos.x][NewPos.y][NewPos.z][NewPos.stack].Type == NONE)
						printf("Moved cell changed to NONE\n");
					}
				else if(Cells[pos.x][pos.y][pos.z][pos.stack].Type == NONE)
					printf("Non-moved cell changed to NONE\n");
				}
			}
		}

	clock_t mystop = clock();
	TimeCellMovement += (mystop-mystart)/(double)CLOCKS_PER_SEC;
	if (VERBOSE) printf("@@@@@@ TIMING @@@@@@: Time spent in CellMovement = %g seconds\n",TimeCellMovement);
}
//-------------------------------------------------------------------
void CellularAutomaton::CellQuiescence()
{
	Location pos;

	for(int i=0;i<ActiveCells.size();i++)
	{
		pos = ActiveCells[i];
		Cells[pos.x][pos.y][pos.z][pos.stack].UpdateQuiescence();
	}
}

//-------------------------------------------------------------------
void CellularAutomaton::CellTypeToFile(char* filename)
{
  FILE* fp = fopen(filename, "w");
  fprintf(fp, "%d, %d, %d, %d\n", XCells, YCells, ZCells, StackCells);

	for(int x=0;x<XCells;x++)
	{
		for(int y=0;y<YCells;y++)
		{
			for(int z=0;z<ZCells;z++)
			{
				for(int stack=0;stack<StackCells;stack++)
				{
					fprintf(fp, "%d", Cells[x][y][z][stack].CellTypeToFile() );
					if(stack!=StackCells-1) fprintf(fp, ", ");
				}
			fprintf(fp, "\n");
			}
		}
	}

  fclose(fp);
}

//-------------------------------------------------------------------
//! Check if the ActiveCells list contains any active cells.
//! If so, exit with error message.
void CellularAutomaton::CheckForNONEInActiveCells()
{
	Location pos;
	for(int i=0;i<ActiveCells.size();i++)
	{
		pos = ActiveCells[i];

		if(Cells[pos.x][pos.y][pos.z][pos.stack].Type==NONE)
		{
			printf("ActiveCell of Type NONE, (%d,%d,%d,%d)\n",pos.x,pos.y,pos.z,pos.stack);
			exit(1);
		}
	}
}

//-------------------------------------------------------------------
//! @brief Create Output Directories.
//!
//! @param outstem Name of the directory in which all outputs are placed.
//! @param parstem The parameter file name, which is used as the directory
//! name to place output for a simulation.
void CellularAutomaton::CreateOutputDirectories(string outstem, string parstem)
{
	std::vector<string> subdirs;
	subdirs.push_back("");
	subdirs.push_back("/subcellular");
	subdirs.push_back("/cells");
	subdirs.push_back("/cells/type");
	subdirs.push_back("/sum");
	subdirs.push_back("/randstate");
	if (param_WritePovRay) subdirs.push_back("/PovRay");
	if (param_WriteOpenInventor) subdirs.push_back("/OpenInventor");

	if (param_isOxygenDiffusible) subdirs.push_back("/oxygen");
	if (param_isVEGFDiffusible) subdirs.push_back("/VEGF");
	if (param_isDrugDiffusible) subdirs.push_back("/drug");
	if (param_isProDrugDiffusible) subdirs.push_back("/prodrug");

	if (param_Vasculature) {
		subdirs.push_back("/Vasculature");
		subdirs.push_back("/Segments");
	}
	CreateOutputDirectoriesFromList(outstem, parstem, subdirs);
}

//-------------------------------------------------------------------
//! @brief Create Output Directories from a vector of directory names
//!
//! @param outstem Name of the directory in which all outputs are placed.
//! @param parstem The parameter file name, which is used as the directory
//! name to place output for a simulation.
//! @todo Eliminate dircount and use an iterator over subdirs instead
void CellularAutomaton::CreateOutputDirectoriesFromList(string outstem, string parstem, std::vector<string> subdirs)
{
	int error, i;
	int dircount = subdirs.size();

	string subdir;

	cout << "************************************************************\n";
	cout << "*** Checking output directory structure.\n";
	cout << "************************************************************\n";

	// try to create directories for output
	for(i=0;i<dircount+1;i++)
 	{
 		if(i==0)
        {
 		// first make sure "output" is there.
	 		subdir=outstem;
			outstem+=parstem;
		}
 		else
 		{
 			subdir=outstem+subdirs[i-1];
 		}
 		// try to make the directory
		error=mkdir(subdir.c_str(),S_IRWXU);
		if(error!=0)
		{
			if (errno==EEXIST)
			{
				// if the directory already exists, check that it has appropriate permissions
				printf("WARNING: Output directory %s already exists\n",subdir.c_str());
				// this structure stores file/directory information
				struct stat buf;
				// get information on the desired output directory
				stat(subdir.c_str(),&buf);
				// if it's a directory, OK, otherwise quit
				if ((buf.st_mode & S_IFMT) == S_IFDIR)
					cout<<subdir<<" is a directory";
				else
				{
					cout<<subdir<<" is not a directory, quitting...\n";
					exit(1);
				}

				if ((buf.st_mode & S_IRWXU) == S_IRWXU)
					cout<<" and is Read/Write/Executable by User\n";
				else
				{
					cout<<" BUT is NOT Read/Write/Executable by User - changing with CHMOD...\n";
					chmod(subdir.c_str(),S_IRWXU);
				}
			}
			else
			{
				printf("ERROR: failed to create directory %s, %i, %s\n",subdir.c_str(),errno,strerror(errno));
				exit(1);
			}
		}
		else
		cout<<"Successfully created output directory "<<subdir.c_str()<<"\n";
	}
}

//-------------------------------------------------------------------
//! Print a textual representation of the CellularAutomaton to stdout
//! @todo Needs to account properly for the stack
void CellularAutomaton::DEBUG_Print()
{
	int i,j,k;

	printf("\n");
	for(k=0;k<ZCells;k++)
	{
	for(i=0;i<XCells;i++)
	{
		for(j=0;j<YCells;j++)
		{
		    if((Cells[i][j][k][0].Type == TIPCELL) || (Cells[i][j][k][1].Type == TIPCELL) || (Cells[i][j][k][2].Type == TIPCELL) || (Cells[i][j][k][3].Type == TIPCELL)) putchar('T');
		    else if((Cells[i][j][k][0].Type == ENDOCELL) || (Cells[i][j][k][1].Type == ENDOCELL) || (Cells[i][j][k][2].Type == ENDOCELL) || (Cells[i][j][k][3].Type == ENDOCELL)) putchar('E');
			//	if(Cells[i][j][k][0].Type == NONE) putchar('.');
				else if(Cells[i][j][k][0].Type == NORMALCELL) putchar('o');
				else if(Cells[i][j][k][0].Type == CANCER && !Cells[i][j][k][0].Quiescent) putchar('X');
				else if(Cells[i][j][k][0].Type == CANCER && Cells[i][j][k][0].Quiescent) putchar('z');
				else if(Cells[i][j][k][0].Type == VESSEL) putchar('*');
				else if(Cells[i][j][k][0].Type == MACROPHAGE) putchar('Z');
				else if(Cells[i][j][k][0].Type == MONOCYTE) putchar('Y');
				else printf("%d", Cells[i][j][k][0].Type);
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
	}
}

//-------------------------------------------------------------------
//! Delete memory allocated to a CellularAutomaton
void CellularAutomaton::DeleteMemory()
{
	if(Cells)
	{
		for(int i=0;i<XCells;i++)
		{
			for(int j=0;j<YCells;j++)
			{
				for (int k=0; k<ZCells; k++)
				{
					delete [] Cells[i][j][k];
				}
				delete [] Cells[i][j];
				delete [] divsum[i][j];
				delete [] macsum[i][j];
				delete [] monosum[i][j];
				delete [] aposum[i][j];
				delete [] diesum[i][j];
				delete [] killsum[i][j];
				delete [] taillocs[i][j];
				delete [] tiplocs[i][j];
				delete [] tipcelloutsum[i][j];
				delete [] drugsum[i][j];
			}
			delete [] Cells[i];
			delete [] divsum[i];
			delete [] macsum[i];
			delete [] monosum[i];
			delete [] aposum[i];
			delete [] diesum[i];
			delete [] killsum[i];
			delete [] taillocs[i];
			delete [] tiplocs[i];
			delete [] tipcelloutsum[i];
			delete [] drugsum[i];
		}
		delete [] Cells;
		delete [] divsum;
		delete [] macsum;
		delete [] monosum;
		delete [] aposum;
		delete [] diesum;
		delete [] killsum;
		delete [] taillocs;
		delete [] tiplocs;
		delete [] tipcelloutsum;
		delete [] drugsum;
	}
}

//-------------------------------------------------------------------
//! Function to empty all Cells in the CA.
//! Calls KillCell at all locations.
//! Need to be sure that XCells, YCells and ZCells are initialised.
void CellularAutomaton::EmptyCACells()
{
	int x, y, z, stack;

	for(x=0;x<XCells;x++)
	{
		for(y=0;y<YCells;y++)
		{
			for(z=0;z<ZCells;z++)
			{
				for(stack=0;stack<StackCells;stack++)
				{
					KillCell(x,y,z,stack);
				}
			}
		}
	}
}

//-------------------------------------------------------------------
//! Calculate the extravasation condition for endothelial tip cells
//!
double CellularAutomaton::ExtravasationConditionTIPCELL(int v, int curx, int cury, int curz)
{
	double extravasation_condition;

	if (Cells[curx][cury][curz][0].NumBonds < 3)
	// to exclude sprouts from 3-nodes
	{
//			extravasation_condition = param_timestep*0.5*param_SproutMax*twopi*Vasc.Vessels[v].Radius*Vasc.Vessels[v].Length*VEGF.GetLevelDirect(curx,cury,curz)/(param_InfluxThreshold[TIPCELL]+VEGF.GetLevelDirect(curx,cury,curz));
		extravasation_condition = param_timestep*0.5*param_SproutMax*VEGF.GetLevelDirect(curx,cury,curz)/(param_InfluxThreshold[TIPCELL]+VEGF.GetLevelDirect(curx,cury,curz));
		//printf("VEGF.GetDirect = %.3f\n",VEGF.GetLevelDirect(curx,cury,curz));
		if(Vasc.Vessels[v].State == COLLAPSEVESSEL)
		{
			extravasation_condition = 0.;
			if (VERBOSE) printf("Trying extravasation from coopted vessel %d. VEGF=%f\n",v,Vasc.Vessels[v].MeanVEGF);
		}

		int minxc = MAX(curx-param_TipSproutRadius, 0);
		int maxxc = MIN(curx+param_TipSproutRadius, XCells-1);
		int minyc = MAX(cury-param_TipSproutRadius, 0);
		int maxyc = MIN(cury+param_TipSproutRadius, YCells-1);
		int minzc = MAX(curz-param_TipSproutRadius, 0);
		int maxzc = MIN(curz+param_TipSproutRadius, ZCells-1);

		int minx = minxc;
		int maxx = maxxc;
		int miny = minyc;
		int maxy = maxyc;
		int minz = minzc;
		int maxz = maxzc;

		if(PERIODIC_X)
		{
			minxc = curx-param_TipSproutRadius;
			minx = MOD(minxc,(XCells));
			maxxc = curx+param_TipSproutRadius;
			maxx = MOD(maxxc,(XCells));
		}
		if(PERIODIC_Y)
		{
			minyc = cury-param_TipSproutRadius;
			miny = MOD(minyc,(YCells));
			maxyc = cury+param_TipSproutRadius;
			maxy = MOD(maxyc,(YCells));
		}
		if(PERIODIC_Z)
		{
			minzc = curz-param_TipSproutRadius;
			minz = MOD(minzc,(ZCells));
			maxzc = curz+param_TipSproutRadius;
			maxz = MOD(maxzc,(ZCells));
		}

		// loop round the neighbourhood checking for existing sprouts
		int tailcount=0;
		int tx,ty,tz;
		for(int txc=minxc;txc<=maxxc;txc++)
		{
			for(int tyc=minyc;tyc<=maxyc;tyc++)
			{
				for(int tzc=minzc;tzc<=maxzc;tzc++)
				{
					tx=txc;
					ty=tyc;
					tz=tzc;

					if(PERIODIC_X) { tx=MOD(txc,(XCells));}
					if(PERIODIC_Y) { ty=MOD(tyc,(YCells));}
					if(PERIODIC_Z) { tz=MOD(tzc,(ZCells));}
					tailcount += taillocs[tx][ty][tz];
				}
			}
		}
		if (tailcount > 0) extravasation_condition = 0;
	}
	else
		extravasation_condition = 0;

	return extravasation_condition;
}

//-------------------------------------------------------------------
//! Checks a single vessel, vessel (given as an index into the Vessels array),
//! for extravasation to the CellularAutomaton CA.
//! This function writes over the contents of cells in the grid.
//!
//! Considers the Start and End Nodes of Vessels[v] and includes the
//! random extravasation of macrophages at each CA location.
//!
//! Parameters are InfluxMax (superseded), InfluxThreshold, ExtravasateCarryCap
void CellularAutomaton::ExtravasationFromVesselSite(int v, int Type, Node nodePos)
{
	int curx = nodePos.x;
	int cury = nodePos.y;
	int curz = nodePos.z;
	Location pos = Location(curx, cury, curz, -1);
	double myrand;

	// calculate probability of extravasation
	double extravasation_condition = 0.0;
	if (Type==MONOCYTE)
		{
		// for macrophage probability should increase with VEGF concentration (and linearly with haematocrit?)
//			extravasation_condition = param_InfluxMax[Type]*twopi*Vessels[v].Radius*Vessels[v].Haematocrit/(param_InfluxThreshold[Type]+twopi*Vessels[v].Radius*Vessels[v].Haematocrit);
		extravasation_condition = Vasc.MonoBloodConcentration*twopi*Vasc.Vessels[v].Radius*Vasc.Vessels[v].Haematocrit*VEGF.GetLevelDirect(curx,cury,curz)/(param_InfluxThreshold[Type]+VEGF.GetLevelDirect(curx,cury,curz));
		}
	else if (Type==MACROPHAGE)
		{
		extravasation_condition = Vasc.MacExtravasationCondition(VEGF, v, MACROPHAGE, curx, cury, curz);
		if (FLOWVERBOSE) printf("MAC extravasation_condition:%f \n",extravasation_condition);
		}
	else if (Type==TIPCELL)
		extravasation_condition = ExtravasationConditionTIPCELL(v, curx, cury, curz);

	// check if there is space
	int filled = StackCount(curx,cury,curz,pos);
	int newmacloc = pos.stack;

	//pick a random number in [0,1)
	myrand = mtrand.rand();
	// if there is space and the die says so, extravasate (make a new cell)
	if ( filled < param_ExtravasateCarryCap[Type] && myrand < extravasation_condition && newmacloc != -1)
		{
		MakeCell(curx,cury,curz,newmacloc,Type);
		pos.stack = newmacloc;
		if (VERBOSE) printf("Extravasation, cell type %d at (%d,%d,%d,%d); myrand = %f, H = %f, cond = %f\n",
			Type, curx, cury, curz, newmacloc, myrand, Vasc.Vessels[v].Haematocrit, extravasation_condition);
		if (Type==MACROPHAGE) macsum[curx][cury][curz]++;
		if (Type==MONOCYTE) monosum[curx][cury][curz]++;
		if (Type==TIPCELL)
			{
			int PVID;
			if (Vasc.AvailableProtoVessels.size() > 0)
			{
				PVID = Vasc.AvailableProtoVessels[0];
				Vasc.AvailableProtoVessels.erase(Vasc.AvailableProtoVessels.begin());
				// remove the previous list entries
				Vasc.ProtoVessels[PVID].xyzstacklist.clear();
			}
			else
			{
				PVID = Vasc.NumProtoVessels++;
				ProtoVessel tempProtoVessel;
				Vasc.ProtoVessels.push_back(tempProtoVessel);
			}

			printf("TIPCELL %d extravasated at (x,y,z,stack)=(%d,%d,%d,%d)\n",PVID,curx,cury,curz,newmacloc);
			Cells[curx][cury][curz][newmacloc].VesselID = PVID;
			Vasc.ProtoVessels[PVID].Connected = false;
			Vasc.ProtoVessels[PVID].Age = 0.0;
			Vasc.ProtoVessels[PVID].xyzstacklist.push_back(pos);

				////////////////////////////////////////
				// The extravasation site should already have a node at (curx,cury)
				////////////////////////////////////////
				// the index of the new node
				if (Vasc.Nodes[Vasc.Vessels[v].Start].x == curx && Vasc.Nodes[Vasc.Vessels[v].Start].y == cury && Vasc.Nodes[Vasc.Vessels[v].Start].z == curz)
					Vasc.ProtoVessels[PVID].TailNode = Vasc.Vessels[v].Start;
				else
					Vasc.ProtoVessels[PVID].TailNode = Vasc.Vessels[v].End;

			// note this means the real ProtoVessel index is one less than taillocs/tiplocs
			taillocs[curx][cury][curz]=PVID+1;
			tiplocs[curx][cury][curz]=PVID+1;
			tipcelloutsum[curx][cury][curz]++;
			}
		}

}

//-------------------------------------------------------------------
//! Cell extravasation to the CellularAutomaton CA.  Calls ExtravasationToCA
//! for every vessel in the Vasculature.  Note that this function writes the
//! vasculature over the current Cells array in the CA, meaning that any
//! elements of Cells containing part of the vasculature have their contents
//! overwritten.
//!
//! Checks a single vessel, vessel (given as an index into the Vessels array),
//! for extravasation to the CellularAutomaton CA. Uses the Start and End Nodes,
//! adding cells of type Type to CA's Cells grid.
//! This function writes over the contents of cells in the grid.
//!
//! Considers the Start and End Nodes of Vessels[v] and includes the
//! random extravasation of macrophages at each CA location.
//!
//! Parameters are InfluxMax (superseded), InfluxThreshold, ExtravasateCarryCap
//! @todo Update this documentation
//! @todo Call with NONE is needed to reproduce Cancer Res results due to extra random number calls
void CellularAutomaton::ExtravasationToCA()
{
	int i,k;
	int curx, cury, curz;
	Location pos;

	for(k=0;k<Vasc.ActiveVessels.size();k++)
	{
		i = Vasc.ActiveVessels[k];
		Node startnode = Vasc.Nodes[Vasc.Vessels[i].Start];
		Node endnode = Vasc.Nodes[Vasc.Vessels[i].End];
		ExtravasationFromVesselSite(i, MACROPHAGE, startnode);
		ExtravasationFromVesselSite(i, MACROPHAGE, endnode);

		ExtravasationFromVesselSite(i, MONOCYTE, startnode);
		ExtravasationFromVesselSite(i, MONOCYTE, endnode);

		ExtravasationFromVesselSite(i, NONE, startnode);
		ExtravasationFromVesselSite(i, NONE, endnode);

		ExtravasationFromVesselSite(i, TIPCELL, startnode);
		ExtravasationFromVesselSite(i, TIPCELL, endnode);
	}
}

//-------------------------------------------------------------------
//! Initialise a CellularAutomaton
//! Allocates space for various arrays and initialises all Diffusibles
void CellularAutomaton::Initialise(int NumXCells, int NumYCells, int NumZCells, int NumStackCells)
{
	int i,j,k,l;
	XCells = NumXCells;
	YCells = NumYCells;
	ZCells = NumZCells;
	StackCells = NumStackCells;

	MaxCellID = 0;

	Cells = new CAElement***[XCells];
	divsum = new int**[XCells];
	macsum = new int**[XCells];
	monosum = new int**[XCells];
	aposum = new int**[XCells];
	diesum = new int**[XCells];
	killsum = new int**[XCells];
	taillocs = new int**[XCells];
	tipcelloutsum = new int**[XCells];
	tiplocs = new int**[XCells];
	drugsum = new double**[XCells];
	TACEParticlessum = new double**[XCells];

  	for(i=0;i<XCells;i++)
	{
		Cells[i] = new CAElement**[YCells];
		divsum[i] = new int*[YCells];
		macsum[i] = new int*[YCells];
		monosum[i] = new int*[YCells];
		aposum[i] = new int*[YCells];
		diesum[i] = new int*[YCells];
		killsum[i] = new int*[YCells];
		taillocs[i] = new int*[YCells];
		tipcelloutsum[i] = new int*[YCells];
		tiplocs[i] = new int*[YCells];
		drugsum[i] = new double*[YCells];
		TACEParticlessum[i] = new double*[YCells];
	    for(j=0;j<YCells;j++)
	    {
			Cells[i][j] = new CAElement*[ZCells];
			divsum[i][j] = new int[ZCells];
			macsum[i][j] = new int[ZCells];
			monosum[i][j] = new int[ZCells];
			aposum[i][j] = new int[ZCells];
			diesum[i][j] = new int[ZCells];
			killsum[i][j] = new int[ZCells];
			taillocs[i][j] = new int[ZCells];
			tiplocs[i][j] = new int[ZCells];
			tipcelloutsum[i][j] = new int[ZCells];
			drugsum[i][j] = new double[ZCells];
			TACEParticlessum[i][j] = new double[ZCells];
			for (k=0; k<ZCells; k++)
            {
				Cells[i][j][k] = new CAElement[StackCells];
            }
	    }
	}

	for(i=0;i<XCells;i++)
	{
		for(j=0;j<YCells;j++)
		{
			for(k=0;k<ZCells;k++)
			{
				divsum[i][j][k] = 0;
				macsum[i][j][k] = 0;
				monosum[i][j][k] = 0;
				aposum[i][j][k] = 0;
				diesum[i][j][k] = 0;
				killsum[i][j][k] = 0;
				// no endothelial sprouts initially
				taillocs[i][j][k] = 0;
				tipcelloutsum[i][j][k] = 0;
				// no endothelial tips initially
				tiplocs[i][j][k] = 0;
				drugsum[i][j][k] = 0.0;
				TACEParticlessum[i][j][k] = 0.0;
				for (l=0;l<StackCells;l++)
				{
					Cells[i][j][k][l].InitialiseElement();
				}
			}
		}
	}

	Oxygen.Initialise(NumXCells, NumYCells, NumZCells);
	VEGF.Initialise(NumXCells, NumYCells, NumZCells);
	Drug.Initialise(NumXCells, NumYCells, NumZCells);
	ProDrug.Initialise(NumXCells, NumYCells, NumZCells);

	VascPropSource.clear();
	for(int i=0;i<NumXCells*NumYCells*NumZCells;i++) VascPropSource.push_back(0.0);
}

//-------------------------------------------------------------------
//! Initialise simulation
void CellularAutomaton::InitialiseSimulation()
{
	string InitialConfigFile;

	// find the suffix of InitialConfigFile
	InitialConfigFile = param_InitialConfigFile;
	// find the last occurence of '.' in InitialConfigFile
	int i = InitialConfigFile.rfind('.', InitialConfigFile.length());
	// if there is a '.', find the suffix, otherwise quit
	string InitialConfigFile_suffix;
	if (i != string::npos) {
		InitialConfigFile_suffix = InitialConfigFile.substr(i+1,InitialConfigFile.length()-i);
		}
	else
		{
		printf("Cannot extract InitialConfigFile_suffix ... will not add any cells.\n");
		}

	// initialise GLOBALTIMESTEP
	GLOBALTIMESTEP = 0;
    printf("%s",InitialConfigFile_suffix.c_str());
    printf("%s",InitialConfigFile_suffix.c_str());

	if (InitialConfigFile_suffix == "ini")
		{
		// initialise random generator with integer seed
		mtrand.seed(param_RandSeed);
		printf(" MersenneRandom initialised with %d, Random = %lu\n", param_RandSeed,mtrand.randInt());
		}
	else if (InitialConfigFile_suffix == "sub")
		{
		SubCellularFromFile(InitialConfigFile);
		// want to replace "subcellular" with "randstate" and ".sub" with ".rnd"
		// to load the state of the random number generator
		// See the example in: http://www.cplusplus.com/reference/string/string/rfind.html
		string InitialRandFile=InitialConfigFile;
		string key ("subcellular");
		i = InitialRandFile.rfind(key);
		if (i != string::npos)
			InitialRandFile.replace(i,key.length(),"randstate");
		key = "sub";
		i = InitialRandFile.rfind(key);
		 printf("InitialRandFile '%s' ... quitting \n", InitialRandFile.c_str());
		if (i != string::npos)
			InitialRandFile.replace(i,key.length(),"rnd");
			printf("InitialRandFile '%s' ... quitting \n", InitialRandFile.c_str());
		string folderInitialRandFile = "input/" + InitialRandFile;
		ifstream stateIn(folderInitialRandFile.c_str());
	//	if( stateIn )
	//		{
	//		     printf("We are here.\n");
	//		stateIn >> mtrand;
	//		stateIn.close();
	//		 printf("We are here.\n");
	//		}
	//	else
	//		{
	//		    printf("FolderInitialRandFile '%s' ... quitting \n", folderInitialRandFile.c_str());
	//		printf("You probably have to rename frameXXXX.rnd -> subframeXXXX.sub \n",InitialRandFile.c_str());
	//		printf("Cannot open randstate stream '%s' ... quitting \n",InitialRandFile.c_str());
	//		exit(1);
	//		}
		}
	else
		{
		printf("InitialConfigFile_suffix is not 'ini' or 'sub' ... will not add any cells.\n");
		mtrand.seed(param_RandSeed);
		}

	if(param_Vasculature)
	{
	    printf("We are here.\n");
		Vasc.ReadVessels(param_InitialVesselFile);
		Vasc.RefineVessels();
		Vasc.CheckNodeConnectivity(taillocs);
		AddVesselsToCA();
		Vasc.CalculateCooptionList();
		if (InitialConfigFile_suffix == "ini")
			UpdateVasculature(); AddVesselsToCA();
	}

	if (InitialConfigFile_suffix == "ini")
	{
		InitialConfigFile = param_InitialConfigFile;
		LoadConfig(InitialConfigFile);
	}

	UpdateDiffusibles();

	if (param_WoundFile!="NONE")
	{
		Wound(param_WoundFile);
		if (param_Vasculature) {
			Vasc.CheckNodeConnectivity(taillocs);
			AddVesselsToCA();
			Vasc.CalculateCooptionList();
			if (InitialConfigFile_suffix == "ini") UpdateVasculature(); AddVesselsToCA();
		}
	}
}

//-------------------------------------------------------------------
//! Kill the Cell at a given location.
void CellularAutomaton::KillCell(int x, int y, int z, int stack)
{
	ActiveCells.erase(std::remove(ActiveCells.begin(), ActiveCells.end(), Location(x,y,z,stack)), ActiveCells.end());
	Cells[x][y][z][stack].KillElement();
}

//-------------------------------------------------------------------
//! Load initial CA configuration from file.
//! Initial CA configuration is specified by an input file
//! defining blocks of cells of various types.
//! @param filename a constant character argument
//!
//! @todo Include check that cell blocks lie in the simulation domain (else exit)
//! @todo Add capability to define circular initial blocks (copy from Macrophage branch)
void CellularAutomaton::LoadConfig(string filename)
{
	FILE* fp = fopen(filename.c_str(), "rt");

	if(!fp)
	// input file is bad so print error, try looking in input/, otherwise exit
	{
		printf("WARNING: Error opening initial cell config file \"%s\" ...\n", filename.c_str());
		filename = "input/" + filename;
		printf("   ... trying \"%s\".\n",filename.c_str());
		fp = fopen(filename.c_str(), "rt");
		if(!fp)
		{
			printf("ERROR: Error opening config file \"%s\" ...\n",filename.c_str());
			exit(1);
		}
	}

	char line[500];
	char type[50];
	char nums[100];
	int startx, endx, starty, endy, startz, endz;

	cout << "************************************************************\n";
	cout << "*** Reading initial config file \"" << filename << "\"\n";
	cout << "************************************************************\n";

	while( fgets(&line[0], 500, fp) )
	{
		startx=-1; endx=-1; starty=-1; endy=-1; startz=-1; endz=-1;

		int i=0;
		char c=line[0];

		// c!=0 checks if c is the null character. This could also be written c!='\0'.
		while(c!=0 && c!='\n' && c!=':')
		{
			type[i++] = c;
			c = line[i];
		}
		type[i]=0;
		if(VERBOSE) printf("CELL BLOCK:   TYPE: %s  ", type);

		while(c!='(' && c!=0 && c!='\n') c=line[i++];

		int j=0;
		c=line[i++];
		while(c!=')' && c!=0 && c!='\n')
		{
			if(c!=' ') nums[j++] = c;
			c=line[i++];
		}
		nums[j]=0;
		sscanf(&nums[0], "%d,%d,%d", &startx, &starty, &startz);
		if(VERBOSE) printf("START:  %d, %d, %d ", startx, starty, startz);

		while(c!='(' && c!=0 && c!='\n') c=line[i++];

		j=0;
		c=line[i++];
		while(c!=')' && c!=0 && c!='\n')
		{
			if(c!=' ') nums[j++] = c;
			c=line[i++];
		}
		nums[j]=0;
		sscanf(&nums[0], "%d,%d,%d", &endx, &endy, &endz);
		if(VERBOSE) printf("END: %d, %d, %d \n", endx, endy, endz);

		if(startx!=-1 && starty!=-1 && startz!=-1)
		{
			int celltype;
			if(!strcmp(&type[0], "cancer")) celltype=CANCER;
			else if(!strcmp(&type[0], "normal")) celltype=NORMALCELL;
			else if(!strcmp(&type[0], "macrophage")) celltype=MACROPHAGE;
			else if(!strcmp(&type[0], "monocyte")) celltype=MONOCYTE;
			else if(!strcmp(&type[0], "none")) celltype=NONE;
			else printf("WARNING: Error processing initial cell config file %s - unknown cell type \"%s\"\n", filename.c_str(), &type[0]);

			if(endx==-1) endx=startx;
			if(endy==-1) endy=starty;
			if(endz==-1) endz=startz;

			for(int x=startx;x<=endx;x++)
			{
				for(int y=starty;y<=endy;y++)
				{
					for(int z=startz;z<=endz;z++)
					{

						if(Cells[x][y][z][0].Type==NONE)
						{
							MakeCell(x, y, z, 0, celltype);
							if(celltype==NORMALCELL || celltype==CANCER)
							{
								Cells[x][y][z][0].RandomisePhase(mtrand.rand());
							}
						}
						if(celltype==NONE) KillCell(x, y, z, 0);
					}
				}
			}
		}
	}
	fclose(fp);
	if(VERBOSE) printf("Closed initial cell config file \"%s\"\n", filename.c_str());

}

//-------------------------------------------------------------------
//! Main step of simulation, CellularAutomaton and Vasculature
void CellularAutomaton::MainStep(string outstem)
{
	RunStep();
	if (param_Vasculature) Vasc.OutputToFile(outstem.c_str(), "Vasculature/debug", 0);
}

//-------------------------------------------------------------------
//! Main step of simulation, CellularAutomaton and Vasculature - reordered
void CellularAutomaton::MainStepReorder(string outstem)
{
	RunStepReorder();
	if (param_Vasculature) Vasc.OutputToFile(outstem.c_str(), "Vasculature/debug", 0);
}

//-------------------------------------------------------------------
//! Main step of simulation, CellularAutomaton and Vasculature
void CellularAutomaton::MainStepConsolidated(string outstem)
{
	RunStepConsolidated();
	if (param_Vasculature) Vasc.OutputToFile(outstem.c_str(), "Vasculature/debug", 0);
}

//-------------------------------------------------------------------
//! Make a daughter cell from a parent
void CellularAutomaton::MakeCellDaughterFromParent(Location& p, Location& d)
{
	MakeCell(d.x, d.y, d.z, d.stack, Cells[p.x][p.y][p.z][p.stack].Type);
	if (VERBOSE) printf("Successful division of Type %d at (%d,%d,%d,%d), daughter at (%d,%d,%d,%d)\n",Cells[p.x][p.y][p.z][p.stack].Type,p.x,p.y,p.z,p.stack,d.x,d.y,d.z,d.stack);
}

//-------------------------------------------------------------------
//! Make a new Cell at a location of a specific type.
//! If the cell is a macrophage or monocyte, randomise its age.
void CellularAutomaton::MakeCell(int x, int y, int z, int stack, int Type)
{
	MaxCellID++;
	int CellID = MaxCellID;
	Cells[x][y][z][stack].MakeNewCell(Type, CellID);
	Cells[x][y][z][stack].RandomiseMacrophageAge(mtrand);
	ActiveCells.push_back(Location(x,y,z,stack));
}
//-------------------------------------------------------------------
//! This simply moves the info from one location to another
//! and also alters just the one entry in ActiveCells within CellMovement
//! Checks tip cells for regression.
//! This function only implements movement once a decision of where to move has been
//! made by CanMove (within CellMovement).
//!
//! Resets the old cell location using CAElement::KillElement?
//!
void CellularAutomaton::MoveCell(Location& OldPos, Location& NewPos)
{
	if (Cells[OldPos.x][OldPos.y][OldPos.z][OldPos.stack].Type == NONE)
    {
        printf("ERROR: Trying to move cell type NONE\n");
		exit(1);
	}

	// Cells(new)=Cells(old) assigns all structural elements of Old to New
	Cells[NewPos.x][NewPos.y][NewPos.z][NewPos.stack] = Cells[OldPos.x][OldPos.y][OldPos.z][OldPos.stack];


	if (Cells[OldPos.x][OldPos.y][OldPos.z][OldPos.stack].Type == TIPCELL)
		// if a tipcell, check for regression
		{
		int PVID = Cells[OldPos.x][OldPos.y][OldPos.z][OldPos.stack].VesselID;
		// check if the 2nd last endothelial cell in the ProtoVessel is at the target (x,y) position.
		// last entry in list is
		int last = Vasc.ProtoVessels[PVID].xyzstacklist.size();
		// provided there is more than 1 entry
		if (last > 1)
        {
			Location lastloc = Vasc.ProtoVessels[PVID].xyzstacklist[last-2];
			// if there IS an ENDOCELL at new position, we must have reversed
			if (NewPos.x==lastloc.x && NewPos.y==lastloc.y && NewPos.z==lastloc.z)
				{
				// Cells(new)=Cells(old) assigns all structural elements of Old to New
				// in particular, this location will inherit the Type TIPCELL
				Cells[NewPos.x][NewPos.y][NewPos.z][lastloc.stack] = Cells[OldPos.x][OldPos.y][OldPos.z][OldPos.stack];
				// no cell at old Tip location
				Cells[OldPos.x][OldPos.y][OldPos.z][OldPos.stack].KillElement();
				// no cell at new Tip location (with the different z co-ord)
				Cells[NewPos.x][NewPos.y][NewPos.z][NewPos.stack].KillElement();
				// remove last entry in ProtoVessel list
				Vasc.ProtoVessels[PVID].xyzstacklist.pop_back();
				// Update Newstack to the new TIPCELL location - not sure if this feeds back...
				NewPos.stack = lastloc.stack;
            }
			else
				{
				// old location must be made into an ENDOCELL
				// note that at the moment these are not ACTIVECELLS
				Cells[OldPos.x][OldPos.y][OldPos.z][OldPos.stack].Type=ENDOCELL;
				// Add the new location to the ProtoVessel list
				Vasc.ProtoVessels[PVID].xyzstacklist.push_back(NewPos);
            }
        }
		else
			{
			// old location must be made into an ENDOCELL
			// note that at the moment these are not ACTIVECELLS
			Cells[OldPos.x][OldPos.y][OldPos.z][OldPos.stack].Type=ENDOCELL;
			// Add the new location to the ProtoVessel list
			Vasc.ProtoVessels[PVID].xyzstacklist.push_back(NewPos);
        }
		//  either way, the tip has moved from [oldx][oldy] to [newx][newy]
		tiplocs[NewPos.x][NewPos.y][NewPos.z] = tiplocs[OldPos.x][OldPos.y][OldPos.z];
		tiplocs[OldPos.x][OldPos.y][OldPos.z] = 0;
    }
	else
		{
		Cells[OldPos.x][OldPos.y][OldPos.z][OldPos.stack].KillElement();
    }

	if (Cells[NewPos.x][NewPos.y][NewPos.z][NewPos.stack].Type==NONE)
		{
		printf("ERROR: Moved cell from (%d,%d,%d,%d) at (%d,%d,%d,%d) has type NONE\n", OldPos.x, OldPos.y, OldPos.z, OldPos.stack, NewPos.x, NewPos.y, NewPos.z, NewPos.stack);
		exit(1);
    }

	if (VERBOSE) printf("Moved cell type %d from (%d,%d,%d,%d) to (%d,%d,%d,%d)\n", Cells[NewPos.x][NewPos.y][NewPos.z][NewPos.stack].Type, OldPos.x, OldPos.y, OldPos.z, OldPos.stack, NewPos.x, NewPos.y, NewPos.z, NewPos.stack);
}

//-------------------------------------------------------------------
//! Random walk test
void CellularAutomaton::RandomWalkTest(int nsteps, string outstem)
{
	int i = 1;
	int j;
	int k;
	char filename[200];
	bool appendsubcellular = true;
	Location pos, NewPos, InitialPos = ActiveCells[0];

	sprintf(filename, "./%s/random-walk.pos", outstem.c_str());
	SingleCellToFile(filename, appendsubcellular);

	mtrand.seed(param_RandSeed);
	while(i<=param_nwrites)
	{
		for(j=1;j<=nsteps;j++)
		{
			pos = ActiveCells[0];
			if(CanMove(pos.x,pos.y,pos.z,pos.stack,NewPos))
				{
				MoveCell(pos,NewPos);
				ActiveCells[0]=NewPos;
				}
	    }
		GLOBALTIMESTEP+=(int)(nsteps*param_timestep);
		SingleCellToFile(filename, appendsubcellular);
		printf("***** End of step %d *****\n",i);
		i++;
	}

	i = 1;
	GLOBALTIMESTEP = 0;
	pos = ActiveCells[0];
	if (!(pos == InitialPos)) {
		MoveCell(pos,InitialPos);
		ActiveCells[0] = InitialPos;
		}

	sprintf(filename, "./%s/random-walk-moore.pos", outstem.c_str());
	SingleCellToFile(filename, appendsubcellular);

	mtrand.seed(param_RandSeed);
	while(i<=param_nwrites)
	{
		for(j=1;j<=nsteps;j++)
		{
			pos = ActiveCells[0];
			if(CanMoveRandomMoore(pos.x,pos.y,pos.z,pos.stack,NewPos))
				{
				MoveCell(pos,NewPos);
				ActiveCells[0]=NewPos;
				}
	    }
		GLOBALTIMESTEP+=(int)(nsteps*param_timestep);
		SingleCellToFile(filename, appendsubcellular);
		printf("***** End of step %d *****\n",i);
		i++;
	}

	i = 1;
	GLOBALTIMESTEP = 0;
	pos = ActiveCells[0];
	if (!(pos == InitialPos)) {
		MoveCell(pos,InitialPos);
		ActiveCells[0] = InitialPos;
		}

	sprintf(filename, "./%s/random-walk-neumann.pos", outstem.c_str());
	SingleCellToFile(filename, appendsubcellular);

	mtrand.seed(param_RandSeed);
	while(i<=param_nwrites)
	{
		for(j=1;j<=nsteps;j++)
		{
			pos = ActiveCells[0];
			if(CanMoveRandomNeumann(pos.x,pos.y,pos.z,pos.stack,NewPos))
				{
				MoveCell(pos,NewPos);
				ActiveCells[0]=NewPos;
				}
	    }
		GLOBALTIMESTEP+=(int)(nsteps*param_timestep);
		SingleCellToFile(filename, appendsubcellular);
		printf("***** End of step %d *****\n",i);
		i++;
	}

}

//-------------------------------------------------------------------
//! Remove vessels with low flow (TWSSsubthreshold > Tdie)
bool CellularAutomaton::RemoveLowFlowVessels()
{
	bool RemovedVessels = false;
	int i,k;

	for(k=0;k<Vasc.ActiveVessels.size();k++)
	{
		i = Vasc.ActiveVessels[k];
		if(Vasc.Vessels[i].TWSSsubthreshold > param_LowFlowDeathTime)
		{
			if (VERBOSE) printf("VVVVVVVVVV Removing low flow vessel from (%d,%d,%d) -> (%d,%d,%d)\n",Vasc.Nodes[Vasc.Vessels[i].Start].x,Vasc.Nodes[Vasc.Vessels[i].Start].y,Vasc.Nodes[Vasc.Vessels[i].Start].z,Vasc.Nodes[Vasc.Vessels[i].End].x,Vasc.Nodes[Vasc.Vessels[i].End].y,Vasc.Nodes[Vasc.Vessels[i].End].z);
			Vasc.RemoveVessel(i, REMOVELOWFLOWVESSELS);
			RemoveVesselFromCA(i, REMOVELOWFLOWVESSELS);
			RemovedVessels = true;
			k--;
		}
	}
	return RemovedVessels;
}

//-------------------------------------------------------------------
//! Remove vessels coopted vessels
bool CellularAutomaton::RemoveCooptedVessels()
{
	bool RemovedVessels = false;
	int i,k;
	int newstate;

	if(VesselRemoval) for(k=0;k<Vasc.ActiveVessels.size();k++)
	{
		i = Vasc.ActiveVessels[k];
		newstate = Vasc.Vessels[i].State;
		if(newstate == COLLAPSEVESSEL && Vasc.Vessels[i].WallStress < param_TWSSCollapse)
		{
			Vasc.RemoveVessel(i, REMOVECOOPTEDVESSELS);
			RemoveVesselFromCA(i, REMOVECOOPTEDVESSELS);
		        printf("Coopted Vessel %d has been removed\n",i);
			RemovedVessels = true;
			k--;
		}
	}
	return RemovedVessels;
}

//-------------------------------------------------------------------
//! Remove vessels collapsed due to proliferative pressure
bool CellularAutomaton::RemoveVesselsCollapsedByProliferatingCells()
{
	bool RemovedVessels = false;
	int i,k;

	for(k=0;k<Vasc.ActiveVessels.size();k++)
	{
		i = Vasc.ActiveVessels[k];
		if(Vasc.Vessels[i].compressedByProliferatingCells == 1)
		{
			if (VERBOSE) printf("VVVVVVVVVV Removing vessel collapsed by proliferating cells from (%d,%d,%d) -> (%d,%d,%d)\n",Vasc.Nodes[Vasc.Vessels[i].Start].x,Vasc.Nodes[Vasc.Vessels[i].Start].y,Vasc.Nodes[Vasc.Vessels[i].Start].z,Vasc.Nodes[Vasc.Vessels[i].End].x,Vasc.Nodes[Vasc.Vessels[i].End].y,Vasc.Nodes[Vasc.Vessels[i].End].z);
			Vasc.RemoveVessel(i, REMOVELOWFLOWVESSELS);
			RemoveVesselFromCA(i, REMOVELOWFLOWVESSELS);
			RemovedVessels = true;
			k--;
		}
	}
	return RemovedVessels;
}

//-------------------------------------------------------------------
//! Remove vessels collapsed due to proliferative pressure
bool CellularAutomaton::RemoveVesselsCollapsedByEmbolisationTherapy()
{
	bool RemovedVessels = false;
	int i,k;


	for(k=0;k<Vasc.ActiveVessels.size();k++)
	{
		i = Vasc.ActiveVessels[k];
		if(Vasc.Vessels[i].collapsedByEmbolisationTherapy == 1)
		{
			// Indicate the positions where there are TACEParticles... (Occupied positions have the value 1.0)
			TACEParticlessum[Vasc.Nodes[Vasc.Vessels[i].Start].x][Vasc.Nodes[Vasc.Vessels[i].Start].y][Vasc.Nodes[Vasc.Vessels[i].Start].z]=1.0;
			TACEParticlessum[Vasc.Nodes[Vasc.Vessels[i].End].x][Vasc.Nodes[Vasc.Vessels[i].End].y][Vasc.Nodes[Vasc.Vessels[i].End].z]=1.0;

			if (VERBOSE) printf("VVVVVVVVVV Removing vessel collapsed by embolisation therapy from (%d,%d,%d) -> (%d,%d,%d)\n",Vasc.Nodes[Vasc.Vessels[i].Start].x,Vasc.Nodes[Vasc.Vessels[i].Start].y,Vasc.Nodes[Vasc.Vessels[i].Start].z,Vasc.Nodes[Vasc.Vessels[i].End].x,Vasc.Nodes[Vasc.Vessels[i].End].y,Vasc.Nodes[Vasc.Vessels[i].End].z);
			Vasc.RemoveVessel(i, REMOVELOWFLOWVESSELS);
			RemoveVesselFromCA(i, REMOVELOWFLOWVESSELS);
			RemovedVessels = true;
			k--;
		}
	}
	return RemovedVessels;
}

//-------------------------------------------------------------------
//! Removes a single ProtoVessel (given as an index into the ProtoVessels array),
//! from the CellularAutomaton CA. This is done using the ProtoVessel's xyzlist,
//! removing cells of type ENDOCELL and TIPCELL from CA's Cells grid.
//! This function writes over the contents of cells in the grid.
//!
void CellularAutomaton::RemoveProtoVesselFromCA(int v)
{
	// the map is stored as Locations
	// in ProtoVessels[v].xyzlist

	Location pos;
	int curx, cury, curz, k;

	// reset the tail to zero
	taillocs[Vasc.ProtoVessels[v].xyzstacklist[0].x][Vasc.ProtoVessels[v].xyzstacklist[0].y][Vasc.ProtoVessels[v].xyzstacklist[0].z]=0;

	for(int i=0;i<Vasc.ProtoVessels[v].xyzstacklist.size();i++)
	{
		pos = Vasc.ProtoVessels[v].xyzstacklist[i];

		// check that it really is an ENDOCELL
		if ((Cells[pos.x][pos.y][pos.z][pos.stack].Type == ENDOCELL) || (Cells[pos.x][pos.y][pos.z][pos.stack].Type == TIPCELL))
		{
			Cells[pos.x][pos.y][pos.z][pos.stack].Type = NONE;
			ActiveCells.erase(std::remove(ActiveCells.begin(), ActiveCells.end(), pos), ActiveCells.end());
		}

	}

}


//-------------------------------------------------------------------
//! Removes a single vessel, vessel (given as an index into the Vessels array),
//! from the CellularAutomaton CA. This is done by
//! removing cells of type VESSEL from CA's Cells grid.
//! This function writes over the contents of cells in the grid.
//!
//! Should also correctly set the NumBonds property.
//! NumBonds has to be updated since Vasculature::RemoveVessel changes the number of bonds
//! by altering NumVesselsIn and NumVesselsOut.
//! NumBonds was set to 3 for Inflows and Outflows in order to not allow sprouts from them,
//! but here we need NumBonds to properly reflect the number of bonds even at inflows and outflows.
//! We do this here, but since Vasculature::AddToCA is called by Vasculature::Update, NumBonds will
//! be reset to 3 for inflows and outflows at this point.
//!
//! @todo Tidy up this part (setting NumBonds), see also Vasculature::AddVesselToCA.
//! @todo At some point we might want to allow new connections to an inflow or outflow.
//! This would require a change so that removing a vessel connecting to an inflow or outflow
//! leaves behind a vessel cell that can be reconnected to.
//!
//! @param CallingFunction is an integer, defined in Definitions.h, which usually is
//! one of REMOVELOWFLOWVESSELS, EXTRAVASATIONFROMVESSEL and ANASTOMOSIS.
//!
void CellularAutomaton::RemoveVesselFromCA(int v, int CallingFunction)
{

	Location pos;
	int curx, cury, curz, k;

	// First update the NumBonds property at either end.
	// This will then show whether we should really remove the cells at either end:
	// If NumBonds > 0 then we should preserve the CA cells at that end and set the VesselID accordingly.

	////////////////////////////////////////////////////
	// START NODE
	////////////////////////////////////////////////////
	// the startnode is a Node
	int startnode = Vasc.Vessels[v].Start;
	curx = Vasc.Nodes[startnode].x;
	cury = Vasc.Nodes[startnode].y;
	curz = Vasc.Nodes[startnode].z;
	// recall that the NumBonds property records the number of vessels connected to a Node
	// and is hence zero if the CA location is NOT a Node
	Cells[curx][cury][curz][0].NumBonds = Vasc.Nodes[startnode].NumVesselsIn + Vasc.Nodes[startnode].NumVesselsOut;
	if (taillocs[curx][cury][curz] > 0) Cells[curx][cury][curz][0].NumBonds++;
	// If NumBonds > 3 something must have gone wrong.
	if (Cells[curx][cury][curz][0].NumBonds > 3)
		printf("Node %d, start of vessel %d at (%d,%d,%d), NumBonds = %d > 3\n",startnode, v, curx, cury, curz, Cells[curx][cury][curz][0].NumBonds);
	// check that it is not a Node which is part of another vessel, unless we are splitting a vessel
	if ((Cells[curx][cury][curz][0].NumBonds == 0) || (CallingFunction == EXTRAVASATIONFROMVESSEL) || (CallingFunction == ANASTOMOSIS) || (taillocs[curx][cury][curz] != 0))
	{
		KillCell(curx,cury,curz,0);
	}
	// check that it's VesselID is not that of the removed Vessel
	else if (v == Cells[curx][cury][curz][0].VesselID)
	{
		if (VERBOSE) printf("In RemoveVesselFromCA for startnode, Cells[%d][%d][%d][0].VesselID == killVID = %d. Need to change this.\n",curx,cury,curz,v);
		for (int i=0;i<Vasc.Nodes[startnode].NumVesselsIn;i++)
		{
			if (Vasc.Nodes[startnode].VesselsIn[i] != v) Cells[curx][cury][curz][0].VesselID = Vasc.Nodes[startnode].VesselsIn[i];
		}
		for (int i=0;i<Vasc.Nodes[startnode].NumVesselsOut;i++)
		{
			if (Vasc.Nodes[startnode].VesselsOut[i] != v) Cells[curx][cury][curz][0].VesselID = Vasc.Nodes[startnode].VesselsOut[i];
		}
		if (v == Cells[curx][cury][curz][0].VesselID)
		{
			printf("In RemoveVesselFromCA for startnode, Cells[%d][%d][%d][0].VesselID == killVID = %d, even after trying to change to another connected vessel.\n",curx,cury,curz,v);
			exit(1);
		}
	}

	////////////////////////////////////////////////////
	// END NODE
	////////////////////////////////////////////////////
	// the endnode is a Node
	int endnode = Vasc.Vessels[v].End;
	curx = Vasc.Nodes[endnode].x;
	cury = Vasc.Nodes[endnode].y;
	curz = Vasc.Nodes[endnode].z;
	Cells[curx][cury][curz][0].NumBonds = Vasc.Nodes[endnode].NumVesselsIn + Vasc.Nodes[endnode].NumVesselsOut;
	if (taillocs[curx][cury][curz] > 0) Cells[curx][cury][curz][0].NumBonds++;
	// If NumBonds > 3 something must have gone wrong.
	if (Cells[curx][cury][curz][0].NumBonds > 3)
		printf("Node %d, end of vessel %d at (%d,%d,%d), NumBonds = %d > 3\n",endnode, v, curx, cury, curz, Cells[curx][cury][curz][0].NumBonds);
	// check that it is not a Node which is part of another vessel, unless we are splitting a vessel, or at a Protovessel tail
	if ((Cells[curx][cury][curz][0].NumBonds == 0) || (CallingFunction == EXTRAVASATIONFROMVESSEL) || (CallingFunction == ANASTOMOSIS) || (taillocs[curx][cury][curz] != 0))
	{
		KillCell(curx,cury,curz,0);
	}
	// check that it's VesselID is not that of the removed Vessel
	else if (v == Cells[curx][cury][curz][0].VesselID)
	{
		if (VERBOSE) printf("In RemoveVesselFromCA for endnode, Cells[%d][%d][%d][0].VesselID == killVID = %d. Need to change this.\n",curx,cury,curz,v);
		for (int i=0;i<Vasc.Nodes[endnode].NumVesselsIn;i++)
		{
			if (Vasc.Nodes[endnode].VesselsIn[i] != v) Cells[curx][cury][curz][0].VesselID = Vasc.Nodes[endnode].VesselsIn[i];
		}
		for (int i=0;i<Vasc.Nodes[endnode].NumVesselsOut;i++)
		{
			if (Vasc.Nodes[endnode].VesselsOut[i] != v) Cells[curx][cury][curz][0].VesselID = Vasc.Nodes[endnode].VesselsOut[i];
		}
		if (v == Cells[curx][cury][curz][0].VesselID)
		{
			printf("In RemoveVesselFromCA for endnode, Cells[%d][%d][%d][0].VesselID == killVID = %d, even after trying to change to another connected vessel.\n",curx,cury,curz,v);
			exit(1);
		}
	}

}

//-------------------------------------------------------------------
//! Run simulation
//! NEED TO AVOID PASSING THE CA AND VASC AROUND ...
void CellularAutomaton::RunSimulation(int nsumsteps, int nsteps, string outstem)
{
	// Calculate initial stats and write initial state to files
	int i = 0;
	CalculateStatistics();
	WriteSummary(outstem.c_str(), i);
	WriteOutput(outstem.c_str(), i);
	if (param_Vasculature) {
		Vasc.OutputToFile(outstem.c_str(), "Vasculature/frame", i);
		Vasc.SegmentLengthsToFile(outstem.c_str(), "Segments/frame", i);
		if (param_WritePovRay) Vasc.OutputToPovRayFile(outstem.c_str(), "PovRay/VascPOV",i,*this);
		if (param_WriteOpenInventor) Vasc.OutputToOpenInventorFile(outstem.c_str(), "OpenInventor/VascOpenInv",i,*this);
	}
	printf("End of initialisation, i=%d, nsumsteps=%d, nsteps=%d\n",i,nsumsteps,nsteps);

	i = 1;
	while(i<=param_nwrites)
	{
		printf("###### GLOBALTIMESTEP = %d ****** MainStep = %s ******\n",GLOBALTIMESTEP,param_MainStep.c_str());

		for(int j=1;j<=nsteps;j++)
		{
			for(int k=1;k<=nsumsteps;k++)
			{
				if (param_MainStep == "reorder")
				{
					MainStepReorder(outstem);
				}
				else if (param_MainStep == "consolidated")
				{
					MainStepConsolidated(outstem);
				}
				else
				{
					MainStep(outstem);
				}
				GLOBALTIMESTEP+=(int)param_timestep;
			}
			if(VERBOSE) DEBUG_Print();

			CalculateStatistics();
			//only matters if i=0 to initialise the file
			WriteSummary(outstem.c_str(), i);
	    }
		WriteOutput(outstem.c_str(), i);
		if (param_Vasculature)
		{
			Vasc.OutputToFile(outstem.c_str(), "Vasculature/frame", i);
			Vasc.SegmentLengthsToFile(outstem.c_str(), "Segments/frame", i);
			if (param_WritePovRay) Vasc.OutputToPovRayFile(outstem.c_str(), "PovRay/VascPOV",i,*this);
			if (param_WriteOpenInventor) Vasc.OutputToOpenInventorFile(outstem.c_str(), "OpenInventor/VascOpenInv",i,*this);
		}
		printf("***** End of steps to spatial output %d *****\n",i);
		i++;
	}
}

//-------------------------------------------------------------------
//! @brief Run a single step of the Cellular Automaton.
//! @todo This function includes a check whether all active cells
//! really are active, which may no longer be necessary.
void CellularAutomaton::RunStep()
{
	clock_t mystart = clock();

	CheckForNONEInActiveCells();

	if (param_isOxygenDiffusible){SetOxygenSourceTerms(); Oxygen.Update(); SetCellOxygens();}

	UpdatetotalMass();
	UpdateCells();
	shuffleActiveCells();
	CellDivision();
	shuffleActiveCells();
	CellMovement();
	UpdateThresholds();

	if (param_isVEGFDiffusible){SetVEGFSourceTerms(); VEGF.Update(); SetExternalVEGF();}
	if (param_isProDrugDiffusible){SetProDrugSourceTerms(); ProDrug.Update();}
	if (param_isDrugDiffusible){SetDrugSourceTerms(); Drug.Update(); SetCellIntercalated();}

	CellQuiescence();
	CellDeath();
	AgeCells(param_timestep);

	clock_t mystop = clock();
	TimeRunStep += (mystop-mystart)/(double)CLOCKS_PER_SEC;
	if (VERBOSE) printf("@@@@@@ TIMING @@@@@@: Time spent in CellularAutomaton::RunStep = %g seconds\n",TimeRunStep);

	if (param_Vasculature) {UpdateVasculature(); AddVesselsToCA();}
}

//-------------------------------------------------------------------
//! @brief Run a single step of the Cellular Automaton.
//! Consolidated order, e.g. with all diffusibles updated together.
//! @todo This function includes a check whether all active cells
//! really are active, which may no longer be necessary.
void CellularAutomaton::RunStepConsolidated()
{
	clock_t mystart = clock();

	CheckForNONEInActiveCells();

	UpdatetotalMass();
	UpdateCells();
	shuffleActiveCells();
	CellDivision();
	shuffleActiveCells();
	CellMovement();
	UpdateThresholds();

	UpdateDiffusibles();

	CellQuiescence();
	CellDeath();
	AgeCells(param_timestep);

	clock_t mystop = clock();
	TimeRunStep += (mystop-mystart)/(double)CLOCKS_PER_SEC;
	if (VERBOSE) printf("@@@@@@ TIMING @@@@@@: Time spent in CellularAutomaton::RunStepConsolidated = %g seconds\n",TimeRunStep);

	if (param_Vasculature) {UpdateVasculature(); AddVesselsToCA();}
}

//-------------------------------------------------------------------
//! @brief Run a single (reordered) step of the Cellular Automaton.
//! @todo This function includes a check whether all active cells
//! really are active, which may no longer be necessary.
//! @todo New daughter cells get OxygenLevel=param_Initial_OxygenLevel which is not updated
//! before the cell cycle is updated by UpdateCells.
void CellularAutomaton::RunStepReorder()
{
	if (param_Vasculature) {UpdateVasculature(); AddVesselsToCA();}

	clock_t mystart = clock();

	CheckForNONEInActiveCells();

	if (param_isVEGFDiffusible){SetVEGFSourceTerms(); VEGF.Update(); SetExternalVEGF();}

	UpdatetotalMass();
	CellMovement();
	UpdateCells();
	UpdateThresholds();

	if (param_isOxygenDiffusible){SetOxygenSourceTerms(); Oxygen.Update(); SetCellOxygens();}
	if (param_isProDrugDiffusible){SetProDrugSourceTerms(); ProDrug.Update();}
	if (param_isDrugDiffusible){SetDrugSourceTerms(); Drug.Update(); SetCellIntercalated();}

	shuffleActiveCells();
	CellDivision();
	shuffleActiveCells();
	CellQuiescence();
	CellDeath();
	AgeCells(param_timestep);

	clock_t mystop = clock();
	TimeRunStep += (mystop-mystart)/(double)CLOCKS_PER_SEC;
	if (VERBOSE) printf("@@@@@@ TIMING @@@@@@: Time spent in CellularAutomaton::RunStepReorder = %g seconds\n",TimeRunStep);
}

//-------------------------------------------------------------------
//! @brief Run a single step of the Cellular Automaton with only a cell cycle model.
//! Uses UpdateCells() to step the cell cycle model forward, but does not actually
//! divide cells, - only uses DivisionReset
void CellularAutomaton::RunCellCycleStep()
{
	Location pos;

	CheckForNONEInActiveCells();

	UpdateCells();
	int i, imax=ActiveCells.size();
	for(i=0;i<imax;i++)
	{
		pos = ActiveCells[i];

		if(Cells[pos.x][pos.y][pos.z][pos.stack].InternalCellDivision())
		{
			//Reset current cell
			Cells[pos.x][pos.y][pos.z][pos.stack].DivisionReset();
		}
	}
	CellQuiescence();
	AgeCells(param_timestep);
}

//-------------------------------------------------------------------
//! Set the internal Intercalated switch according to the drug distribution.
//! Loop through active cells and set Intercalated = true if drug > threshold.
//!
void CellularAutomaton::SetCellIntercalated()
{
	int NumActiveCells = ActiveCells.size();
	Location pos;

	for(int i=0;i<NumActiveCells;i++)
	{
		pos = ActiveCells[i];
		if (Drug.GetLevelDirect(pos.x, pos.y, pos.z) > param_DrugKillThreshold) Cells[pos.x][pos.y][pos.z][pos.stack].Intercalated = true;
	}
}

//-------------------------------------------------------------------
//! Set the internal OxygenLevel according to the Diffusible distribution.
//! Loop through active cells and set OxygenLevel = Oxygen.GetLevelDirect.
//!
//! @todo Could this be part of the Diffusible itself, by acting on CA.
void CellularAutomaton::SetCellOxygens()
{
	int NumActiveCells = ActiveCells.size();
	Location pos;

	for(int i=0;i<NumActiveCells;i++)
	{
		pos = ActiveCells[i];
		Cells[pos.x][pos.y][pos.z][pos.stack].OxygenLevel = Oxygen.GetLevelDirect(pos.x, pos.y, pos.z);
	}
}

//-------------------------------------------------------------------
//! Set the ExternalVEGF according to the Diffusible distribution.
//! Loop through active cells and set ExternalVEGF = VEGF.GetLevelDirect.
//!
void CellularAutomaton::SetExternalVEGF()
{
	int NumActiveCells = ActiveCells.size();
	Location pos;

	for(int i=0;i<NumActiveCells;i++)
	{
		pos = ActiveCells[i];
		Cells[pos.x][pos.y][pos.z][pos.stack].ExternalVEGF = VEGF.GetLevelDirect(pos.x, pos.y, pos.z);
	}
}

//------------------------------------------------------------------- OVERRIDE - NEW
//! Set the internal DrugLevel according to the Diffusible distribution.
//! Loop through active cells and set DrugLevel = Drug.GetLevelDirect.
//!
void CellularAutomaton::SetCellDrug()
{
    int NumActiveCells = ActiveCells.size();
    Location pos;

    for(int i=0; i<NumActiveCells; i++)
    {
        pos = ActiveCells[i];
        Cells[pos.x][pos.y][pos.z][pos.stack].DrugLevel = Drug.GetLevelDirect(pos.x, pos.y, pos.z);
    }
}

//-------------------------------------------------------------------
//! Set the Proportional source term at the grid element xgrid, ygrid.
void CellularAutomaton::SetVascPropSource(int xgrid, int ygrid, int zgrid, double value)
{
//	VascPropSource.SetElement(ELEMENT(xgrid,ygrid,zgrid), value);
	VascPropSource[ELEMENT(xgrid,ygrid,zgrid)] = value;
}

//-------------------------------------------------------------------
void CellularAutomaton::shuffleActiveCells()
{
	int last;
	int size = ActiveCells.size();
	Location pos;

	if (VERBOSE) printf("Starting CellularAutomaton::shuffleActiveCells of %d Members\n", size);

	for (last = size; last > 1; last--)
	{
		int RANDSTATE = mtrand.randInt(last-1);
		pos = ActiveCells[RANDSTATE];
		ActiveCells[RANDSTATE] = ActiveCells[last - 1];
		ActiveCells[last - 1] = pos;
	}
}// end shuffleCells( )

//-------------------------------------------------------------------
//! For large transitions Probabilities it is important that the directions are shuffled,
//! if not the cells move only in one direction because Pij>1
void CellularAutomaton::shuffleLocationsAndProbabilities(Location Test[], double Pij[], const int NumTests)
{
	int last;
	int size = NumTests;
	Location pos;
	double P;
#if 0
	printf("Before ShuffleDirectionsAndProbabilities\n");
	for (int i=0; i<NumTests; ++i)
	{
		printf("Location: %d  %d  %d    Probability= %.3f\n",Test[i].x,Test[i].y,Test[i].z,Pij[i]);
	}
#endif
	if (VERBOSE) printf("Starting CellularAutomaton::shuffleDirectionsAndProbabilities of %d Members\n", size);

	for (last = size; last > 1; last--)
	{
		//RANDSTATE = rand() % last; // in the range from 0 to last-1: http://www.cplusplus.com/reference/clibrary/cstdlib/rand.html
		int RANDSTATE = mtrand.randInt(last-1);
		//           pos = ActiveCells[randomNum];
		pos = Test[RANDSTATE];
		P=Pij[RANDSTATE];
		//           ActiveCells[randomNum] = ActiveCells[last - 1];
		Test[RANDSTATE] = Test[last - 1];
		Pij[RANDSTATE] = Pij[last - 1];
		//           ActiveCells[last - 1] = pos;
		Test[last - 1] = pos;
		Pij[last -1] = P;
	}
#if 0
	printf("After ShuffleDirectionsAndProbabilities\n");
	for (int i=0; i<NumTests; ++i)
	{
		printf("Location: %d  %d  %d    Probability= %.3f\n",Test[i].x,Test[i].y,Test[i].z,Pij[i]);
	}
	std::cin.get();
#endif
}// end shuffleDirectionsAndProbabilities( )


//-------------------------------------------------------------------
//! Count the number of cells at (x,y,z) and update emptySites such that the
//! stack location is empty (if there is such an empty space)
int CellularAutomaton::StackCount(int x, int y, int z, Location &emptySites)
{
	int filled=0;
	for (int k=0; k<StackCells; k++)
	{
		// count cells at test element
		if(Cells[x][y][z][k].Type!=NONE)
		// if this location is not empty
		{
			filled++;
		}
		else
		// if it is empty, set this to be the new position
		// I think this always puts new cell at top of the stack
		{
			emptySites.stack=k;
		}
	}
	return filled;
}

//-------------------------------------------------------------------
//! Function to read cell type and subcellular information from output files.
//! Creates new cells of the appropriate type at those locations.
//!
//! Checks that none of the locations (x,y,z) lie outside the declared range
void CellularAutomaton::SubCellularFromFile(string filename)
{
	ifstream fp;
	string folderFilename = "input/" + filename;
	fp.open(folderFilename.c_str(),  std::ifstream::in);
	int Type;
	int x, y, z, stack;
	int CellID;
	int i = 0;
	std::vector<int> aci;
	std::vector<Location> acl;

	// file is bad so print error and exit
	if( !fp.good() )
		{
		cout << "ERROR: Error opening subcellular file " << filename << "\n";
		exit(1);
		}
	else
		{
			cout << "************************************************************\n";
			cout << "*** Reading subcellular file " << filename << "\n";
			cout << "************************************************************\n";
		}

	char buf[256];

	while( fp.getline(buf, sizeof(buf)) )
		{
		if(VERBOSE) printf("Reading line %d\n", i++);
		istringstream iss(buf);
		int CellID;
		int celltype;
		int ac;
		iss >> x >> y >> z >> stack >> CellID >> celltype;
 		Cells[x][y][z][stack].MakeNewCell(celltype, CellID);
		iss.seekg ( 0, std::ios::beg );
		iss >> x >> y >> z >> stack >> Cells[x][y][z][stack] >> ac;
		// NEED TO RECONSTRUCT ActiveCells...
		// ac is the index into ActiveCells from the input data
		// aci is a vector of those indices
		// acl is a vector of their Locations
		// Then we should construct ActiveCells such that
		// ActiveCells[i] = acl[aci.IsMember(i)];
		if (ac != -1) {
			MaxCellID++;
			aci.push_back(ac);
			acl.push_back(Location(x,y,z,stack));
			}
		}

	// NEED TO RECONSTRUCT ActiveCells...
	for (i=0;i<aci.size();i++) {
		std::vector<int>::iterator it = find (aci.begin(), aci.end(), i);
		if ( it != aci.end()) {
			int j = it - aci.begin();
			ActiveCells.push_back(acl[j]);
			}
		else {
			printf("Cannot find ActiveCells member %d in Loaded Cells. Exiting.\n", i);
			exit(1);
		}
	}

	fp.close();

	// Remove any ENDOCELLs or TIPCELLs
	Location pos;
	for (i=0;i<ActiveCells.size();i++)
	  {
	    pos = ActiveCells[i];
	    if (Cells[pos.x][pos.y][pos.z][pos.stack].Type == ENDOCELL || Cells[pos.x][pos.y][pos.z][pos.stack].Type == TIPCELL)
	      {
		KillCell(pos.x,pos.y,pos.z,pos.stack);
		i--;
	      }
	  }
	for(int x=0;x<XCells;x++)
	{
		for(int y=0;y<YCells;y++)
		{
			for(int z=0;z<ZCells;z++)
			{
				for(int stack=0;stack<StackCells;stack++)
				{
				  if (Cells[x][y][z][stack].Type == ENDOCELL || Cells[x][y][z][stack].Type == TIPCELL)
				    {
				      printf("Killing EC at (%d,%d,%d,%d) &&**&&**&&**&&**&&**&&\n",x,y,z,stack);
				      KillCell(x,y,z,stack);
				    }
				}
			}
		}
	}
}

//-------------------------------------------------------------------
//! This should write the location and the data for any CAElement at that location.
//! When running a cell cycle test (CellCycleTest), and we want a single file of data
//! for a single cell, this should be called with append=true.
//!
//! @todo Should we bother writing NONE's to file?
//! This really depends on what we do with the data.
void CellularAutomaton::SubCellularToFile(char* filename, bool append)
{
	ofstream subcell;

	int xx;
	if (append && GLOBALTIMESTEP != 0) {
		subcell.open(filename, ios::app);
	}
	else {
		subcell.open(filename, ios::trunc);
	}
	for(int x=0;x<XCells;x++)
	{
		for(int y=0;y<YCells;y++)
		{
			for(int z=0;z<ZCells;z++)
			{
			    for(int stack=0; stack<StackCells;stack++)
			    {
                    if (append) {
                        xx = GLOBALTIMESTEP;
                    }
                    else {
                        xx = x;
                    }
                    // find location in ActiveCells
					std::vector<Location>::iterator it = std::find(ActiveCells.begin(), ActiveCells.end(), Location(x,y,z,stack));
					int curActiveCell = it - ActiveCells.begin();
					if (it == ActiveCells.end()) curActiveCell = -1;
                    subcell << xx << " " << y << " " << z << " " << stack << " " << Cells[x][y][z][stack] << " " << curActiveCell << "\n";
			    }
			}
		}
	}

	subcell.close();
}

//-------------------------------------------------------------------
//! This should write the location and the data for any CAElement at that location.
//! When running a cell cycle test (CellCycleTest), and we want a single file of data
//! for a single cell, this should be called with append=true.
//!
void CellularAutomaton::SingleCellToFile(char* filename, bool append)
{
	ofstream subcell;

	int xx;
	if (append && GLOBALTIMESTEP != 0) {
		subcell.open(filename, ios::app);
	}
	else {
		subcell.open(filename, ios::trunc);
	}

	for(int i=0;i<ActiveCells.size();i++)
	{
		Location pos = ActiveCells[i];
		if (Cells[pos.x][pos.y][pos.z][pos.stack].Type != NONE)
		{
			subcell << GLOBALTIMESTEP << " " << pos.x << " " << pos.y << " " << pos.z << " " << pos.stack << "\n";
		}
	}

	subcell.close();
}

//-------------------------------------------------------------------
//! Write current and cumulative spatial data to file
void CellularAutomaton::sumToFile(char* filename)
{
	FILE* fp = fopen(filename, "w");
	for(int x=0;x<XCells;x++)
	{
		for(int y=0;y<YCells;y++)
		{
			for(int z=0;z<ZCells;z++)
			{
				fprintf(fp, "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %g, %d",
				divsum[x][y][z], macsum[x][y][z], monosum[x][y][z], 0, 0, 0, 0,
				0, 0, aposum[x][y][z], diesum[x][y][z], killsum[x][y][z], taillocs[x][y][z], tiplocs[x][y][z],
				drugsum[x][y][z], tipcelloutsum[x][y][z]);
				fprintf(fp, "\n");
			}
		}
	}
	fclose(fp);
}

//-------------------------------------------------------------------
//! If there is a VESSEL at a lattice site, find next nearest neighbours.
//! Used when checking for space to divide, and updating thresholds
//!
//! @todo Should probably not do this - but eliminating this would modify Cancer Res results
int CellularAutomaton::TransVesselLocations(int x, int y, int z, int txc, int tyc, int tzc, Location (&test)[3])
{
	int Numtests;
	int tx, ty, tz;
	// Making sure we don't look outside the grid
	int minxc = MAX(x-param_DivisionRadius, 0);
	int maxxc = MIN(x+param_DivisionRadius, XCells-1);
	int minyc = MAX(y-param_DivisionRadius, 0);
	int maxyc = MIN(y+param_DivisionRadius, YCells-1);
	int minzc = MAX(z-param_DivisionRadius,0);
	int maxzc = MIN(z+param_DivisionRadius, ZCells-1);

	int minx = minxc;
	int maxx = maxxc;
	int miny = minyc;
	int maxy = maxyc;
	int minz = minzc;
	int maxz = maxzc;

	if(PERIODIC_X)
	{
		minxc = x-param_DivisionRadius;
		minx = MOD(minxc,(XCells));
		maxxc = x+param_DivisionRadius;
		maxx = MOD(maxxc,(XCells));
	}
	if(PERIODIC_Y)
	{
		minyc = y-param_DivisionRadius;
		miny = MOD(minyc,(YCells));
		maxyc = y+param_DivisionRadius;
		maxy = MOD(maxyc,(YCells));
	}
	if(PERIODIC_Z)
	{
		minzc = z-param_DivisionRadius;
		minz = MOD(minzc,(ZCells));
		maxzc = z+param_DivisionRadius;
		maxz = MOD(maxzc,(ZCells));
	}

	// Do we consider a 2D or 3D domain?
	bool is3Dx = true;
	bool is3Dy = true;
	bool is3Dz = true;
	if (minxc==maxxc) is3Dx=false;
	if (minyc==maxyc) is3Dy=false;
	if (minzc==maxzc) is3Dz=false;

	tx=txc;
	ty=tyc;
	tz=tzc;

	if(PERIODIC_X) { tx=MOD(txc,(XCells));}
	if(PERIODIC_Y) { ty=MOD(tyc,(YCells));}
	if(PERIODIC_Z) { tz=MOD(tzc,(ZCells));}

	if(Cells[tx][ty][tz][0].Type==VESSEL)
		    {
		       	if(tx==minx && ty==miny && tz==minz && is3Dx && is3Dy && is3Dz)
				{
					test[0].x = minx-1; test[0].y=miny;   test[0].z=minz;
					test[1].x = minx;   test[1].y=miny-1; test[1].z=minz;
					test[2].x = minx;   test[2].y=miny;   test[2].z=minz-1;
					Numtests=3;
				}
				else if(tx==maxx && ty==miny && tz==minz && is3Dx && is3Dy && is3Dz)
				{
					test[0].x = maxx+1; test[0].y=miny;     test[0].z=minz;
					test[1].x = maxx;   test[1].y=miny-1;   test[1].z=minz;
					test[2].x = maxx;   test[2].y=miny;     test[2].z=minz-1;
					Numtests=3;
				}
				else if(tx==minx && ty==maxy && tz==minz && is3Dx && is3Dy && is3Dz)
				{
					test[0].x = minx-1; test[0].y=maxy;     test[0].z=minz;
					test[1].x = minx;   test[1].y=maxy+1;   test[1].z=minz;
					test[2].x = minx;   test[2].y=maxy;     test[2].z=minz-1;
					Numtests=3;
				}
				else if(tx==minx && ty==miny && tz==maxz && is3Dx && is3Dy && is3Dz)
				{
					test[0].x = minx-1; test[0].y=miny;     test[0].z=maxz;
					test[1].x = minx;   test[1].y=miny-1;   test[1].z=maxz;
					test[2].x = minx;   test[2].y=miny;     test[2].z=maxz+1;
					Numtests=3;
				}
				else if(tx==maxx && ty==maxy && tz==minz && is3Dx && is3Dy && is3Dz)
				{
					test[0].x = maxx+1; test[0].y=maxy;     test[0].z=minz;
					test[1].x = maxx;   test[1].y=maxy+1;   test[1].z=minz;
					test[2].x = maxx;   test[2].y=maxy;     test[2].z=minz-1;
					Numtests=3;
				}
				else if(tx==maxx && ty==miny && tz==maxz && is3Dx && is3Dy && is3Dz)
				{
					test[0].x = maxx+1; test[0].y=miny;     test[0].z=maxz;
					test[1].x = maxx;   test[1].y=miny-1;   test[1].z=maxz;
					test[2].x = maxx;   test[2].y=miny;     test[2].z=maxz+1;
					Numtests=3;
				}
				else if(tx==minx && ty==maxy && tz==maxz && is3Dx && is3Dy && is3Dz)
				{
					test[0].x = minx-1; test[0].y=maxy;     test[0].z=maxz;
					test[1].x = minx;   test[1].y=maxy+1;   test[1].z=maxz;
					test[2].x = minx;   test[2].y=maxy;     test[2].z=maxz+1;
					Numtests=3;
				}
				else if(tx==maxx && ty==maxy && tz==maxz && is3Dx && is3Dy && is3Dz)
				{
					test[0].x = maxx+1; test[0].y=maxy;     test[0].z=maxz;
					test[1].x = maxx;   test[1].y=maxy+1;   test[1].z=maxz;
					test[2].x = maxx;   test[2].y=maxy;     test[2].z=maxz+1;
					Numtests=3;
				}
				else if(tx==minx && ty==miny && is3Dx && is3Dy)
				{
                    test[0].x = minx-1; test[0].y=miny;     test[0].z=tz;
					test[1].x = minx;   test[1].y=miny-1;   test[1].z=tz;
					Numtests=2;
				}
				else if(tx==minx && ty==maxy && is3Dx && is3Dy)
				{
                    test[0].x = minx-1; test[0].y=maxy;     test[0].z=tz;
					test[1].x = minx;   test[1].y=maxy+1;   test[1].z=tz;
					Numtests=2;
				}
				else if(tx==minx && tz==minz && is3Dx && is3Dz)
				{
                    test[0].x = minx-1; test[0].y=ty;     test[0].z=minz;
					test[1].x = minx;   test[1].y=ty;     test[1].z=minz-1;
					Numtests=2;
				}
				else if(tx==minx && tz==maxz && is3Dx && is3Dz)
				{
                    test[0].x = minx-1; test[0].y=ty;     test[0].z=maxz;
					test[1].x = minx;   test[1].y=ty;     test[1].z=maxz+1;
					Numtests=2;
				}
				else if(tx==maxx && ty==miny && is3Dx && is3Dy)
				{
                    test[0].x = maxx+1; test[0].y=miny;   test[0].z=tz;
					test[1].x = maxx;   test[1].y=miny-1; test[1].z=tz;
					Numtests=2;
				}
				else if(tx==maxx && ty==maxy && is3Dx && is3Dy)
				{
                    test[0].x = maxx+1; test[0].y=maxy;     test[0].z=tz;
					test[1].x = maxx;   test[1].y=maxy+1;   test[1].z=tz;
					Numtests=2;
				}
				else if(tx==maxx && tz==minz && is3Dx && is3Dz)
				{
                    test[0].x = maxx+1; test[0].y=ty;       test[0].z=minz;
					test[1].x = maxx;   test[1].y=ty;       test[1].z=minz-1;
					Numtests=2;
				}
				else if(tx==maxx && tz==maxz && is3Dx && is3Dz)
				{
                    test[0].x = maxx+1; test[0].y=ty;       test[0].z=maxz;
					test[1].x = maxx;   test[1].y=ty;       test[1].z=maxz+1;
					Numtests=2;
				}
				else if(ty==miny && tz==minz && is3Dy && is3Dz)
				{
                    test[0].x = tx;   test[0].y=miny-1;     test[0].z=minz;
					test[1].x = tx;   test[1].y=miny;       test[1].z=minz-1;
					Numtests=2;
				}
				else if(ty==miny && tz==maxz && is3Dy && is3Dz)
				{
                    test[0].x = tx;   test[0].y=miny-1;     test[0].z=maxz;
					test[1].x = tx;   test[1].y=miny;       test[1].z=maxz+1;
					Numtests=2;
				}
				else if(ty==maxy && tz==minz && is3Dy && is3Dz)
				{
                    test[0].x = tx;   test[0].y=maxy+1;     test[0].z=minz;
					test[1].x = tx;   test[1].y=maxy;       test[1].z=minz-1;
					Numtests=2;
				}
				else if(ty==maxy && tz==maxz && is3Dy && is3Dz)
				{
                    test[0].x = tx;   test[0].y=maxy+1;     test[0].z=maxz;
					test[1].x = tx;   test[1].y=maxy;       test[1].z=maxz+1;
					Numtests=2;
				}
				else if(tx==minx && is3Dx)
				{
					test[0].x = minx-1; test[0].y=ty;   test[0].z=tz;
					Numtests=1;
				}
				else if(tx==maxx && is3Dx)
				{
					test[0].x = maxx+1; test[0].y=ty;   test[0].z=tz;
					Numtests=1;
				}
				else if(ty==miny && is3Dy)
				{
					test[0].x = tx; test[0].y=miny-1;   test[0].z=tz;
					Numtests=1;
				}
				else if(ty==maxy && is3Dy)
				{
					test[0].x = tx; test[0].y=maxy+1;   test[0].z=tz;
					Numtests=1;
				}
				else if(tz==minz && is3Dz)
				{
					test[0].x = tx; test[0].y=ty;   test[0].z=minz-1;
					Numtests=1;
				}
				else if(tz==maxz && is3Dz)
				{
					test[0].x = tx; test[0].y=ty;   test[0].z=maxz+1;
					Numtests=1;
				}
				else
				{
					test[0].x=tx; test[0].y=ty;     test[0].z=tz;
					Numtests=1;
				}
		    }
		  else
		    {
		      Numtests=1;
		      test[0].x = tx; test[0].y = ty; test[0].z =tz;
		    }
	return Numtests;
}
//-------------------------------------------------------------------
void CellularAutomaton::UpdateCells()
//! Call UpdateCell for each member of ActiveCells.
{
	clock_t mystart = clock();

	//Update all the cells that contain, erm, cells
	int NumActiveCells = ActiveCells.size();
	Location pos;

	if (VERBOSE) printf("Starting CellularAutomaton::UpdateCells\n");

	for(int i=0;i<NumActiveCells;i++)
	{
		pos = ActiveCells[i];
		if(Cells[pos.x][pos.y][pos.z][pos.stack].Type == NONE)
		printf("Trying to update NONE, (%d,%d,%d,%d)\n",pos.x,pos.y,pos.z,pos.stack);
		Cells[pos.x][pos.y][pos.z][pos.stack].UpdateCell();
	}

	clock_t mystop = clock();
	TimeUpdateCells += (mystop-mystart)/(double)CLOCKS_PER_SEC;
	if (VERBOSE) printf("@@@@@@ TIMING @@@@@@: Time spent in UpdateCells = %g seconds\n",TimeUpdateCells);
}

//-------------------------------------------------------------------
//! Update the Vasculature.
//! Allows updating with reference to CA properties such as VEGF.
//!
//! @todo Check the effect of newflow true/false (choice between CalculateVascularFlows_Nodal and CalculateVascularFlows_NC)
void CellularAutomaton::UpdateVasculature()
{
	double MaxRadiusChange=param_RadiusChangeTolerance+1;
	bool exceededMaxHaematocritIterations=false;
	bool tempChangedVessels;
	int iterations=0;
	clock_t mystart;
	clock_t mystop;

	Vasc.UpdateMacBloodConcentration();
	Vasc.UpdateMonoBloodConcentration();
	Vasc.UpdateProdrugBloodConcentration();
	Vasc.UpdateDrugBloodConcentration();

	// returns true if there are new vessels
	mystart = clock();
 	ExtravasationToCA();
	mystop = clock();
	if (VERBOSE) printf("FINISHED Vasculature::ExtravasationToCA, elapsed time = %g seconds\n",(mystop-mystart)/(double)CLOCKS_PER_SEC);

	// returns true if there are new vessels
	mystart = clock();
 	tempChangedVessels = Anastomosis();
 	if (tempChangedVessels) Vasc.ChangedVesselsThisStep=true;
	mystop = clock();
	if (VERBOSE) printf("FINISHED Vasculature::Anastomosis, elapsed time = %g seconds\n",(mystop-mystart)/(double)CLOCKS_PER_SEC);

	CalculateVascularDiffusibles();

	Vasc.UpdateVesselTWSSsubthreshold();
	Vasc.AgeVessels();

	tempChangedVessels = RemoveLowFlowVessels();
	if (tempChangedVessels) Vasc.ChangedVesselsThisStep=true;

	// Code for Embolisation therapy
	for(int i =1; i<param_Embolisationtimes.size(); ++i)
	{
		if(GLOBALTIMESTEP == param_Embolisationtimes[i])
		{
			printf("%d. Embolisation Intervall at time %d with radius %.3f\n",i, param_Embolisationtimes[i], param_Embolisationradii[i]);
			Vasc.EmbolisationTherapy(50, 50, 0, param_Embolisationradii[i]); // Holger !!!!! change hard coded position...
		};
		tempChangedVessels = RemoveVesselsCollapsedByEmbolisationTherapy();
		if (tempChangedVessels) Vasc.ChangedVesselsThisStep=true;
	}

	if(param_isRemoveVesselsCollapsedByProliferatingCells)
	{
		tempChangedVessels = RemoveVesselsCollapsedByProliferatingCells();
		if (tempChangedVessels) Vasc.ChangedVesselsThisStep=true;
	}

	tempChangedVessels = RemoveCooptedVessels();
 	if (tempChangedVessels) Vasc.ChangedVesselsThisStep=true;
	if (Vasc.ChangedVesselsThisStep) Vasc.CheckNodeConnectivity(taillocs);

	Vasc.ChangedVesselsThisStep = false;
	UpdateVesselCooptionStates();

	if (VERBOSE) printf("Starting Vasculature flow calculation\n");
	clock_t t1=clock();

	while(MaxRadiusChange>param_RadiusChangeTolerance && iterations<param_MaxVascularIterations)
	{
		if (VERBOSE) printf("Vascular Iteration %d, MaxRadiusChange=%g\n", iterations, MaxRadiusChange);
		Vasc.UpdateVesselResistances();

		if (iterations==0) Vasc.CalculateVascularFlows_NCFirst();
		else Vasc.CalculateVascularFlows_NC();

		Vasc.UpdateVesselWallStress();

		exceededMaxHaematocritIterations = Vasc.CalculateVesselHaematocrit();

		MaxRadiusChange = CalculateVascularRadii();
		iterations++;
	}
	Vasc.FreeSuperLUMemory();

	clock_t t2=clock();
	if (VERBOSE) {
		printf("%.4lf seconds of flow processing\n", (t2-t1)/(double)CLOCKS_PER_SEC);
		printf("Time spent in UpdateVesselResistances = %g seconds\n",Vasc.UVRtime/(double)CLOCKS_PER_SEC);
		printf("Time spent in CalculateVascularFlows = %g seconds\n",Vasc.CVFtime/(double)CLOCKS_PER_SEC);
		printf("Time spent in UpdateVesselWallStress = %g seconds\n",Vasc.UVWtime/(double)CLOCKS_PER_SEC);
		printf("Time spent in CalculateVesselHaematocrit = %g seconds\n",Vasc.CVHtime/(double)CLOCKS_PER_SEC);
		printf("Time spent in CalculateVesselRadii = %g seconds\n",Vasc.CVRtime/(double)CLOCKS_PER_SEC);
	}
	if(exceededMaxHaematocritIterations) printf("WARNING: Exceeded MaxHaematocritIterations\n");
	if(iterations==param_MaxVascularIterations) printf("WARNING: Exceeded MaxVascularIterations\n");
	else if(VERBOSE) printf("Vasculature converged in %d iterations\n", iterations);
}

//-------------------------------------------------------------------
//! Update the State of each vessel to:
//! NORMALVESSEL, ANGIOVESSEL or COLLAPSEVESSEL,
//! depending on the surrounding cell types (NORMALCELLs or CANCERs)
//! and the VEGF level
void CellularAutomaton::UpdateVesselCooptionStates()
{
	int i,k;
	int NumActiveVessels = Vasc.ActiveVessels.size();
	int newstate;

	for(k=0;k<NumActiveVessels;k++)
	{
		i = Vasc.ActiveVessels[k];
		newstate = NORMALVESSEL;
		if (CalculateVesselCooption(i))
		{
			if(Vasc.Vessels[i].MeanVEGF<param_cooptionVc)
			{
				// Collapse
				newstate = COLLAPSEVESSEL;
			}
			else
			{
				// ``Angiogenesis''
				newstate = ANGIOVESSEL;
			}
		}
		if (VERBOSE)
			if(Vasc.Vessels[i].State != newstate)
				printf("Vessel %d state change: %d -> %d\n",i,Vasc.Vessels[i].State,newstate);
		Vasc.Vessels[i].State = newstate;
	}
}

//-------------------------------------------------------------------
//! @brief LoadParameters for each Diffusible
void CellularAutomaton::LoadDiffusibleParameters(const string parfile)
{
	Oxygen.LoadParameters(parfile.c_str());
	VEGF.LoadParameters(parfile.c_str());
	ProDrug.LoadParameters(parfile.c_str());
	Drug.LoadParameters(parfile.c_str());
}

//-------------------------------------------------------------------
//! @brief Run Update for all Diffusibles and propagate to Cells where necessary.
void CellularAutomaton::UpdateDiffusibles()
{
	if (param_isOxygenDiffusible){SetOxygenSourceTerms(); Oxygen.Update(); SetCellOxygens();}
	if (param_isVEGFDiffusible){SetVEGFSourceTerms(); VEGF.Update(); SetExternalVEGF();}
	if (param_isProDrugDiffusible){SetProDrugSourceTerms(); ProDrug.Update();}
	if (param_isDrugDiffusible){SetDrugSourceTerms(); Drug.Update(); SetCellIntercalated();}
}
//-------------------------------------------------------------------
//! Checks elements surrounding each living cell for similar/different
//! living cells and adjusts the thresholds (for p53 induced apoptosis,
//! quiescence entry etc.) for the cells.
//! Should be called before CellDeath function.
void CellularAutomaton::UpdateThresholds()
{
	int NumActiveCells = ActiveCells.size();
	Location pos;

	if (VERBOSE) printf("Starting CellularAutomaton::UpdateThresholds\n");

	for(int i=0;i<NumActiveCells;i++)
	{
		pos = ActiveCells[i];
		int x = pos.x;
		int y = pos.y;
		int z = pos.z;
		int stack = pos.stack;

		int numsame = 0;
		int outof = 0;
		int k;

		if(STENCIL==LOCAL)
		{
			// loop round the test locations
			for(k=0;k<StackCells;k++)
			{
				if(Cells[x][y][z][k].Type==Cells[x][y][z][stack].Type) numsame++;
				if((Cells[x][y][z][k].Type==Cells[x][y][z][stack].Type) || (Cells[x][y][z][k].Type==CANCER)) outof++;
			}

			if((double)numsame/(double)outof > param_UseHigherThresholdsRatio)
				Cells[x][y][z][stack].p53Threshold=param_p53ThresholdHigh[Cells[x][y][z][stack].Type];
			else
			{
				Cells[x][y][z][stack].p53Threshold=param_p53ThresholdLow[Cells[x][y][z][stack].Type];
				if (VERBOSE) printf("LOCAL: Setting threshold low: %d, %d, %f, %f, %f\n",numsame,outof,(double)numsame/(double)outof,param_UseHigherThresholdsRatio,Cells[x][y][z][stack].p53Threshold);
			}
		}
		if(STENCIL==TRANSVESSEL || outof==1)
		{
			outof = 0;
			numsame = 0;
			int minxc = MAX(x-param_DivisionRadius, 0);
			int maxxc = MIN(x+param_DivisionRadius, XCells-1);
			int minyc = MAX(y-param_DivisionRadius, 0);
			int maxyc = MIN(y+param_DivisionRadius, YCells-1);
			int minzc = MAX(z-param_DivisionRadius, 0);
			int maxzc = MIN(z+param_DivisionRadius, ZCells-1);

			if(PERIODIC_X)
			{
			  minxc = x-param_DivisionRadius;
			  maxxc = x+param_DivisionRadius;
			}
			if(PERIODIC_Y)
			{
			  minyc = y-param_DivisionRadius;
			  maxyc = y+param_DivisionRadius;
			}
			if(PERIODIC_Z)
			{
			  minzc = z-param_DivisionRadius;
			  maxzc = z+param_DivisionRadius;
			}

			// Do we consider a 2D or 3D domain?
			bool is3Dx = true;
			bool is3Dy = true;
			bool is3Dz = true;
			if (minxc==maxxc) is3Dx=false;
			if (minyc==maxyc) is3Dy=false;
			if (minzc==maxzc) is3Dz=false;

			Location test[3];
			int Numtests;
			bool found = false;
			int tx,ty,tz;

			for(int txc=minxc;txc<=maxxc;txc++)
			{
				for(int tyc=minyc;tyc<=maxyc;tyc++)
				{
					for(int tzc=minzc;tzc<=maxzc;tzc++)
					{
						tx=txc;
						ty=tyc;
						tz=tzc;

						if(PERIODIC_X) { tx=MOD(txc,(XCells));}
						if(PERIODIC_Y) { ty=MOD(tyc,(YCells));}
						if(PERIODIC_Z) { tz=MOD(tzc,(ZCells));}

						// if test location is a vessel, up to three alternative test
						// locations are used, hence Numtests is set to 1, 2 or 3
						Numtests = TransVesselLocations(x, y, z, tx, ty, tz, test);

						for(int t=0;t<Numtests;t++)
						{
							if(PERIODIC_X) { test[t].x=MOD(test[t].x,(XCells));}
							if(PERIODIC_Y) { test[t].y=MOD(test[t].y,(YCells));}
							if(PERIODIC_Z) { test[t].z=MOD(test[t].z,(ZCells));}

							for(k=0;k<StackCells;k++)
							{
								if(test[t].x>=0 && (test[t].x<XCells || !is3Dx) && test[t].y>=0 && (test[t].y<YCells || !is3Dy) && test[t].z>=0 && (test[t].z<ZCells || !is3Dz) && !((test[t].x==x || !is3Dx) && (test[t].y==y || !is3Dy) && (test[t].z==z|| !is3Dz) && k==stack))
								{
									if(Cells[test[t].x][test[t].y][test[t].z][k].Type==Cells[x][y][z][stack].Type) numsame++;
									if((Cells[test[t].x][test[t].y][test[t].z][k].Type==Cells[x][y][z][stack].Type) || (Cells[test[t].x][test[t].y][test[t].z][k].Type==CANCER)) outof++;
								}
							}
						}

					}
				}
			}
			if(outof==0 || (double)numsame/(double)outof > param_UseHigherThresholdsRatio)
			{
				Cells[x][y][z][stack].p53Threshold=param_p53ThresholdHigh[Cells[x][y][z][stack].Type];
			}
			else
			{
				Cells[x][y][z][stack].p53Threshold=param_p53ThresholdLow[Cells[x][y][z][stack].Type];
			}
		}
	}
}

//-------------------------------------------------------------------
void CellularAutomaton::UpdatetotalMass()
//! Update the totalMass property for each member of ActiveCells.
//! This is to allow cell-cycle progression to depend on the mass at
//! each lattice site. Not used yet.
//! @todo - more efficiently would not revisit an (x,y) location more than once
{
	//Update all the cells that contain, erm, cells
	int NumActiveCells = ActiveCells.size();
	Location pos;

	if (VERBOSE) printf("Starting CellularAutomaton::UpdatetotalMass\n");

	for(int i=0;i<NumActiveCells;i++)
	{
		pos = ActiveCells[i];
		double totalMass=0;
		for (int k=0; k<StackCells; k++)
		{
			totalMass+=Cells[pos.x][pos.y][pos.z][k].GetMass();
		}
		Cells[pos.x][pos.y][pos.z][pos.stack].totalMass = totalMass;
	}
}

//-------------------------------------------------------------------
void CellularAutomaton::WriteOutput(const char* outstem, int i)
{
	char filename[200];
	bool appendsubcellular = false;

	sprintf(filename, "./%s/subcellular/frame%04d.sub", outstem, i);
	SubCellularToFile(filename, appendsubcellular);
	sprintf(filename, "./%s/cells/type/frame%04d.ca", outstem, i);
	CellTypeToFile(filename);

	if (param_isOxygenDiffusible) Oxygen.DumpToFile(outstem, i);
	if (param_isVEGFDiffusible) VEGF.DumpToFile(outstem, i);
	if (param_isDrugDiffusible) Drug.DumpToFile(outstem, i);
	if (param_isProDrugDiffusible) ProDrug.DumpToFile(outstem, i);

	sprintf(filename, "./%s/sum/frame%04d.sum", outstem, i);
	sumToFile(filename);


	sprintf(filename, "./%s/randstate/frame%04d.rnd", outstem, i);
	ofstream stateOut( filename );
	if( stateOut )
	{
		stateOut << mtrand;
		stateOut.close();
	}

}

//-------------------------------------------------------------------
//! Write summary data, produced by CalculateStatistics, to
//! output data directory.
//!
//! @todo Currently clears the file if i=0, this should be amended to
//! allow appending subsequent statistics, e.g. from a restart.
void CellularAutomaton::WriteSummary(const char* outstem, int i)
{
	char filename[200];

	ofstream summary;
	sprintf(filename, "./%s/summary.txt", outstem);
	if (i==0)
		summary.open(filename, ios::trunc);
	else
		summary.open(filename, ios::app);
	summary << GLOBALTIMESTEP << " " << TotalNormal << " " << TotalCancer << " " << TotalQuiescent;
	summary << " " << TotalMacs << " " << TotalMonos << " " << TotalVessel;
	summary << " " << "0" << " " << "0" << " " << TotalEmpty;
	summary << " " << MeanDrug << " " << MeanOxygen << " " << MeanProDrug << " " << MeanVEGF;
	summary << "\n";
	summary.close();
}

//-------------------------------------------------------------------
//! Load initial wound configuration from file.
//! Initial wound configuration is specified by an input file
//! defining rectangular wound region.
//! @param filename a constant character argument
void CellularAutomaton::Wound(string filename)
{
	FILE* fp = fopen(filename.c_str(), "rt");

	if(!fp)
	// input file is bad so print error, try looking in input/, otherwise exit
	{
		printf("WARNING: Error opening wound file \"%s\" ...\n", filename.c_str());
		filename = "input/" + filename;
		printf("   ... trying \"%s\".\n",filename.c_str());
		fp = fopen(filename.c_str(), "rt");
		if(!fp)
		{
			printf("ERROR: Error opening wound file \"%s\" ...\n",filename.c_str());
			exit(1);
		}
	}

	char line[500];
	char type[50];
	char nums[100];
	int startx, endx, starty, endy, startz, endz;

    cout << "************************************************************\n";
	cout << "*** Reading wound file \"" << filename << "\"\n";
	cout << "************************************************************\n";

	while( fgets(&line[0], 500, fp) )
	{
		startx=-1; endx=-1; starty=-1; endy=-1; startz=-1; endz=-1;

		int i=0;
		char c=line[0];

		while(c!=0 && c!='\n' && c!=':')
		{
			type[i++] = c;
			c = line[i];
		}
		type[i]=0;
		while(c!='(' && c!=0 && c!='\n') c=line[i++];

		int j=0;
		c=line[i++];
		while(c!=')' && c!=0 && c!='\n')
		{
			if(c!=' ') nums[j++] = c;
			c=line[i++];
		}
		nums[j]=0;
		sscanf(&nums[0], "%d,%d,%d", &startx, &starty, &startz);
		if(VERBOSE) printf("RECTANGULAR WOUND STARTING AT:  %d, %d, %d ", startx, starty, startz);

		while(c!='(' && c!=0 && c!='\n') c=line[i++];

		j=0;
		c=line[i++];
		while(c!=')' && c!=0 && c!='\n')
		{
			if(c!=' ') nums[j++] = c;
			c=line[i++];
		}
		nums[j]=0;
		sscanf(&nums[0], "%d,%d,%d", &endx, &endy, &endz);
		if(VERBOSE) printf("RECTANGULAR WOUND ENDING AT: %d, %d, %d \n", endx, endy, endz);

		if(startx!=-1 && starty!=-1 && startz!=-1)
		{

			if(endx==-1) endx=startx;
			if(endy==-1) endy=starty;
			if(endz==-1) endz=startz;

			for(int x=startx;x<=endx;x++)
			{
				for(int y=starty;y<=endy;y++)
				{
					for(int z=0;z<ZCells;z++)
					{
						for(int stack=0;stack<StackCells;stack++)
						{
							if(Cells[x][y][z][stack].Type!=NONE && Cells[x][y][z][stack].Type!=VESSEL) KillCell(x, y, z, stack);
							else if(Cells[x][y][z][stack].Type==VESSEL)
							{
								while(Cells[x][y][z][stack].Type==VESSEL)
								{
									int VID = Cells[x][y][z][stack].VesselID;
									Vasc.RemoveVessel(VID,WOUND);
									RemoveVesselFromCA(VID, WOUND);
									if (VERBOSE) printf("In Wound removed vessel\n");
								}
								Vasc.ChangedVesselsThisStep = true;
							}
						}
					}
				}
			}
		}
}
fclose(fp);
}
