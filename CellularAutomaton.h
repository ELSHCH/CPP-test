//! @file CellularAutomaton.h
//! @brief Defines the CellularAutomaton class

#ifndef CELLULARAUTOMATON_H
#define CELLULARAUTOMATON_H

#include "CAElement.h"
#include "Diffusible.h"
#include "DiffusibleOxygen.h"
#include "DiffusibleVEGF.h"
#include "DiffusibleDrug.h"
#include "DiffusibleProDrug.h"

#include "DataStructures.h"
#include "Definitions.h"
#include "MatrixWrap.h"
#include "Parameters.h"
#include "Vasculature.h"
#include "MersenneTwister.h"
#include "RandomWrap.h"

#include <errno.h> //!< needed for CreateOutputDirectories
#include <sys/stat.h> //!< needed for CreateOutputDirectories

#define MOD(r,s) (r%s+s*(r%s<0))

//! @brief Cellular Automaton class
//!
//! This class defines a data structure for the cellular automaton layer
//! of the model.  It consists of a three dimensional
//! array of CAElements that contain the cell data, along with functions to
//! make the automaton functional.  It also contains Diffusible objects
//! representing e.g. oxygen and VEGF concentrations over the spatial
//! area of the cellular automaton.
class CellularAutomaton
{
private:
	double TimeUpdateCells; //!< Total time spent on UpdateCells()

	//! @brief Adds the pressure at each node of the vascular system to the CA.
	void AddVesselPressureToCA();

	bool CalculateVesselCooption(int v);

	bool CanDivideDirect4(int x, int y, int z, int stack, Location&);
	bool CanDivideLocal(int x, int y, int z, int stack, Location&);
	bool CanDivideRandom(int x, int y, int z, int stack, Location&);
	bool CanDivideTransVessel(int x, int y, int z, int stack, Location&);
	bool CanMoveRandomMoore(int x, int y, int z, int stack, Location&);
	bool CanMoveRandomNeumann(int x, int y, int z, int stack, Location&);
	
	//! @brief Remove a ProtoVessel from the CellularAutomaton CA.
	void RemoveProtoVesselFromCA(int vessel);
	//! @brief Returns transvessel locations to test for division, thresholds, etc. 
	int TransVesselLocations(int x, int y, int z, int tx, int ty, int tz, Location (&testSites)[3]);

	//! @brief Randomly shuffle the list of Locations and Probabilities.
	void shuffleLocationsAndProbabilities(Location Test[], double Pij[], const int NumTests);
	//! @brief Write position of a single cell to file
	void SingleCellToFile(char* filename, bool append);

public:
	int TotalCancer; //!< Total number of CANCER cells
	int TotalQuiescent; //!< Total number of QUIESCENT CANCER cells
	int TotalNormal; //!< Total number of NORMALCELL cells
	int TotalNormalPasteur; //!< Total number of PASTEUR NORMAL cells
	int TotalCancerWarburg; //!< Total number of WARBURG CANCER cells
    int TotalCancerWarburgQuiescent; //!< Total number of QUIESCENT WARBURG CANCER  cells
	int TotalMacs; //!< Total number of MACROPHAGE cells
	int TotalMonos; //!< Total number of MONOCYTE cells
	int TotalVessel; //!< Total number of VESSEL cells
	int TotalEmpty; //!< Total number of Empty sites
	double MeanDrug; //!< Mean Drug across the tissue region
	double MeanOxygen; //!< Mean Oxygen across the tissue region
	double MeanProDrug; //!< Mean ProDrug across the tissue region
	double MeanVEGF; //!< Mean VEGF across the tissue region

	double TimeRunStep; //!< Total time spent on the whole of RunStep()
	double TimeCellMovement; //!< Total time spent on CellMovement()

//! @brief The random number generator object
//!
//! We use std::rand(), via the interface in RandomWrap.h,
//! or the implementation of the Mersenne Twister in MersenneTwister.h
//!
//! See: http://www-personal.umich.edu/~wagnerr/MersenneTwister.html
//! For other examples see:
//! 	http://www.agner.org/random/
//! 	http://www.codeproject.com/KB/recipes/SimpleRNG.aspx
//! 	http://www.johndcook.com/cpp_TR1_random.html
//! 	http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html
#ifdef STDRAND
	vasctumRand mtrand;
#else
	MTRand mtrand;
#endif
	//! @brief 4-D array of cells at each location.
	//!
	//! The 4-D array (3-space, one stack to allow multiple cells to occupy
	//! a location) of CAElements containing the information about the
	//! type and state of cells at each location.  This is public and may be
	//! addressed as Cells[x][y][z][stack].
	CAElement**** Cells;
	ChosenVasculature Vasc;
	//! @brief X dimension of the cell array in elements
	//!
	//! Declared public. Should not be modified except by internal functions.
	//! @attention What does this mean?
	int XCells;
	//! @brief Y dimension of the cell array in elements
	//!
	//! Declared public. Should not be modified except by internal functions.
	//! @attention What does this mean?
	int YCells;
	//! @brief Z dimension of the cell array in elements
	//!
	//! Declared public. Should not be modified except by internal functions.
	//! @attention What does this mean?
	int ZCells;
	//! @brief Stack dimension of the cell array (x,y,z)
	//!
	//! Declared public. Should not be modified except by internal functions.
	//! @attention What does this mean?
	int StackCells;


	int MaxCellID; //!< The maximum CellID across all elements of Cells.
	//! @brief Array to count number of divisions at each (x,y,z) location
	int*** divsum;
	//! @brief Array to count number of macs extravasating at each (x,y,z) location
	int*** macsum;
	//! @brief Array to count number of monos extravasating at each (x,y,z) location
	int*** monosum;
	//! @brief Array to count number of cells dying at each (x,y,z) location (APOPTOSIS)
	int*** aposum;
	//! @brief Array to count number of cells dying at each (x,y,z) location (NO SPACE TO DIVIDE)
	int*** diesum;
	//! @brief Array to count number of cells killed by drug at each (x,y,z) location
	int*** killsum;
	//! @brief Array to map locations where endothelial sprouts emerge
	int*** taillocs;
	//! @brief Array to store cumulative tip cell sprouting
	int*** tipcelloutsum;
	//! @brief Array to map current locations of endothelial tips
	int*** tiplocs;
	//! @brief Array to store distribution of accumulated drug
	double*** drugsum;
	//! @brief Array to store distribution of accumulated TACE particles
	double*** TACEParticlessum;

	std::vector<double> VascPropSource; //!< Contribution to source terms from vessels
	void SetVascPropSource(int xgrid, int ygrid, int zgrid, double value);		//Set the Constant source term at one of the grid points

	DiffusibleOxygen Oxygen; //!< Oxygen field
	DiffusibleVEGF VEGF; //!< VEGF field
	DiffusibleDrug Drug; //!< Drug field
	DiffusibleProDrug ProDrug; //!< ProDrug field

	//! @brief List of active CAElement locations
	//!
	//! A List (see DataStructures) of the automaton Locations
	//! (see DataStructures) containing living cells.
	//! This list is used in updates and must be maintained by functions
	//! that create/kill cells.
	std::vector<Location> ActiveCells;

	//! @brief CellularAutomaton constructor
	CellularAutomaton();
	//! @brief CellularAutomaton destructor
	~CellularAutomaton();

	//! @brief Add the Diffusible terms at a single Vessel lattice site
	virtual void AddVesselSiteDiffusibleToCA(int vessel, int x, int y, int z);

	//! @brief Add a single Vessel lattice site
	void AddVesselSiteToCA(int vessel, int x, int y, int z);
	//! Add vessels to the underlying CellularAutomaton.
	void AddVesselsToCA();
	//! @brief Add a single vessel (given as an index into the Vessels array), to the CellularAutomaton CA.
	void AddVesselToCA(int vessel);

	//! @brief This function ages all cells by TimeStep
	void AgeCells(double TimeStep);

	//! @brief Find new vessel connections.
	bool Anastomosis();

	//! @brief Calculate summary statistics during or after simulation.
	virtual void CalculateStatistics();

	//! @brief Calculate Diffusible concentrations across the Vasculature
	virtual void CalculateVascularDiffusibles();

	virtual double CalculateVascularRadii();
	
	bool CanDivide(int x, int y, int z, int stack, Location&);
	bool CanMove(int x, int y, int z, int stack, Location&);

	//! Cell cycle model test
	void CellCycleTest(int nsteps, string outstem);

	virtual void CellDeath();
	//! @brief Check the internal state of a cell to see if it should try to divide.
	virtual void CellDivision();
	//! @brief Move all Active cells in the CA.
	virtual void CellMovement();
	void CellQuiescence();


	void CheckForNONEInActiveCells();

	//! Create Output Directories
	virtual void CreateOutputDirectories(string outstem, string parstem);
	//! Create Output Directories from a list
	void CreateOutputDirectoriesFromList(string outstem, string parstem, std::vector<string> subdirs);

	void DEBUG_Print();
	
	virtual void DeleteMemory();
	
	double ExtravasationConditionTIPCELL(int v, int curx, int cury, int curz);

	virtual void ExtravasationFromVesselSite(int vessel, int Type, Node nodePos);
	//! @brief Allow extravasation of cells from vessels.
	virtual void ExtravasationToCA();

	virtual void Initialise(int NumXCells, int NumYCells, int NumZCells, int NumStackCells);
	void InitialiseSimulation();
	
	//! @brief Empty all Cells in a CellularAutomaton.
	void EmptyCACells();
	//! @brief Load initial CA configuration from file.
	void LoadConfig(string filename);
	virtual void LoadDiffusibleParameters(const string parfile);
	void KillCell(int x, int y, int z, int stack);
	
	//! Main step of simulation, CellularAutomaton and Vasculature
	void MainStep(string outstem);
	
	//! Main step of simulation, CellularAutomaton and Vasculature - reordered
	void MainStepReorder(string outstem);
	
	//! Main step of simulation, CellularAutomaton and Vasculature - consolidated
	void MainStepConsolidated(string outstem);
	
	//! @brief Make cell based on parent and daughter cell locations 
	virtual void MakeCellDaughterFromParent(Location& parent, Location& daughter);
	void MakeCell(int x, int y, int z, int stack, int Type);

	//! @brief Move an individual Cell in the CA. Checks tip cells for regression.
	void MoveCell(Location& OldPos, Location& NewPos);

	//! Random walk test
	void RandomWalkTest(int nsumsteps, string outstem);

	//! @brief Removes vessels that have had WSS too low for too long.
	bool RemoveLowFlowVessels();
	//! @brief Remove vessels that have been coopted.
	bool RemoveCooptedVessels();
	//! @brief Remove vessels collapsed by proliferating cells.
	bool RemoveVesselsCollapsedByProliferatingCells();
	//! @brief Remove vessels collapsed by embolisation therapy.
	bool RemoveVesselsCollapsedByEmbolisationTherapy();
	//! @brief Remove a vessel from the CellularAutomaton CA.
	virtual void RemoveVesselFromCA(int v, int CallingFunction);

	//! Run simulation
	void RunSimulation(int nsumsteps, int nsteps, string outstem);

	//! @brief Run a single step of the Cellular Automaton 
	virtual void RunStep();
	//! @brief Run a single step of the Cellular Automaton - altered ordering of processing 
	virtual void RunStepReorder();
	//! @brief Run a single step of the Cellular Automaton - consolidated order (e.g. all Diffusibles together) 
	virtual void RunStepConsolidated();
	//! @brief Run a single step of the Cellular Automaton with only a cell cycle model. 
	virtual void RunCellCycleStep();

	//! @brief Returns cell count at a site to test for division, movement, etc. 
	int StackCount(int x, int y, int z, Location &emptySites);

	//! @brief Run Update for all Diffusibles and propgate to Cells where necessary.
	virtual void UpdateDiffusibles();

	virtual void UpdateVasculature();
	//! @brief Read cell type data from file
	void CellTypeFromFile(char* filename);
	//! @brief Write cell type data to file
	void CellTypeToFile(char* filename);

	//! @brief Randomly shuffle the list of active cells.
	void shuffleActiveCells();

	void SubCellularFromFile(string filename);
	//! @brief Write subcellular data to file
	void SubCellularToFile(char* filename, bool append);
	//! @brief Output location summaries - divisions, extravasation, etc
	void sumToFile(char* filename);

	//! @brief Update the totalMass property of all living cells (CAElements).
	void UpdatetotalMass();

	//! @brief Run the UpdateCell function of all living cells (CAElements).
	void UpdateCells();
	//! @brief Updates the State of every vessel (i.e. NORMALVESSEL, ANGIOVESSEL or COLLAPSEVESSEL).
	void UpdateVesselCooptionStates();
	//! @brief Adjust internal thresholds based on surrounding cells.
	void UpdateThresholds();

	//! @brief Write output to files.
	virtual void WriteOutput(const char* outstem, int i);
	//! @brief Write summary statistics to files.
	virtual void WriteSummary(const char* outstem, int i);
	//! @brief create a wound from file
	void Wound(string filename);

	//! @brief Set the source terms for Oxygen
	virtual void SetOxygenSourceTerms(); //!< may be overridden for derived models
	//! @brief Set the source terms for Drug
	void SetDrugSourceTerms();
	//! @brief Set the source terms for ProDrug
	void SetProDrugSourceTerms();
	//! @brief Set the source terms for VEGF
	void SetVEGFSourceTerms();

	//! @brief Set the internal OxygenLevel according to the Diffusible distribution.
	void SetCellOxygens();
	//! @brief Set the ExternalVEGF according to the Diffusible distribution.
	void SetExternalVEGF();
	//! @brief Set the internal Intercalated switch according to the drug distribution.
	void SetCellIntercalated();
	//! @brief Set the internal DrugLevel according to the Diffusible distribution.
	void SetCellDrug();
};

#endif
