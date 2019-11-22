#pragma once

#include "Common.h"
#include "Grid.h"
//#include "DCPQuery.h"
#include "Particles.h"

#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <dirent.h>

//////////////////////////////////////////////////////////////////////////
////Particle fluid simulator
template<int d> class ExplosionSim
{
	using VectorD = Vector<real, d>; using VectorDi = Vector<int, d>;
public:
	////// Data Structure Variables
	Grid<d> grid;
	Array<VectorD> velocities;							//Velocity    on grid cells
	Array<real>    densities;							//Density     on grid cells
	Array<real>	   temps;								//Temperature on grid cells
	Array<real>    pressures;							//Pressure    on grid cells
	Array<real>    div_vels;							//Divergence  on grid cells
	Array<real>	   vorticities;							//Vortices    on grid cells

	Array<Array<real>> controlPaths;

	////// Other Constants and Variables

	// Grid and Control Paths
	int node_num;
	int n_per_dim = 32;
	real w = 3;							//Width of sweeping region- m (AC)

	// Physics
	const real gamma = 1.4;				//Ratio of specific heats at ambient pressure (no units)
	const real temp_amb = 25;			//Ambient temperature- Celcius
	const real p_0 = 1;					//Ambient pressure at sea level- atm
	real explosionTemp = 1000000;		//Temperature at start of explosion
	real P_p = 30;						//Peak overpressure- atm (AC- arbitarily chosen)
	real P_m = -10;						//Minimum negative pressure (AC)
	real t_d = 7 / 100000;				//Time to reach p_0 (AC)
	real t_n = 3 * t_d;					//Time to reach p_0 from t_d (AC)
	real b = 4;							//Decreasing coefficient (AC)
	real density_0 = 1.1839;			//Ambient density at STP (AC)
	const real p_scale = 0.5;        // scale the pressure by pressure velocity with this constant (randomly picked)

	//time
	real time;							//Current time passed (seconds)
	real dt = 0.01;						//Timestep (seconds)
	const int dim = 3;					//Number of dimensions

	// Simulation
	bool explosion_done;				//Determines whether or not we are only doing process 5

	virtual void Initialize()
	{
		//Init Grid
		VectorDi cell_counts = VectorDi::Ones() * n_per_dim;
		real dx = (real)1. / (real)n_per_dim;
		VectorD domain_min = VectorD::Zero();
		grid.Initialize(cell_counts, dx, domain_min);

		//Resize arrays for all grid nodes
		node_num = grid.node_counts.prod();
		velocities.resize(node_num, VectorD::Zero());
		densities.resize(node_num, density_0);
		temps.resize(node_num, temp_amb);
		pressures.resize(node_num, p_0);
		div_vels.resize(node_num, 0);
		vorticities.resize(node_num, 0);

		//TODO: Initiating particles
	}

	virtual void Advection(real dt)
	{
	}

	virtual void Projection()
	{
	}

	virtual void Set_Boundary_Conditions()
	{
	}

	virtual void Advance(const real dt)
	{
		time += dt;
		if (!explosion_done)
		{
			//TODO: Calculate Density Opacity
			//TODO: Get nurbs curve value - getDensityFront, getPressureFront
			//Get which cells need to be looked at rn - Sweep Region
			//TODO: Density and Temperature in Sweep Region
			//TODO: Make fuel particles
			//TODO: Calculate propagation velocity magnitude of density
			//TODO: Create vortex particles
			//TODO: Amplify pressure fields by detonation-y stuff
			//TODO: Vortical effects
		}
		Advection(dt);
		Projection();
	}

	virtual void SweepRegion(const VectorD& pos, Array<int>& cells){

		// pos is center of circle
		// w is radius of circle
		// TODO: change distance from sphere to circle

		for (int i = 0; i < cells.size(); i++) {

			real distance = 0;
			Vector3i currentCellCoordinates = Coord(i);
			for (int j = 0; j < d; j++) {

				distance += sqrt(pos[j] * pos[j] - currentCellCoordinates[j] * currentCellCoordinates[j]);
			}

			if (distance < w) {

				// Allocate uniform density value to each grid square here based on densityOpacityCurve
				densities[i] = densityOpacityCurve(time);

				// Allocate user-specified uniform temperature value, based upon how much time has passed, to each grid square
				temps[i] = getCurrentExplosionTemp();

				// Allocate direction velocity to each grid square according to u_d(G(g)) = V(t_i)t(g), where g is a grid square,
				// G(g) is the set of grid square in the region, and t(g) is the unit tangent vector derived from g on the flow control path
				// V(t_i) is the same for each grid square in the region

				// TODO: calculate the following variable
				Vector3i unitVectorTangentToControlPath;
				velocities[i] = densityPropagationCurve(time, temps[i]) * unitVectorTangentToControlPath;
			}
		}
	}

	virtual void PreProcessing(){//string dir_name){
		// 	Intialize Variables
		Array<std::string> file_names;

		// open through directory, find array of text files
		read_directory("/Users/student/explosions/NURBS", file_names);

		//for every textfile, call read_from_file, add to array
		for(int i = 0; i<file_names.size(); i++){
			controlPaths.push_back(read_from_file(file_names[i]));
		}
		std::ofstream outfile ("test.txt");

		for(int i = 0; i<controlPaths.size(); i++){
			outfile<<"-------------------------"<<std::endl;
			for(int j = 0; j<controlPaths[i].size(); j++){
				outfile<<string(controlPaths[i][j])<<std::endl;
			}
		}
		outfile.close();


	}


/////////////////////////////////// HELPER FUNCTIONS ///////////////////////////////////
protected:
	//get front of curve point

	/////////////////////////////////// Control Path Functions ///////////////////////////////////

	//get the front of the density wave
	 real getDensityFront(const real t, const real dt, const real temperature, const Array<real> path){
		// variables
		real curTime = 0;
		real curDistance = 0;
		real sumDistance = 0;
		real curIndex = 1;
		// sum the distance along the path
		while(curTime<t){
			sumDistance += densityPropagationCurve(curTime, temperature)*dt;
			curTime += dt;
		}
		// find the point along the path corresponding to that distance
		while(curDistance<sumDistance && (curIndex)<path.size()){
			curDistance += distanceNd(path[curIndex], path[curIndex-1]);
			curIndex++;
		}
		return path[curIndex];
	}

	// get the front of the pressure wave
	real getPressureFront(const real t, const real dt, const real temperature, const Array<real> path){
	 // variables
	 real curTime = 0;
	 real curDistance = 0;
	 real sumDistance = 0;
	 real curIndex = 1;
	 // sum the distance along the path
	 while(curTime<t){
		 sumDistance += pressurePropagationCurve(curTime, temperature)*dt;
		 curTime += dt;
	 }
	 // find the point along the path corresponding to that distance
	 while(curDistance<sumDistance && (curIndex)<path.size()){
		 curDistance += distanceNd(path[curIndex], path[curIndex-1]);
		 curIndex++;
	 }
	 return path[curIndex];
 }

	inline real scale_by_distance(const real value, const real dist_trav, const real total_length){
		return (value*total_length) / (dist_trav);
	}

	//find the length of a control path
	real find_path_length(Array<real> path){
		real pathLength = 0.0;
		for(int i = 1; i<path.size(); i++){
			pathLength += distanceNd(path[i], path[i-1]);
		}
		return pathLength;
	}

/////////////////////////////////// Reading Directories and Files ///////////////////////////////////

	void read_directory(const std::string& name, Array<std::string> v)
	{
			DIR* dirp = opendir(name.c_str());
			struct dirent *dp;
			while ((dp = readdir(dirp)) != NULL) {
				size_t len = strlen(dp->d_name);
				if(len > 4 && strcmp(dp->d_name + len - 4, ".txt") == 0){
						v.push_back(dp->d_name);
				}
			}
			closedir(dirp);
	}

	//read points from file and create control paths
	Array<real> read_from_file(std::string file_path){ //what params?
		std::ifstream inFile;

		Array<real> grid_indices;

		//read from file given file path
		inFile.open(file_path);
		std::string data;

		//read from file, get string
		if (!inFile) {
			 std::cout << "Unable to open file";
			 exit(1); // terminate with error
	 	}
		for(std::string line; std::getline(inFile, line); ) {
			Array<real> pointReals;
			Vector3 point;

			std::istringstream lineIterator(line);
			Array<std::string> parsedLine(std::istream_iterator<std::string>{lineIterator},
				std::istream_iterator<std::string>());

			for(int i = 0; i<parsedLine.size(); i++){
				pointReals.push_back((real)parsedLine[i]);
			}
			point = Vector3(pointReals.data());

			VectorDi coordinates = grid.Cell_Coord(point);

			real index = Idx(coordinates);

			grid_indices.push_back(index);

		}
		//process string, get array of points

		// turn array of points into grid points

		// return array

	}

	/////////////////////////////////// Physics Helper Functions ///////////////////////////////////

	// Returns pressure based on current time
	real pressureMagnitudeCurve(real t) {

		real output = p_0;
		if (t <= t_d) {
			output += P_p * (1 - t / t_d) * exp((-b*t) / t_d);
		}
		else if (t <= t_d + t_n / 2) {
			output -= (2 * (P_m / t_n) * (t - t_d));
		}
		else if (t <= t_d + t_n) {
			output -= 2 * (P_m / t_n) * (t_d + t_n - t);
		}
		else {
			// do nothing, output should just be p_0
		}

		return output;
	}

	// Returns density based on current time
	real densityOpacityCurve(real t) {

		real output = (gamma + 1) * pressureMagnitudeCurve(t) + 2 * gamma * p_0;
		output /= ((gamma - 1) * pressureMagnitudeCurve(t) + 2 * gamma * p_0);
		output *= density_0;

		return output;
	}

	// Returns velocity magnitude v_p(t)
	real pressurePropagationCurve(real t, real temperature) {

		real c = getSpeedOfSoundInAir(temperature);
		return c * sqrt(1 + ((gamma + 1) * (pressureMagnitudeCurve(t) / (2 * gamma * p_0))));
	}

	// returns velocity magnitude v(t)
	real densityPropagationCurve(real t, real temperature) {

		real c = getSpeedOfSoundInAir(temperature);
		real output = 1 / sqrt(1 + (gamma + 1) * pressureMagnitudeCurve(t) / (2 * gamma * p_0));
		output *= c;
		output *= pressureMagnitudeCurve(t);
		output /= gamma;
		output /= p_0;

		return output;
	}

	// Returns speed of sound for a given temperature in m/s
	real getSpeedOfSoundInAir(real temperature) {
		return 331.5 + 0.6 * temperature;
	}

	// Returns uniform temperature value based on time elapsed
	real getCurrentExplosionTemp() {
		return explosionTemp - 100000 * time;
	}

/////////////////////////////////// Interpolation and Distance ///////////////////////////////////

	VectorD Interpolate(const Array<VectorD>& u, VectorD& pos)
	{
		//// define variables
		const real one_over_dx = 1.0/grid.dx;
		////clamp pos
		for (int i = 0; i < d; i++) {
			if (pos[i] < grid.domain_min[i])pos[i] = grid.domain_min[i];
			else if (pos[i] > grid.domain_max[i])pos[i] = grid.domain_max[i];
		}

		VectorDi cell = ((pos - grid.domain_min) * one_over_dx).template cast<int>();
		VectorD frac = (pos * one_over_dx - cell.template cast<real>());
		return Interpolate_Helper(cell, frac, u);
	}

	////2D bi-linear interpolation
	Vector2 Interpolate_Helper(const Vector2i& cell, const Vector2& frac, const Array<Vector2>& u)
	{
		return ((real)1 - frac[0]) * ((real)1 - frac[1]) * u[grid.Node_Index(cell)]
			+ frac[0] * ((real)1 - frac[1]) * u[grid.Node_Index(Vector2i(cell[0] + 1, cell[1]))]
			+ ((real)1 - frac[0]) * frac[1] * u[grid.Node_Index(Vector2i(cell[0], cell[1] + 1))]
			+ frac[0] * frac[1] * u[grid.Node_Index(Vector2i(cell[0] + 1, cell[1] + 1))];
	}

	////3D tri-linear interpolation
	Vector3 Interpolate_Helper(const Vector3i& cell, const Vector3& frac, const Array<Vector3>& u)
	{
		return lerp(
			lerp(
				lerp(u[Idx(cell)],
					u[Idx(Vector3i(cell[0] + 1, cell[1], cell[2]))],
					frac[0]),
				lerp(u[Idx(Vector3i(cell[0], cell[1] + 1, cell[2]))],
					u[Idx(Vector3i(cell[0] + 1, cell[1] + 1, cell[2]))],
					frac[0]),
				frac[1]),
				lerp(u[Idx(Vector3i(cell[0], cell[1], cell[2] + 1)],
				u[Idx(Vector3i(cell[0] + 1, cell[1], cell[2] + 1))],
				frac[0]),
				lerp(u[Idx(Vector3i(cell[0], cell[1] + 1, cell[2] + 1))],
					u[Idx(Vector3i(cell[0] + 1, cell[1] + 1, cell[2] + 1))],
					frac[0]),
				frac[1]),
			frac[2]);
	}

	template <class T>
	T lerp(T val1, T val2, real v) { return v * val2 + (1 - v) * val1; }


	inline real distanceNd(VectorD& pos1, VectorD& pos2){
		real distance = 0.0;
		for(int i = 0; i<d; i++){
			distance += (pos1[i] - pos2[i]) * (pos1[i] - pos2[i]);
		}
		return sqrt(distance);
	}


/////////////////////////////////// Grid Helper Functions ///////////////////////////////////

	////return the node index given its coordinate
	int Idx(const Vector3i& node_coord) const
	{
		return grid.Node_Index(node_coord);
	}

	////return the coordinate given its index
	VectorDi Coord(const int node_index) const
	{
		return grid.Node_Coord(node_index);
	}

	////return the node position given its index
	VectorD Pos(const int node_index) const
	{
		return grid.Node(node_index);
	}

	////check if a node is on the boundary of the grid
	////given its coordinate or index
	bool Bnd(const Vector3i& node_coord) const
	{
		for (int i = 0; i < d; i++) {
			if (node_coord[i] == 0 || node_coord[i] == grid.node_counts[i] - 1)
				return true;
		}
		return false;
	}
	bool Bnd(const int node_index) const
	{
		return Bnd(Coord(node_index));
	}
};
