#pragma once

#include "Common.h"
#include "Grid.h"
#include "Particles.h"

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <dirent.h>

//////////////////////////////////////////////////////////////////////////
////Particle fluid simulator
template <int d>
class ExplosionSim
{
	using VectorD = Vector<real, d>;
	using VectorDi = Vector<int, d>;

public:
	////// Data Structure Variables
	Grid<d> grid;
	Array<VectorD> velocities; //Velocity    on grid cells
	Array<real> densities;	 //Density     on grid cells
	Array<real> temps;		   //Temperature on grid cells
	Array<real> pressures;	 //Pressure    on grid cells
	Array<real> div_vels;	  //Divergence  on grid cells
	Array<real> vorticities;   //Vortices    on grid cells

	Array<Array<int>> controlPaths;
	Array<Array<real>> st_times;

	////// Other Constants and Variables

	// Grid and Control Paths
	int node_num;
	int n_per_dim = 32;
	real w = 3; //Width of sweeping region- m (AC)
	int cur_index = 0;

	// Physics
	//based on velocity of 7000 m/s
	// t is scaled by T and length of control path
	// so assign start time using dt?
	const real gamma = 1.4;		  //Ratio of specific heats at ambient pressure (no units)
	const real temp_amb = 25;	 //Ambient temperature- Celcius
	const real p_0 = 1;			  //Ambient pressure at sea level- atm
	real explosionTemp = 1000000; //Temperature at start of explosion
	real P_p = 30;				  //Peak overpressure- atm (AC- arbitarily chosen)
	real P_m = -10;				  //Minimum negative pressure (AC)
	real t_d = 7 / 100000;		  //Time to reach p_0 (AC)
	real t_n = 3 * t_d;			  //Time to reach p_0 from t_d (AC)
	real b = 4;					  //Decreasing coefficient (AC)
	real density_0 = 1.1839;	  //Ambient density at STP (AC)
	const real p_scale = 0.5;	 // scale the pressure by pressure velocity with this constant (randomly picked)
	real q = 500.0;

	//time
	real time;		   //Current time passed (seconds)
	real dt = 0.01;	//Timestep (seconds)
	const int dim = 3; //Number of dimensions

	// Simulation
	bool explosion_done; //Determines whether or not we are only doing process 5

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
		////Advection using the semi-Lagrangian method
		Array<VectorD> u_copy = velocities;
		for (int i = 0; i < node_num; i++) {

			velocities[i] = VectorD::Zero();

			VectorD position = Pos(Coord(i)) - (velocities[i] * dt / 2);
			velocities[i] = Interpolate(u_copy, position);
			VectorD position_1 = Pos(Coord(i)) - velocities[i] * dt;
			velocities[i] = Interpolate(u_copy, position_1);
		}
	}

	virtual void Projection()
	{
		real dx = grid.dx;
		real dx2 = grid.dx * grid.dx;

		////Projection step 1: calculate the velocity divergence on each node
		////Read this sample code to learn how to access data with the node index and coordinate
		std::fill(div_vels.begin(), div_vels.end(), (real)0);
		for (int i = 0; i < node_num; i++) {

			if (Bnd(i))continue;		////ignore the nodes on the boundary
			VectorDi node = Coord(i);
			div_vels[i] = (real)0;

			for (int j = 0; j < d; j++) {

				VectorD u_1 = velocities[Idx(node - VectorDi::Unit(j))];
				VectorD u_2 = velocities[Idx(node + VectorDi::Unit(j))];
				div_vels[i] += (u_2[j] - u_1[j]) / (2 * dx);
			}
		}

		////Projection step 2: solve the Poisson's equation -lap p= div u
		////using the Gauss-Seidel iterations
		std::fill(pressures.begin(), pressures.end(), (real)0);
		for (int iter = 0; iter < 40; iter++) {
			for (int i = 0; i < node_num; i++) {

				if (Bnd(i))continue;		////ignore the nodes on the boundary
				VectorDi node = Coord(i);

				real temp = 0;
				for (int j = 0; j < d; j++) {

					temp += pressures[Idx(node + VectorDi::Unit(j))];
					temp += pressures[Idx(node - VectorDi::Unit(j))];

				}
				pressures[i] = (-div_vels[i] + (1 / dx2) * temp) / (4 / dx2);
			}
		}

		////Projection step 3: correct velocity with the pressure gradient
		for (int i = 0; i < node_num; i++) {
			if (Bnd(i))continue;		////ignore boundary nodes
			VectorDi node = Coord(i);
			VectorD grad_p = VectorD::Zero();

			for (int j = 0; j < d; j++) {
				grad_p[j] = (1 / (2 * dx)) * (pressures[Idx(node + VectorDi::Unit(j))] -
			     pressures[Idx(node - VectorDi::Unit(j))]);
			}
			velocities[i] -= grad_p;
		}
	}

/* NOTE: change to only work on the paths?? Probably*/
	void Vorticity_Confinement(const real dt, Array<Array<int>> sweep_cells)
	{
		real dx=grid.dx;
              ////Vorticity confinement step 1: update vorticity
		std::fill(vorticities.begin(),vorticities.end(),(real)0);
		for(int m = 0; m<sweep_cells.size(); m++)
		{
		  for(int j=0;j<sweep_cells[m].size();j++)
			{
				int i = sweep_cells[m][j];
				if(Bnd(i)){continue;}             ////ignore boundary nodes
		    VectorDi node=Coord(i);

				vorticities[i] = updateVorticity(dt, dx, node);

		  }

			////Vorticity confinement step 2: update N = (grad(|vor|)) / |grad(|vor|)|
			Array<VectorD> N(node_num,VectorD::Zero());
			for(int k=0;k<sweep_cells[m].size();k++)
			{
				int i = sweep_cells[m][k];
				if(Bnd(i)){continue;}             ////ignore boundary nodes
				VectorDi node=Coord(i);
				N[i]=VectorD::Zero();

		    // Calculate divergence of vorticity (same as projection step 1 code)
				for(int j=0;j<d;j++)
				{
		    	real vor_1=vorticities[Idx(node-VectorDi::Unit(j))];
		    	real vor_2=vorticities[Idx(node+VectorDi::Unit(j))];
		    	N[i][j] = (vor_2 - vor_1);
		    }

	 			//Normalize
		    N[i].normalize();

		  }

		  ////Vorticity confinement step 3: calculate confinement force and use it to update velocity
		  real vor_conf_coef=(real)4;
			for(int j=0;j<sweep_cells[m].size();j++)
	 	 	{
	 		  int i = sweep_cells[m][j];
		  	if(Bnd(i)){continue;}             ////ignore boundary nodes
		    VectorD f=vor_conf_coef*dx*Cross(N[i],vorticities[i]);
		    velocities[i]+=f*dt;     ////we don't have mass by assuming density=1
		  }
		}
	}
	inline Vector3 Cross(const Vector3& a,const Vector3& b) const
	{return Vector3( (a[1] * b[2] - a[2] * b[1]), (a[0] * b[2] - a[2] * b[0]), (a[0] * b[1] - a[1] * b[0]));}




	// Probably don't need this- Bnd should do this
	virtual void Set_Boundary_Conditions()
	{
	}

	virtual void Advance(const real dt)
	{
		time += dt;
		cur_index += 1;

		Array<Array<int>> sweepRegions;

		if (!explosion_done)
		{
			//for each control path
			for(int i = 0; i<controlPaths.size(); i++)
			{
				//TODO: Calculate Density Opacity

				//for each node behind the current one
				///////// FIX THIS LOGIC

					//find relative time to front
					int grid_cell = controlPaths[i][cur_index];
					Array<int> sweepRegion;
					SweepRegion(grid_cell, sweepRegion, controlPaths[i]);
					sweepRegions.push_back(sweepRegion);
					//add pressure here
					for(int node = 0; node<controlPaths[i].size(); node++)
					{
				}
			}

		}

			//Get which cells need to be looked at rn - Sweep Region
			//TODO:- Density and Temperature in Sweep Region
			//TODO: Make fuel particles
			//TODO:- Calculate propagation velocity magnitude of density
			//TODO:- Create vortex particles
			//TODO:- Amplify pressure fields by detonation-y stuff
			//TODO: Vortical effects

		Vorticity_Confinement(dt, sweepRegions);
		Advection(dt);
		Projection();
	}

	virtual void drasticPressureChange(real t, real temperature) {
		real M = pressurePropagationCurve(t, temperature) / getSpeedOfSoundInAir(temperature);
		real a = -(M * M * gamma);
		real b = 1 + 2 * M * M * gamma + a;
		real c = M * M * gamma - a * a - (4 * gamma) / (gamma + 1) - 2 * q * a;

		// n is the change rate for 1 / density
		real n = (8 * M * M * M * M * gamma * gamma) + (4 * M * M * gamma * gamma * a * a * gamma) -
			(16 * M * M * gamma * gamma) / (gamma + 1) - (8 * M * M * gamma * a * q) +
			(4 * M * M * gamma) + (a * a) + (2 * a) + 1;
	}

	virtual void SweepRegion(int index, Array<int> &cells, Array<int> &path)
	{
		//eq_pos_time = start time for index. with the current time step it should be 0

		// use dot product to scale time values.
		VectorD pos = Pos(Coord(index));
		for (int i = 0; i < cells.size(); i++)

			real distance = 0;
			Vector3i currentCellCoordinates = Coord(index);
			VectorD tangentVector = findTangent(pos, path);
			if (abs((pos - currentCellCoordinates).normalize().dot(tangentVector)) <= (grid.dx * grid.dx / 4))
			{
				//since time should be 0 relative to start point, this should give relative time to position
				real cell_eq_time = (((pos - currentCellCoordinates).dot(tangentVector))/grid.dx)*dt;


				cells.push_back(index);
				// Allocate uniform density value to each grid square here based on densityOpacityCurve
				densities[index] = densityOpacityCurve(cell_eq_time);

				// Allocate user-specified uniform temperature value, based upon how much time has passed, to each grid square
				temps[index] = getCurrentExplosionTemp();

				// Allocate direction velocity to each grid square according to u_d(G(g)) = V(t_i)t(g), where g is a grid square,
				// G(g) is the set of grid square in the region, and t(g) is the unit tangent vector derived from g on the flow control path
				// V(t_i) is the same for each grid square in the region

				// TODO: calculate the following variable
				// is supposed to be a real vector
				Vector3 unitVectorTangentToControlPath;
				velocities[i] = densityPropagationCurve(cell_eq_time, temps[i]) * tangentVector;
			}
		}
	}

	virtual void PreProcessing()
	{ //string dir_name){
		// 	Intialize Variables
		Array<std::string> file_names;

		// open through directory, find array of text files
		read_directory("/Users/student/explosions/NURBS/", file_names);

		//for every textfile, call read_from_file, add to array
		for (int i = 0; i < file_names.size(); i++)
		{
			controlPaths.push_back(read_from_file(file_names[i]));
		}


		//Assign times
		for(int i=0; i<controlPaths.size(); i++){
			real assign_time = 0.0;
			Array<real> times_to_push;
			for(int j=0; j<controlPaths[i].size(); j++){
				times_to_push.push_back(assign_time);
				assign_time += dt;
			}
			st_times.push_back(times_to_push);
		}

			//////TESTING
			std::ofstream outfile("test.txt");
			for (int i = 0; i < controlPaths.size(); i++)
			{
				outfile << "-------------------------" << std::endl;
				for (int j = 0; j < controlPaths[i].size(); j++)
				{
					outfile << std::to_string(controlPaths[i][j]) << " - " << std::to_string(st_times[i][j]) << std::endl;
				}
			}
			outfile.close();
			//////TESTING
		}

	///////////////////////////////////// HELPER FUNCTIONS ///////////////////////////////////////
protected:
	//get front of curve point

	/////////////////////////////////// Control Path Functions ///////////////////////////////////
	//find the length of a control path
	real find_path_length(const Array<int> &path)
	{
		real pathLength = 0.0;
		for (int i = 1; i < path.size(); i++)
		{
			pathLength += distanceNd(path[i], path[i - 1]);
		}
		return pathLength;
	}

  //shouldn't need this at all
	real getStartTimeFromIndex(const Array<int> &path, int index){
		return st_times[findInVector<int>(path, index)];
	}

	// give array and index of the array you are looking at
	// wait we need unqiue positions on the control path, so this functiom has to look for unique positions
	VectorD findTangent(int index, const Array<int> &path)
	{
		int node = path[index];
		int j = index;
		while (j>1 && j<path.size()-1 && path[j]==node){
			j--;
		}
		VectorD	pos1 = Pos(Coord(path[j]));
		j = index;
		while (j>1 && j<path.size()-1 && path[j]==node){
			j++;
		}
  	VectorD	pos2 = Pos(Coord(path[j]));
		VectorD result = (pos2-pos1).normalize();
		return result;
	}
	/////////////////////////////////// Reading Directories and Files ///////////////////////////////////

	void read_directory(const std::string &name, Array<std::string> &v)
	{
		DIR *dirp = opendir(name.c_str());
		struct dirent *dp;
		while ((dp = readdir(dirp)) != NULL)
		{
			size_t len = strlen(dp->d_name);
			if (len > 4 && strcmp(dp->d_name + len - 4, ".txt") == 0)
			{
				v.push_back(name + dp->d_name);
			}
		}
		closedir(dirp);
	}

	//read points from file and create control paths
	Array<int> read_from_file(std::string file_path)
	{ //what params?
		std::ifstream inFile;

		Array<int> grid_indices;

		//read from file given file path
		inFile.open(file_path);
		std::string data;

		//read from file, get string
		if (!inFile)
		{
			std::cout << "Unable to open file" << std::endl;
			exit(1); // terminate with error
		}

		for (std::string line; std::getline(inFile, line);)
		{
			Array<real> pointReals;
			VectorD point;
			std::istringstream in(line);
			for (std::string numb; std::getline(in, numb, ',');)
			{
				char *end;
				real push = std::strtof(numb.c_str(), &end);
				pointReals.push_back(push);
			}

			point = VectorD(pointReals.data());

			VectorDi coordinates = grid.Cell_Coord(point);

			int index = Idx(coordinates);
			grid_indices.push_back(index);
		}
		inFile.close();

		return grid_indices;
	}

	/////////////////////////////////// Physics Helper Functions ///////////////////////////////////

	Vector3i calculateColor(real temperature) {

		Vector3i returnColor;
		real adjustedTemp = temperature / 100;

		// Calculate red
		if (adjustedTemp <= 66) {
			returnColor[0] = 255;
		}
		else {
			returnColor[0] = adjustedTemp - 60;
			returnColor[0] = 329.698727446 * pow(returnColor[0], -0.1332047592);
			
			if (returnColor[0] < 0) {
				returnColor[0] = 0;
			}
			else if (returnColor[0] > 255) {
				returnColor[0] = 255;
			}
			
		}

		// Calculate green
		if (adjustedTemp <= 66) {
			returnColor[1] = adjustedTemp;
			returnColor[1] = 99.4708025861 * log(returnColor[1]) - 161.1195681661;

			if (returnColor[1] < 0) {
				returnColor[1] = 0;
			}
			else if (returnColor[1] > 255) {
				returnColor[1] = 255;
			}
		}
		else {
			returnColor[1] = adjustedTemp - 60;
			returnColor[1] = 288.1221695283 * pow(returnColor[1], -0.0755148492);
			
			if (returnColor[1] < 0) {
				returnColor[1] = 0;
			}
			else if (returnColor[1] > 255) {
				returnColor[1] = 255;
			}
		}

		// Calculate blue
		if (adjustedTemp >= 66) {
			returnColor[2] = 255;
		}
		else if (adjustedTemp <= 19){
			returnColor[2] = 0;
		}
		else {
			returnColor[2] = adjustedTemp - 10;
			returnColor[2] = 138.5177312231 * log(returnColor[2]) - 305.0447927307;
			if (returnColor[2] < 0) {
				returnColor[2] = 0;
			}
			else if (returnColor[2] > 255) {
				returnColor[2] = 255;
			}
		}

		returnColor /= 255;

		return returnColor;
	}

	// Returns pressure based on current time
	real pressureMagnitudeCurve(real t){
		real output = p_0;
		if (t <= t_d)
		{
			output += P_p * (1 - t / t_d) * exp((-b * t) / t_d);
		}
		else if (t <= t_d + t_n / 2)
		{
			output -= (2 * (P_m / t_n) * (t - t_d));
		}
		else if (t <= t_d + t_n)
		{
			output -= 2 * (P_m / t_n) * (t_d + t_n - t);
		}
		else
		{
			// do nothing, output should just be p_0
		}

		return output;
	}

	// Returns density based on current time
	real densityOpacityCurve(real t)
	{
		real output = (gamma + 1) * pressureMagnitudeCurve(t) + 2 * gamma * p_0;
		output /= ((gamma - 1) * pressureMagnitudeCurve(t) + 2 * gamma * p_0);
		output *= density_0;

		return output;
	}

	// Returns velocity magnitude v_p(t)
	real pressurePropagationCurve(real t, real temperature)
	{
		real c = getSpeedOfSoundInAir(temperature);
		return c * sqrt(1 + ((gamma + 1) * (pressureMagnitudeCurve(t) / (2 * gamma * p_0))));
	}

	// returns velocity magnitude v(t)
	real densityPropagationCurve(real t, real temperature)
	{

		real c = getSpeedOfSoundInAir(temperature);
		real output = 1 / sqrt(1 + (gamma + 1) * pressureMagnitudeCurve(t) / (2 * gamma * p_0));
		output *= c;
		output *= pressureMagnitudeCurve(t);
		output /= gamma;
		output /= p_0;

		return output;
	}

	// Returns speed of sound for a given temperature in m/s
	real getSpeedOfSoundInAir(real temperature)
	{
		return 331.5 + 0.6 * temperature;
	}

	// Returns uniform temperature value based on time elapsed
	real getCurrentExplosionTemp()
	{
		return explosionTemp - 100000 * time;
	}

	VectorD updateVorticity(const real dt, const real dx, VectorDi node){

		VectorD ip= velocities[Idx(node+VectorDi::Unit(0))];
		VectorD in= velocities[Idx(node-VectorDi::Unit(0))];
		VectorD jp= velocities[Idx(node+VectorDi::Unit(1))];
		VectorD jn= velocities[Idx(node-VectorDi::Unit(1))];
		VectorD kp= velocities[Idx(node+VectorDi::Unit(2))];
		VectorD kn= velocities[Idx(node-VectorDi::Unit(2))];

		Vector3 result = Vector3((jp[2]-jn[2])-(kp[1]-kn[1]),
														(kp[0]-kn[0])-(ip[2]-in[2]),
														(ip[1]-in[1])-(jp[0]-jn[0]));

		return((1.0/(2.0*dx)) * result);


	}

	/////////////////////////////////// Interpolation and Distance ///////////////////////////////////

	VectorD Interpolate(const Array<VectorD> &u, VectorD &pos)
	{
		//// define variables
		const real one_over_dx = 1.0 / grid.dx;
		////clamp pos
		for (int i = 0; i < d; i++)
		{
			if (pos[i] < grid.domain_min[i])
				pos[i] = grid.domain_min[i];
			else if (pos[i] > grid.domain_max[i])
				pos[i] = grid.domain_max[i];
		}

		VectorDi cell = ((pos - grid.domain_min) * one_over_dx).template cast<int>();
		VectorD frac = (pos * one_over_dx - cell.template cast<real>());
		return Interpolate_Helper(cell, frac, u);
	}

	////2D bi-linear interpolation
	Vector2 Interpolate_Helper(const Vector2i &cell, const Vector2 &frac, const Array<Vector2> &u)
	{
		return ((real)1 - frac[0]) * ((real)1 - frac[1]) * u[grid.Node_Index(cell)] + frac[0] * ((real)1 - frac[1]) * u[grid.Node_Index(Vector2i(cell[0] + 1, cell[1]))] + ((real)1 - frac[0]) * frac[1] * u[grid.Node_Index(Vector2i(cell[0], cell[1] + 1))] + frac[0] * frac[1] * u[grid.Node_Index(Vector2i(cell[0] + 1, cell[1] + 1))];
	}

	////3D tri-linear interpolation
	Vector3 Interpolate_Helper(const Vector3i &cell, const Vector3 &frac, const Array<Vector3> &u)
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
			lerp(u[Idx(Vector3i(cell[0], cell[1], cell[2] + 1))],
				 u[Idx(Vector3i(cell[0] + 1, cell[1], cell[2] + 1))],
				 frac[0]),
			lerp(u[Idx(Vector3i(cell[0], cell[1] + 1, cell[2] + 1))],
				 u[Idx(Vector3i(cell[0] + 1, cell[1] + 1, cell[2] + 1))],
				 frac[0]),
			frac[1],
			frac[2]);
	}

	template <class T>
	T lerp(T val1, T val2, real v) { return v * val2 + (1 - v) * val1; }

	inline real distanceNd(const int index1, const int index2)
	{
		VectorD pos1 = Pos(Coord(index1));
		VectorD pos2 = Pos(Coord(index2));
		real distance = 0.0;
		for (int i = 0; i < d; i++)
		{
			distance += (pos1[i] - pos2[i]) * (pos1[i] - pos2[i]);
		}
		return sqrt(distance);
	}

	/////////////////////////////////// Grid Helper Functions ///////////////////////////////////

	////return the node index given its coordinate
	int Idx(const Vector3i &node_coord) const
	{
		return grid.Node_Index(node_coord);
	}

	////return the coordinate given its index
	VectorDi Coord(const int node_index) const
	{
		return grid.Node_Coord(node_index);
	}

	////return the node position given its index
	VectorD Pos(const VectorDi& node) const
	{
		return grid.Node(node);
	}

	////check if a node is on the boundary of the grid
	////given its coordinate or index
	bool Bnd(const Vector3i &node_coord) const
	{
		for (int i = 0; i < d; i++)
		{
			if (node_coord[i] == 0 || node_coord[i] == grid.node_counts[i] - 1)
				return true;
		}
		return false;
	}
	bool Bnd(const int node_index) const
	{
		return Bnd(Coord(node_index));
	}

	////////////////////////////////// Other Helper Functions ///////////////////////////////////

	/*
	Generic function to find an element in vector and also its position.
	It returns a pair of bool & int i.e.
	bool : Represents if element is present in vector or not.
	int : Represents the index of element in vector if its found else -1
	*/
	template < typename T> int findInVector(const std::vector<T>  & vecOfElements, const T  & element)
	{
		int result;
		// Find given element in vector
		auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);
		if (it != vecOfElements.end()){ result = distance(vecOfElements.begin(), it);}
		else{ result = -1;}
		return result;
  }

};
