#pragma once

#include "Common.h"
#include "Grid.h"
#include "DCPQuery.h"
#include "Particles.h"

//////////////////////////////////////////////////////////////////////////
////Particle fluid simulator
template<int d> class ExplosionSim
{
	using VectorD = Vector<real, d>; using VectorDi = Vector<int, d>;
public:
	////////// Data Structure Variables
	Grid<d> grid;
	Array<VectorD> velocities;							//Velocity    on grid cells
	Array<real>    densities;							//Density     on grid cells
	Array<real>	   temps;								//Temperature on grid cells
	Array<real>    pressures;							//Pressure    on grid cells
	Array<real>    div_vels;							//Divergence  on grid cells
	Array<real>	   vorticities;							//Vortices    on grid cells

	////////// Other Constants and Variables

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
	real heatFlux = 286478.9;			//Heat flux (J/m^2)

	//time
	real time;							//Current time passed (seconds)
	real dt = 0.01;						//Timestep (seconds)
	const int d = 3;					//Number of dimensions

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
		densities.resize(node_num, d_amb);
		temps.resize(node_num, temp_amb);
		pressures.resize(node_num, pres_atm);
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
			//TODO: Get nurbs curve value
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
		
		for (int i = 0; i < cells.Size(); i++) {

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


////////// Helper functions
protected:
	//get front of curve point
	VectorD get_distance_traveled_chord(const real t, const real dt, const real temperature){
		// call velocity from 0 to t by dt
	}

	inline real scale_by_distance(const real value, const real dist_trav, const real total_length){

	}

	real find_path_length(Array<VectorD>){

	}

	void read_from_file(){
		
	}

	// Returns pressure based on current time
	real pressureMagnitudeCurve(real t) {

		real output = p_0;
		if (t <= t_d) {
			output += P_p * (1 - t / t_d) * exp((-bt) / t_d);
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
		return c * sqrt(1 + ((gamma + 1) * (pressureMagnitudeCurve(t) / (2 * gamma * p_0)));
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

	VectorD Interpolate(const Array<VectorD>& u, VectorD& pos)
	{
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
				lerp(u[Idx(cell),
					u[Idx(Vector3i(cell[0] + 1, cell[1], cell[2]))],
					frac[0]),
				lerp(u[Idx(Vector3i(cell[0], cell[1] + 1, cell[2]))],
					u[Idx(Vector3i(cell[0] + 1, cell[1] + 1, cell[2]))],
					frac[0]),
				frac[1]),
				lerp(u[Idx(Vector3i(cell[0], cell[1], cell[2] + 1),
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

	/*
	/////
	/////
	Grid helper functions: REVIEW NECESSARY
	/////
	/////
	/////
	/////
	*/
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