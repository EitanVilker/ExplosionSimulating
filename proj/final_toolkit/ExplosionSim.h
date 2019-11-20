#pragma once

#include "Common.h"
#include "Grid.h"

//////////////////////////////////////////////////////////////////////////
////Particle fluid simulator
template<int d> class ExplosionSim
{
	using VectorD = Vector<real, d>; using VectorDi = Vector<int, d>;
public:
	Grid<d> grid;
	Array<VectorD> velocities;							//Velocity    on grid cells
	Array<real>    densities;							//Density     on grid cells
	Array<real>	   temps;								//Temperature on grid cells
	Array<real>    pressures;							//Pressure    on grid cells
	Array<real>    div_vels;							//Divergence  on grid cells
	Array<real>	   vorticities;							//Vortices    on grid cells

	// Constants
	int node_num;
	int n_per_dim = 32;
	const real d_amb = 1.4;				//Ratio of specific heats at ambient pressure (no units)
	const real temp_amb = 25;			//Ambient temperature- Celcius
	const real p_0 = 1;					//Ambient pressure- atm
	real P_p = 1;						//Peak overpressure- need to determine appropriate value (NTD)
	real P_m = -1;						//Minimum negative pressure (NTD)
	real t_d = 3;						//Time to reach p_0 (NTD)
	real t_n = 3;						//Time to reach p_0 from t_d (NTD)
	real b;								//Decreasing coefficient (NTD)
	real density_0;						//Ambient density (NTD)- paper mentions that this is specifically set by the user

	real time;							//current time
	bool explosion_done;				

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

		//TODO: Initing particles
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

	virtual void SweepRegion(const VectorD& pos, Array<int>& cells) 
	{
		//TODO: get the sweep region around pos and put into cells
	}

	////Helper functions
protected:

	// This curve takes as input various predetermined times and coefficients (seems somewhat arbitrary)
	// Everything but t should be a predetermined value that is global to ExplosionSim
	// Returns pressure based on current time
	virtual real pressureMagnitudeCurve(real t) {

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
	real densityOpacityCurve(real t, real gamma) {

		real output = (gamma + 1) * pressureMagnitudeCurve(t) + 2 * gamma * p_0;
		output /= ((gamma - 1) * pressureMagnitudeCurve(t) + 2 * gamma * p_0);
		output *= density_0;

		return output;
	}

	// Returns velocity magnitude v_p(t)
	real pressurePropagationCurve(real t, real gamma, real temperature) {

		real c = getSpeedOfSoundInAir(temperature);
		return c * sqrt(1 + ((gamma + 1) * (pressureMagnitudeCurve(t) / (2 * gamma * p_0)));
	}

	// Returns speed of sound for a given temperature in m/s
	real getSpeedOfSoundInAir(real temperature) {

		return 331.5 + 0.6 * temperature;
	}

	// returns velocity magnitude v(t)
	real densityPropagationCurve(real t, real c, real gamma) {

		real output = 1 / sqrt(1 + (gamma + 1) * pressureMagnitudeCurve(t) / (2 * gamma * p_0));
		output *= c;
		output *= pressureMagnitudeCurve(t);
		output /= gamma;
		output /= p_0;

		return output;
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
		return Vector3::Zero();
		////your implementation here
	}
};