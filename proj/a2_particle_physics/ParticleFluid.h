//#####################################################################
// Particle Fluid (SPH)
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#ifndef __ParticleFluid_h__
#define __ParticleFluid_h__
#include "Common.h"
#include "Particles.h"
#include "ImplicitGeometry.h"

//////////////////////////////////////////////////////////////////////////
////Kernel function
template<int d> class Kernel
{using VectorD=Vector<real,d>;
public:
	////precomputed coefs;
	real h;
	real coef_Wspiky;
	real coef_dWspiky;
	real coef_Wvis;
	real coef_d2Wvis;
	real pi=3.1415927;

	void Precompute_Coefs(real _h)
	{
		h=_h;
		coef_Wspiky=15.0/(pi*pow(h,6));
		coef_dWspiky=-45.0/(pi*pow(h,6));
		coef_Wvis=2*pi*pow(h,3);
		coef_d2Wvis=45.0/(pi*pow(h,6));
	}

	////Kernel Spiky
	real Wspiky(const VectorD& xji)
	{
		real r=xji.norm();
		if(r>=0&&r<=h){return 15.0/(pi*pow(h,6))*pow(h-r,3);}
		else{return 0;}
	}
	VectorD gradientWspiky(const VectorD& v){
		real r=v.norm();
		if(r<= h&&r>0){return -45.0/(pi*pow(h,6))*pow(h-r,2)*v/r;}
		else{return VectorD::Zero();}
	}

	////Kernel viscosity
	real Wvis(const VectorD& xji){
		real r=xji.norm();
		if(r>=0&&r<=h){return 15.0/(2*pi*pow(h,3))*((-pow(r,3)/(2*pow(h,3))+r*r/(h*h)+h/(2*r)-1));}
		else{return 0;}
	}
	real laplacianWvis(const VectorD& v){
		real r=v.norm();
		if(r<=h&&r>0){return 45.0/(pi*pow(h,6))*(h-r);}
		else{return 0;}
	}
};

//////////////////////////////////////////////////////////////////////////
////Spatial hashing
template<int d> class SpatialHashing
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;
public:
	real dx=1.;	////grid cell size
	Hashtable<VectorDi,Array<int> > voxels;

	void Update_Voxels(const Array<VectorD>& points)
	{Clear_Voxels();for(int i=0;i<(int)points.size();i++)Add_Point(i,points[i]);}

	void Clear_Voxels(){voxels.clear();}

	bool Add_Point(const int point_idx,const VectorD& point_pos)
	{
		VectorDi cell=Cell_Coord(point_pos);
		auto iter=voxels.find(cell);
		if(iter==voxels.end())iter=voxels.insert(std::make_pair(cell,Array<int>())).first;
		Array<int>& bucket=iter->second;
		bucket.push_back(point_idx);
		return true;
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): find all the neighboring particles within the "kernel_radius" around "pos" and record their indices in "nbs", the position of the particles are given in "points"
	////You need to traverse all the 3^d neighboring cells in the background grid around the cell occupied by "pos", and then check the distance between each particle in each neighboring cell and the given "pos"
	////Use the helper function Cell_Coord to get the cell coordinates for a given "pos"
	////Use the helper function Nb_R to get the cell coordinates of the ith neighboring cell around the cell "coord"
	bool Find_Nbs(const VectorD& pos,const Array<VectorD>& points,const real kernel_radius,/*returned result*/Array<int>& nbs) const
	{
		auto cell = Cell_Coord(pos);
		for (int i = 0; i < pow(3, d); i++)
		{
			auto it = voxels.find(Nb_R(cell, i));
			if (it == voxels.end()) continue;
			for (auto j : it->second) 
			{
				if ((points[j] - pos).norm() <= kernel_radius) nbs.emplace_back(j);
			}
		}
		return nbs.size()>0;
	}

protected:	////Helper functions
	VectorDi Cell_Coord(const VectorD& pos) const
	{VectorD coord_with_frac=(pos)/dx;return coord_with_frac.template cast<int>();}
	Vector2i Nb_R(const Vector2i& coord,const int index) const
	{assert(index>=0&&index<9);int i=index/3;int j=index%3;return coord+Vector2i(-1+i,-1+j);}
	Vector3i Nb_R(const Vector3i& coord,const int index) const
	{assert(index>=0&&index<27);int i=index/9;int m=index%9;int j=m/3;int k=m%3;return coord+Vector3i(-1+i,-1+j,-1+k);}
};

//////////////////////////////////////////////////////////////////////////
////Particle fluid simulator
template<int d> class ParticleFluid
{using VectorD=Vector<real,d>;
public:
	Particles<d> particles;
	Array<Array<int> > neighbors;
	SpatialHashing<d> spatial_hashing;
	Kernel<d> kernel;

	real kernel_radius=(real).8;			////kernel radius
	real pressure_density_coef=(real)1e1;	////pressure-density-relation coefficient, used in Update_Pressure()
	real density_0=(real)10.;				////rest density, used in Update_Pressure()
	real viscosity_coef=(real)1e1;			////viscosity coefficient, used in Update_Viscosity_Force()
	real kd=(real)1e2;						////stiffness for environmental collision response
	VectorD g=VectorD::Unit(1)*(real)-1.;	////gravity
	
	////Environment objects
	Array<ImplicitGeometry<d>* > env_objects;

	virtual void Initialize()
	{
		kernel.Precompute_Coefs(kernel_radius);
	}

	virtual void Update_Neighbors()
	{
		spatial_hashing.Clear_Voxels();
		spatial_hashing.Update_Voxels(particles.XRef());

		neighbors.resize(particles.Size());
		for(int i=0;i<particles.Size();i++){
			Array<int> nbs;
			spatial_hashing.Find_Nbs(particles.X(i),particles.XRef(),kernel_radius,nbs);
			neighbors[i]=nbs;}
	}

	virtual void Advance(const real dt)
	{
		for(int i=0;i<particles.Size();i++){
			particles.F(i)=VectorD::Zero();}

		Update_Neighbors();
		Update_Density();
		Update_Pressure();
		Update_Pressure_Force();
		Update_Viscosity_Force();
		Update_Body_Force();
		Update_Boundary_Collision_Force();

		for(int i=0;i<particles.Size();i++){
			particles.V(i)+=particles.F(i)/particles.D(i)*dt;
			particles.X(i)+=particles.V(i)*dt;}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): update the density (particles.D(i)) of each particle based on the kernel function (Wspiky)
	void Update_Density()
	{
		for (int particle = 0; particle < particles.Size(); particle++)
		{
			particles.D(particle) = 0;
			for (auto neighbor : neighbors[particle])
			{
				auto xji = particles.X(neighbor) - particles.X(particle);
				particles.D(particle) += particles.M(neighbor) * kernel.Wspiky(xji);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): update the pressure (particles.P(i)) of each particle based on its current density (particles.D(i)) and the rest density (density_0)
	void Update_Pressure()
	{
		for (int particle = 0; particle < particles.Size(); particle++)
		{
			particles.P(particle) = pressure_density_coef * (particles.D(particle) - density_0);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): compute the pressure force for each particle based on its current pressure (particles.P(i)) and the kernel function gradient (gradientWspiky), and then add the force to particles.F(i)
	void Update_Pressure_Force()
	{
		for (int particle = 0; particle < particles.Size(); particle++)
		{
			for (auto neighbor : neighbors[particle])
			{
				auto xji = particles.X(particle) - particles.X(neighbor);
				particles.F(particle) -= (particles.P(neighbor) + particles.P(particle)) / 2 *
					particles.M(neighbor) / particles.D(neighbor) * kernel.gradientWspiky(xji);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): compute the viscosity force for each particle based on its current velocity difference (particles.V(j)-particles.V(i)) and the kernel function Laplacian (laplacianWvis), and then add the force to particles.F(i)
	void Update_Viscosity_Force()
	{
		for (int particle = 0; particle < particles.Size(); particle++)
		{
			for (auto neighbor : neighbors[particle])
			{
				auto xji = particles.X(particle) - particles.X(neighbor);
				particles.F(particle) += viscosity_coef * 
					(particles.V(neighbor) - particles.V(particle)) *
					particles.M(neighbor) / particles.D(neighbor) * kernel.laplacianWvis(xji);
			}
		}
	}

	void Update_Body_Force()
	{
		for(int i=0;i<particles.Size();i++){
			particles.F(i)+=particles.D(i)*g;}	
	}

	void Update_Boundary_Collision_Force()
	{
		for(int i=0;i<particles.Size();i++){
			for(int j=0;j<env_objects.size();j++){
				real phi=env_objects[j]->Phi(particles.X(i));
				if(phi<particles.R(i)){
					VectorD normal=env_objects[j]->Normal(particles.X(i));
					particles.F(i)+=normal*kd*(particles.R(i)-phi)*particles.D(i);}}}
	}
};

#endif
