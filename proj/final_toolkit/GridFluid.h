#ifndef __GridFluid_h__
#define __GridFluid_h__
#include "Common.h"
#include "Grid.h"

//////////////////////////////////////////////////////////////////////////
////Particle fluid simulator
template<int d> class GridFluid
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;
public:
	Grid<d> grid;
	Array<VectorD> u;
	
	SparseMatrixT A;
	VectorX p;
	VectorX div_u;

	int node_num;
	real one_over_dx;

	virtual void Initialize()
	{
		int n=32;
		VectorDi cell_counts=VectorDi::Ones()*n;
		real dx=(real)1./(real)n;
		one_over_dx=(real)1/dx;
		real one_over_dx_sq=one_over_dx*one_over_dx;
		VectorD domain_min=VectorD::Unit(0)*5;
		grid.Initialize(cell_counts,dx,domain_min);
		node_num=grid.node_counts.prod();
		u.resize(node_num,VectorD::Unit(0)*(real).01);

		A.resize(node_num,node_num);

		p.resize(node_num);div_u.fill((real)0);
		div_u.resize(node_num);div_u.fill((real)0);
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
		Advection(dt);
		Projection();
	}

	////Helper functions
protected:
	VectorD Interpolate(const Array<VectorD>& u,VectorD& pos)
	{
		////clamp pos
		for(int i=0;i<d;i++){
			if(pos[i]<grid.domain_min[i])pos[i]=grid.domain_min[i];
			else if(pos[i]>grid.domain_max[i])pos[i]=grid.domain_max[i];}

		VectorDi cell=((pos-grid.domain_min)*one_over_dx).template cast<int>();
		VectorD frac=(pos*one_over_dx-cell.template cast<real>());
		return Interpolate_Helper(cell,frac,u);
	}

	////2D bi-linear interpolation
	Vector2 Interpolate_Helper(const Vector2i& cell,const Vector2& frac,const Array<Vector2>& u)
	{
		return ((real)1-frac[0])*((real)1-frac[1])*u[grid.Node_Index(cell)]
			+frac[0]*((real)1-frac[1])*u[grid.Node_Index(Vector2i(cell[0]+1,cell[1]))]
			+((real)1-frac[0])*frac[1]*u[grid.Node_Index(Vector2i(cell[0],cell[1]+1))]
			+frac[0]*frac[1]*u[grid.Node_Index(Vector2i(cell[0]+1,cell[1]+1))];
	}

	////3D tri-linear interpolation
	Vector3 Interpolate_Helper(const Vector3i& cell,const Vector3& frac,const Array<Vector3>& u)
	{
		return Vector3::Zero();	
		////your implementation here
	}
};

#endif
