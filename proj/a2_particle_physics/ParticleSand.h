//#####################################################################
// Particle Sand (DEM)
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################

#ifndef __ParticleSand_h__
#define __ParticleSand_h__
#include "Common.h"
#include "Particles.h"
#include "ImplicitGeometry.h"

template<int d> class ParticleSand
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using MatrixD=Matrix<real,d>;
public:
	Particles<d> particles;
	real ks=(real)2e2;		////stiffness for the collision force
	real kd=(real).5e1;		////damping for the collision force
	VectorD g=VectorD::Unit(1)*(real)-1.;	////gravity

	Array<ImplicitGeometry<d>* > env_objects;	////list of implicit geometries describing the environment, by default it has one element, a circle with its normals pointing inward (Bowl)
	Array<Vector2i> particle_particle_collision_pairs;
	Array<Vector2i> particle_environment_collision_pairs;
	
	virtual void Advance(const real dt)
	{
		////Clear forces on particles
		for(int i=0;i<particles.Size();i++){
			particles.F(i)=VectorD::Zero();}

		////Accumulate body forces
		for(int i=0;i<particles.Size();i++){
			particles.F(i)+=particles.M(i)*g;}

		Particle_Environment_Collision_Detection();
		Particle_Environment_Collision_Response();
		Particle_Particle_Collision_Detection();
		Particle_Particle_Collision_Response();

		for(int i=0;i<particles.Size();i++){
			particles.V(i)+=particles.F(i)/particles.M(i)*dt;
			particles.X(i)+=particles.V(i)*dt;}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): detect collision between particles and env_objects (implicit surface) and record the detected collisions in particle_environment_collision_pairs
	////env_objects is a list of implicit geometries, by default there is only one geometry (the bowl) in the list
	////Each element in particle_environment_collision_pairs is a Vector2i, with the first index for the particle and the second index for the env_objects
	virtual void Particle_Environment_Collision_Detection()
	{
		particle_environment_collision_pairs.clear();
		for (int implicit_id = 0; implicit_id < env_objects.size(); implicit_id++)
		{
			for (int particle_id = 0; particle_id < particles.Size(); particle_id++)
			{
				if (env_objects[implicit_id]->Phi(particles.X(particle_id)) < particles.R(particle_id))
				{
					particle_environment_collision_pairs.emplace_back(Vector2i(particle_id, implicit_id));
				}
			}
		}
	}
		
	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): compute the penalty-based collision force for the particles that are colliding with the env_objects
	////The collision response force consists of a spring force and a damping force
	virtual void Particle_Environment_Collision_Response()
	{
		for(int pair_idx=0;pair_idx<particle_environment_collision_pairs.size();pair_idx++){
			int i=particle_environment_collision_pairs[pair_idx][0];	////particle index
			int j=particle_environment_collision_pairs[pair_idx][1];	////env_objects index

			VectorD n_ij = -env_objects[j]->Normal(particles.X(i));

			//Add Spring force
			particles.F(i)+= ks * n_ij * (env_objects[j]->Phi(particles.X(i)) - particles.R(i));
			//Add Damping force
			particles.F(i) += kd * n_ij.dot(-particles.V(i)) * n_ij;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): find all the pairs of particles that are colliding each other and record the detected pairs in particle_particle_collision_pairs
	////Each element in particle_particle_collision_pairs is a Vector2i specifying the indices of the two colliding particles
	virtual void Particle_Particle_Collision_Detection()
	{
		particle_particle_collision_pairs.clear();

		for (int i = 0; i < particles.Size() - 1; i++) 
		{
			for (int j = i + 1; j < particles.Size(); j++)
			{
				if ((particles.X(i) - particles.X(j)).norm() < (particles.R(i) + particles.R(j)))
				{
					particle_particle_collision_pairs.emplace_back(Vector2i(i, j));
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): compute penalty-based collision forces for pairs of colliding particles in particle_particle_collision_pairs and add the forces to particle.F 
	////The collision response force for each pair consists of a spring force and a damping force
	virtual void Particle_Particle_Collision_Response()
	{
		for(int pair_idx=0;pair_idx<particle_particle_collision_pairs.size();pair_idx++){
			int i=particle_particle_collision_pairs[pair_idx][0];	////the first particle index in the pair
			int j=particle_particle_collision_pairs[pair_idx][1];	////the second particle index in the pair

			VectorD n_ij = (particles.X(j) - particles.X(i)).normalized();
			
			VectorD spring_force = ks * ((particles.X(j) - particles.X(i)).norm() - (particles.R(i) + particles.R(j))) * n_ij;
			VectorD dampin_force = kd * (particles.V(j) - particles.V(i)).dot(n_ij) * n_ij;

			particles.F(i) += spring_force + dampin_force;
			particles.F(j) -= spring_force + dampin_force;
		}

	}
};

#endif
