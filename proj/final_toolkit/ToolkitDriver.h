#ifndef __GridFluidDriver_h__
#define __GridFluidDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "Particles.h"
#include "OpenGLMesh.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLParticles.h"

template<int d> class ToolkitDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	real dt=.02;
	GridFluid<d> fluid;

	OpenGLTriangleMesh* opengl_tri_mesh=nullptr;						////bunny
	OpenGLSegmentMesh* opengl_vectors=nullptr;							////vector field
	Array<OpenGLSphere*> opengl_spheres;								////spheres
	OpenGLPolygon* opengl_polygon=nullptr;

public:
	virtual void Initialize()
	{
		fluid.Initialize();
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		////initialize bunny
		opengl_tri_mesh=Add_Interactive_Object<OpenGLTriangleMesh>();
		auto tri_mesh=&opengl_tri_mesh->mesh;
		std::string file_name="bunny.txt";
		File::Read_Text_From_File(file_name,*tri_mesh);
		Set_Polygon_Mode(opengl_tri_mesh,PolygonMode::Fill);
		Set_Shading_Mode(opengl_tri_mesh,ShadingMode::Lighting);
		opengl_tri_mesh->Set_Data_Refreshed();
		opengl_tri_mesh->Initialize();

		////initialize vector field (write all vectors as segments)
		opengl_vectors=Add_Interactive_Object<OpenGLSegmentMesh>();
		opengl_vectors->mesh.elements.resize(fluid.node_num);
		opengl_vectors->mesh.vertices->resize(fluid.node_num*2);
		for(int i=0;i<fluid.node_num;i++){
			opengl_vectors->mesh.elements[i]=Vector2i(i*2,i*2+1);
			(*opengl_vectors->mesh.vertices)[i*2]=V3(fluid.grid.Node(fluid.grid.Node_Coord(i)));
			VectorD pos2=fluid.grid.Node(fluid.grid.Node_Coord(i))+fluid.u[i];
			(*opengl_vectors->mesh.vertices)[i*2+1]=V3(pos2);}
		opengl_vectors->Set_Data_Refreshed();
		opengl_vectors->Initialize();

		////initialize sphere
		for(int i=0;i<3;i++){
			OpenGLSphere* opengl_sphere=Add_Interactive_Object<OpenGLSphere>();
			opengl_sphere->pos=Vector3::Unit(0)*((real).5*(real)i+(real)2);
			opengl_sphere->radius=(real).02;
			Set_Color(opengl_sphere,OpenGLColor(.0+.5*(real)i,1.,.0,1.));
			opengl_sphere->Set_Data_Refreshed();
			opengl_sphere->Initialize();
			opengl_spheres.push_back(opengl_sphere);}

		////initialize polygon
		opengl_polygon=Add_Interactive_Object<OpenGLPolygon>();
		opengl_polygon->vtx.push_back(Vector3::Zero());
		opengl_polygon->vtx.push_back(Vector3::Unit(0)*(real)2);
		opengl_polygon->vtx.push_back(Vector3::Unit(0)*(real)2+Vector3::Unit(1)*(real)1);
		opengl_polygon->vtx.push_back(Vector3::Unit(1)*(real)1);
		Set_Color(opengl_polygon,OpenGLColor(.0,1.,1.,1.));
		Set_Line_Width(opengl_polygon,4.f);
		opengl_polygon->Set_Data_Refreshed();
		opengl_polygon->Initialize();

		////set lights
		auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
		OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
		OpenGLUbos::Update_Lights_Ubo();
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		////sync data for vectors
		for(int i=0;i<fluid.node_num;i++){
			VectorD pos2=fluid.grid.Node(fluid.grid.Node_Coord(i))+fluid.u[i];
			(*opengl_vectors->mesh.vertices)[i*2+1]=V3(pos2);}
		opengl_vectors->Set_Data_Refreshed();

		////update and sync data for spheres
		for(int i=0;i<opengl_spheres.size();i++){
			opengl_spheres[i]->pos+=Vector3::Unit(1)*.1;
			opengl_spheres[i]->Set_Data_Refreshed();}
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		fluid.Advance(dt);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}

	////Keyboard interaction
	virtual void Initialize_Common_Callback_Keys()
	{
		OpenGLViewer::Initialize_Common_Callback_Keys();
		Bind_Callback_Key('a',&Keyboard_Event_A_Func,"press A");
	}

	virtual void Keyboard_Event_A()
	{
		std::cout<<"A pressed"<<std::endl;
	}
	Define_Function_Object(ToolkitDriver,Keyboard_Event_A);

protected:
	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2){return Vector3(v2[0],v2[1],.0);}
	Vector3 V3(const Vector3& v3){return v3;}
};
#endif