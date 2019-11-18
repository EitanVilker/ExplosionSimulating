//#####################################################################
// Opengl viewer driver
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __OpenGLViewerDriver_h__
#define __OpenGLViewerDriver_h__
#include "OpenGLViewer.h"
#include "OpenGLMesh.h"
#include "OpenGLParticles.h"

class OpenGLViewerMassSpring : public OpenGLViewer
{using Base=OpenGLViewer;
public:
	virtual void Initialize_Data()
	{

		////draw particles
		auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");
		Set_Visibility(p,'x',true);
		Set_Color(p,OpenGLColor::Red());
		Set_Point_Size(p,6.f);

		////draw fixed particles
		auto psi_D=Add_Object<OpenGLParticles<Particles<3> > >("psi_D");
		Set_Visibility(psi_D,'x',true);
		Set_Color(psi_D,OpenGLColor::Yellow());
		Set_Point_Size(psi_D,12.f);
	
		////draw triangle surface
		auto tri=Add_Object<OpenGLTriangleMesh>("triangle_mesh");
		Set_Polygon_Mode(tri,PolygonMode::Fill);
		Set_Visibility(tri,'t',true);
		if(tri){
			tri->recomp_vtx_normal=true;
			tri->use_mat=true;
			tri->mat.mat_dif=glm::vec4(.1f,0.3f,1.f,1.f);
			tri->mat.mat_spec=glm::vec4(.1f,0.3f,1.f,1.f);}
		
		////draw segment mesh
		auto seg=Add_Object<OpenGLColoredSegmentMesh>("segment_mesh");
		Set_Visibility(seg,'s',true);
		Set_Line_Width(seg,4.f);
		Set_Color(seg,OpenGLColor(1.f,1.f,1.f,1.f));

		////set lights
		auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
		OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
		OpenGLUbos::Update_Lights_Ubo();
	}
};

#endif
