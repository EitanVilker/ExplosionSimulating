//#####################################################################
// Opengl interactive object
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __OpenGLMarkerObjects_h__
#define __OpenGLMarkerObjects_h__
#include "glm.hpp"
#include "Mesh.h"
#include "OpenGLObject.h"
#ifdef USE_IMGUI
#include "OpenGLUiImgui.h"
#include "OpenGLViewer.h"
#include "OpenGLWindow.h"
#endif

class OpenGLBackground : public OpenGLObject
{public:typedef OpenGLObject Base;
	Box<2> box;	
	std::string tex_name;
	real depth;
	bool use_fbo_tex;

	OpenGLBackground();

	void Set_Box(const Vector2& min_corner,const Vector2& max_corner){box=Box<2>(min_corner,max_corner);}
	void Set_Texture(const std::string& _tex_name){use_vtx_tex=true;tex_name=_tex_name;}
	void Set_Depth(const real _depth){depth=_depth;}
	void Set_Fbo(){}
	virtual void Initialize();
	virtual void Display() const;
};

class OpenGLAxes : public OpenGLObject
{public:typedef OpenGLObject Base;
	real axis_length=(real).5;
	bool use_2d_display=false;

	OpenGLAxes(){name="axes";polygon_mode=PolygonMode::Wireframe;}

	virtual void Initialize();
	virtual void Display() const;
};

class OpenGLPoint : public OpenGLObject
{public:typedef OpenGLObject Base;using Base::color;
	Vector3 pos=Vector3::Zero();
	GLfloat point_size=16.f;

	OpenGLPoint(){name="point";color=OpenGLColor::Red();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
};

class OpenGLTriangle : public OpenGLObject
{public:typedef OpenGLObject Base;using Base::color;using Base::line_width;
	ArrayF<Vector3,3> vtx;

	OpenGLTriangle(){name="triangle";color=OpenGLColor::Red();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
};

class OpenGLPolygon : public OpenGLObject
{public:typedef OpenGLObject Base;using Base::color;using Base::line_width;
	Array<Vector3> vtx;

	OpenGLPolygon(){name="polygon";color=OpenGLColor::Blue();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
};

class OpenGLCircle : public OpenGLObject
{public: typedef OpenGLObject Base;
	Vector3 pos=Vector3::Zero();
	real radius=(real).1;
	glm::mat4 model=glm::mat4(1.f);

	int n=16;
	Array<Vector3> vtx;

	OpenGLCircle(){name="circle";color=OpenGLColor::Green();polygon_mode=PolygonMode::Fill;}

	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
	virtual void Update_Model_Matrix();
};

class OpenGLSolidCircle : public OpenGLObject
{public: typedef OpenGLObject Base;
	Vector3 pos=Vector3::Zero();
	real radius=(real).1;
	glm::mat4 model=glm::mat4(1.f);

	int n=16;
	Array<Vector3> vtx;

	OpenGLSolidCircle(){name="circle";color=OpenGLColor::Green();polygon_mode=PolygonMode::Fill;}

	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
	virtual void Update_Model_Matrix();
};

class OpenGLMarkerTriangleMesh : public OpenGLObject
{public:typedef OpenGLObject Base;using Base::color;
	TriangleMesh<3> mesh;
	glm::mat4 model=glm::mat4(1.f);

	OpenGLMarkerTriangleMesh(){name="interactive_triangle_mesh";color=OpenGLColor::Green();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
	virtual void Update_Model_Matrix(){}

protected:
	virtual void Update_Mesh_Data_To_Render();
};

class OpenGLSphere : public OpenGLMarkerTriangleMesh
{public:typedef OpenGLMarkerTriangleMesh Base;
	using Base::color;using Base::mesh;using Base::model;
	Vector3 pos=Vector3::Zero();
	real radius=(real).1;

	OpenGLSphere(){name="sphere";color=OpenGLColor::Red();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Model_Matrix();
};

//////////////////////////////////////////////////////////////////////////
////UI Object
#ifdef USE_IMGUI
class OpenGLUI : public OpenGLObject
{public:typedef OpenGLObject Base;using Base::name;using Base::Update_Polygon_Mode;
	Box<2> box;
	int mouse_x=-1,mouse_y=-1;
	int left=0,right=0,middle=0,scroll=0;
	bool is_mouse_hovering_window=false;
	OpenGLViewer* opengl_viewer=nullptr;

	OpenGLUI(){name="basic_ui";polygon_mode=PolygonMode::Fill;}

	virtual void Display() const
    {
		Update_Polygon_Mode();
		ImguiGlut::New_Frame_Imgui(mouse_x,mouse_y,0,right,middle,scroll);
		
		ImGui::SetNextWindowSize(ImVec2(256,256));
		ImGui::Begin("UI");
		ImGui::Text("Basic UI");
		ImGui::End();

		ImGui::Render();
    }

	virtual bool Mouse_Drag(int x,int y,int w,int h)
	{
		mouse_x=x;mouse_y=h-y;
		is_mouse_hovering_window=ImguiGlut::Update_Mouse_Position(mouse_x,mouse_y);
		return is_mouse_hovering_window;
	}

	virtual bool Mouse_Click(int _left,int _right,int _middle,int x,int y,int w,int h)
	{	
		left=_left;
		right=_right;middle=_middle;mouse_x=x;mouse_y=h-y;
		is_mouse_hovering_window=ImguiGlut::Update_Mouse_Position(mouse_x,mouse_y);
		return is_mouse_hovering_window;
	}

	virtual bool Mouse_Scroll(int wheel,int direction,int x,int y)
	{scroll=direction;return is_mouse_hovering_window;}
};


static const int text_size=64;
static char input_text[text_size]={0};

class OpenGLUICommand : public OpenGLUI
{public:typedef OpenGLUI Base;using Base::name;Base::Update_Polygon_Mode;using Base::visible;
	Array<std::string> commands;
	OpenGLUICommand():Base(){name="ui_command";}

	bool toggle_exit=false;

	virtual void Display() const
    {
		if(!visible)return;

		Update_Polygon_Mode();
		ImguiGlut::New_Frame_Imgui(mouse_x,mouse_y,0,right,middle,scroll);

		float w=460;float h=42;
		ImGui::SetNextWindowSize(ImVec2(w,h));
		ImGui::SetNextWindowPos(ImVec2(Win_Width()*.5f-w*.5f,Win_Height()*0.9f));
		ImGui::Begin("Command Window",0,ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoResize);
		ImGui::SetWindowFontScale(2.1f);
		ImGui::PushItemWidth(w-60);
		ImGui::PushStyleColor(ImGuiCol_Text,ImVec4(1.f,1.f,0.f,1.f));
		ImGui::PushStyleColor(ImGuiCol_PopupBg,ImVec4(.1f,1.f,.1f,1.f));
		ImGui::SetKeyboardFocusHere();
		bool is_done=ImGui::InputText("Cmd",input_text,text_size,ImGuiInputTextFlags_EnterReturnsTrue);
		ImGui::PopStyleColor(2);
		ImGui::PopItemWidth();
		ImGui::End();
		ImGui::Render();

		if(is_done){for(int i=0;i<text_size;i++){input_text[i]=0;}}
		if(toggle_exit){
			auto non_const_ptr=const_cast<OpenGLUICommand*>(this);
			non_const_ptr->visible=false;
			non_const_ptr->toggle_exit=false;}
    }

	virtual bool Keyboard(unsigned char key,int x,int y,bool down)
	{
		if(key==(unsigned char)13&&down&&!visible){visible=true;return false;}
		if(!visible){return false;}
		////Enter and get command
		if(key==(unsigned char)13&&down&&visible){
			std::string cmd=Get_Command();
			std::cout<<"Parse cmd: "<<cmd<<std::endl;
			commands.push_back(cmd);
			if(opengl_viewer!=nullptr)opengl_viewer->Toggle_Command(cmd);
			if(cmd=="exit"){toggle_exit=true;}}
		////Process input key
		ImguiGlut::Update_Keyboard(key,x,y,down);
		return true;
	}

	virtual bool Keyboard_Special(int key,int x,int y,bool down)
	{
		if(!visible){return false;}
		ImguiGlut::Update_Keyboard_Special(key,x,y,down);
		return true;
	}

	std::string Get_Command() const
	{
		std::string cmd;
		int i=0;for(;i<text_size;i++){if(input_text[i]==0)break;}
		cmd.append(&input_text[0],i);
		return cmd;
	}
};
#endif	//USE_IMGUI

#endif
