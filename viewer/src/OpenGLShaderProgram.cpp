//#####################################################################
// Opengl shader program
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include <iostream>
#include "OpenGLBufferObjects.h"
#include "OpenGLShaderProgram.h"

//////////////////////////////////////////////////////////////////////////
////OpenGLShaders
namespace OpenGLShaders{
#define To_String(S) #S

//////////////////////////////////////////////////////////////////////////
////predefined uniform blocks and functions
const std::string version=To_String(
~version 330 core\n
);

const std::string material=To_String(
uniform vec4 mat_amb=vec4(1.f);
uniform vec4 mat_dif=vec4(1.f,1.f,1.f,1.f);
uniform vec4 mat_spec=vec4(1.f);
uniform vec4 mat_shinness=vec4(32.f,0.f,0.f,0.f);
);

const std::string phong_dl_func=To_String(
vec3 phong_dl(int i,vec3 norm)
{
	vec4 amb=lt[i].amb*mat_amb;
	vec3 light_dir=lt[i].dir.xyz;
    float dif_coef=max(dot(norm,-light_dir),0.);
    vec4 dif=dif_coef*lt[i].dif*mat_dif;
	vec4 color=amb+dif;return color.rgb;
}
);

const std::string phong_pl_func=To_String(
vec3 phong_pl(int i,vec3 pos,vec3 norm)
{
	vec4 amb=lt[i].amb*mat_amb;
	vec3 light_dir=lt[i].pos.xyz-pos;float dis=length(light_dir);light_dir=light_dir/dis;
    float dif_coef=max(dot(norm,light_dir),0.f);
    vec4 dif=dif_coef*mat_dif;
	vec3 view_dir=normalize(position.xyz-pos);
	vec3 half_dir=normalize(light_dir+view_dir);
	float spec_coef=pow(max(dot(norm,half_dir),0.f),mat_shinness[0]);
	vec4 spec=spec_coef*lt[i].spec*mat_spec;
	
	vec4 color=amb+dif+spec;
	float atten_coef=1.f/(lt[i].atten[0]+lt[i].atten[1]*dis+lt[i].atten[2]*dis*dis);
	color*=atten_coef;
	
	return color.rgb;
}
);

const std::string phong_sl_func=To_String(
vec3 phong_sl(int i,vec3 pos,vec3 norm)
{
	vec4 amb=lt[i].amb*mat_amb;
	vec3 light_dir=lt[i].pos.xyz-pos;float dis=length(light_dir);light_dir=light_dir/dis;
	float theta=dot(light_dir,-lt[i].dir.xyz);
	float spot_coef=clamp((theta-lt[i].r[1])/lt[i].r[2],0.f,1.f);

    float dif_coef=max(dot(norm,light_dir),0.f);
    vec4 dif=dif_coef*mat_dif;
	vec3 view_dir=normalize(position.xyz-pos);
	vec3 half_dir=normalize(light_dir+view_dir);
	float spec_coef=pow(max(dot(norm,half_dir),0.f),mat_shinness[0]);
	vec4 spec=spec_coef*lt[i].spec*mat_spec;

	vec4 color=amb+(dif+spec)*spot_coef;
	float atten_coef=1.f/(lt[i].atten[0]+lt[i].atten[1]*dis+lt[i].atten[2]*dis*dis);
	color*=atten_coef;

	return color.rgb;
}
);

const std::string phong_dl_fast_func=To_String(
vec3 phong_dl_fast(vec3 norm)
{
    float dif_coef=abs(dot(norm,vec3(1.f,1.f,1.f)));
    vec4 dif=dif_coef*vec4(.5f)*mat_dif+vec4(.1f);
	return dif.rgb;
}
);

//////////////////////////////////////////////////////////////////////////
////vtx shader
const std::string vpos_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
}														
);

const std::string vpos_model_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
uniform mat4 model=mat4(1.0f);
void main()												
{
	gl_Position=pvm*model*vec4(pos.xyz,1.f);
}														
);

const std::string vpos_model_vnormal_vfpos_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 normal;
out vec3 vtx_normal;
out vec3 vtx_frg_pos;
void main()												
{
	gl_Position=pvm*model*vec4(pos.xyz,1.f);
	vtx_normal=vec3(normal);
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
}														
);

const std::string vcolor_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
out vec4 vtx_color;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_color=v_color;
}														
);

const std::string vclip_vfpos_vtx_shader=To_String(
~include version;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
out vec3 vtx_frg_pos;
void main()
{
	gl_Position=vec4(pos.xyz,1.f);
	vtx_frg_pos=pos.xyz;
}
);

const std::string vnormal_vfpos_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 normal;
out vec3 vtx_normal;
out vec3 vtx_frg_pos;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_normal=vec3(normal);
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
}
);

const std::string psize_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
uniform float point_size=1.f;
void main()												
{
	gl_PointSize=point_size;
	gl_Position=pvm*vec4(pos.xyz,1.f);
}														
);

const std::string vcolor_ortho_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
out vec4 vtx_color;
void main()												
{
	gl_Position=ortho*model*vec4(pos.xy,1.f,1.f);
	vtx_color=v_color;
}														
);

const std::string vfpos_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
out vec3 vtx_frg_pos;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
}
);

//////////////////////////////////////////////////////////////////////////
////frg shader
const std::string none_frg_shader=To_String(
~include version;
void main(){}
);

const std::string ucolor_frg_shader=To_String(
~include version;
uniform vec4 color=vec4(1.f,1.f,0.f,1.f);
out vec4 frag_color;
void main()								
{										
    frag_color=color;	  
}										
);

const std::string gcolor_frg_shader=To_String(
~include version;
uniform vec4 color=vec4(1.f,1.f,0.f,1.f);
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()								
{ 
	vec3 c2=vec3(1.f,1.f,1.f);
	float m=abs(vtx_frg_pos.x);
	vec3 c=mix(c2,color.xyz,m*m);
	frag_color=vec4(c,1.f);
}										
);

const std::string vcolor_frg_shader=To_String(
~include version;
in vec4 vtx_color;
out vec4 frag_color;
void main()								
{										
    frag_color=vtx_color;
}										
);

const std::string vnormal_vfpos_lt_frg_shader=To_String(
~include version;
~include material;
~include camera;
~include lights;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
out vec4 frag_color;
~include phong_dl_func;
~include phong_pl_func;
~include phong_sl_func;
void main()
{
    vec3 normal=normalize(vtx_normal);
	vec3 color=mat_amb.rgb*amb.rgb;
	for(int i=0;i<lt_att[0];i++){
		vec3 c0=vec3(0.f);
		switch(lt[i].att[0]){
		case 0:{c0=phong_dl(i,normal);}break;
		case 1:{c0=phong_pl(i,vtx_frg_pos,normal);}break;
		case 2:{c0=phong_sl(i,vtx_frg_pos,normal);}break;}
		color+=c0;}
	frag_color=vec4(color,1.f);
}
);
}

const std::string vnormal_vfpos_dl_fast_frg_shader=To_String(
~include version;
~include material;
~include lights;
~include phong_dl_fast_func;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()
{
    vec3 norm=normalize(vtx_normal);
	vec3 color=phong_dl_fast(norm);
	frag_color=vec4(color,1.f);
}
);

using namespace OpenGLShaders;

//////////////////////////////////////////////////////////////////////////
////OpenGLShaderProgram

void OpenGLShaderProgram::Initialize(const std::string& vtx_shader_input,const std::string& frg_shader_input)
{vtx_shader=vtx_shader_input;frg_shader=frg_shader_input;compiled=false;
vtx_id=0;frg_id=0;geo_id=0;prg_id=0;
use_geo=false;geo_input_type=GL_POINTS;geo_output_type=GL_TRIANGLE_STRIP;}

void OpenGLShaderProgram::Initialize(const std::string& vtx_shader_input,const std::string& frg_shader_input,
    const std::string& _geo_shader_input,GLenum _geo_input_type,GLenum _geo_output_type,int _max_geo_vtx_output)
{Initialize(vtx_shader_input,frg_shader_input);
use_geo=true;geo_shader=_geo_shader_input;geo_input_type=_geo_input_type;
geo_output_type=_geo_output_type;max_geo_vtx_output=_max_geo_vtx_output;}
	
void OpenGLShaderProgram::Set_Uniform(const std::string& name,GLint value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform1i(location,value);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,GLfloat value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform1f(location,value);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,Vector2f value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform2f(location,value[0],value[1]);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,Vector3f value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform3f(location,value[0],value[1],value[2]);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,Vector4f value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform4f(location,value[0],value[1],value[2],value[3]);}

void OpenGLShaderProgram::Set_Uniform(const std::string& name,glm::vec2 value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform2f(location,value[0],value[1]);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,glm::vec3 value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform3f(location,value[0],value[1],value[2]);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,glm::vec4 value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform4f(location,value[0],value[1],value[2],value[3]);}

void OpenGLShaderProgram::Set_Uniform_Array(const std::string& name,GLsizei count,const GLint* value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform1iv(location,count,value);}
void OpenGLShaderProgram::Set_Uniform_Array(const std::string& name,GLsizei count,const GLfloat* value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform1fv(location,count,value);}

void OpenGLShaderProgram::Set_Uniform_Matrix4f(const std::string& name,const GLfloat* value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniformMatrix4fv(location,1,GL_FALSE,value);}
void OpenGLShaderProgram::Set_Uniform_Vec4f(const std::string& name,const GLfloat* value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform4f(location,value[0],value[1],value[2],value[3]);}

void OpenGLShaderProgram::Set_Uniform_Mat(const Material* mat)
{Set_Uniform("mat_amb",mat->mat_amb);Set_Uniform("mat_dif",mat->mat_dif);Set_Uniform("mat_spec",mat->mat_spec);Set_Uniform("mat_shinness",mat->mat_shinness);}

void OpenGLShaderProgram::Bind_Uniform_Block(const std::string& name,const GLuint binding_point)
{GLuint location=glGetUniformBlockIndex(prg_id,name.c_str());glUniformBlockBinding(prg_id,location,binding_point);}

void OpenGLShaderProgram::Bind_Texture2D(const std::string& name,GLuint tex_id,GLint tex_unit)
{GLuint location=glGetUniformLocation(prg_id,name.c_str());glActiveTexture(GL_TEXTURE0+tex_unit);
glBindTexture(GL_TEXTURE_2D,tex_id);glUniform1i(location,tex_unit);}

void OpenGLShaderProgram::Begin(){if(!compiled)Compile();glUseProgram(prg_id);}
void OpenGLShaderProgram::End(){glUseProgram(0);}

bool OpenGLShaderProgram::Compile()
{
	if(compiled)return true;

	vtx_id=glCreateShader(GL_VERTEX_SHADER);
	const char* vtx_shader_string=vtx_shader.c_str();
	GLint vtx_string_length=(GLint)vtx_shader.length()+1;
	glShaderSource(vtx_id,1,&vtx_shader_string,&vtx_string_length);
	glCompileShader(vtx_id);
	GLint vtx_compile_status;
	glGetShaderiv(vtx_id,GL_COMPILE_STATUS,&vtx_compile_status);
	if(vtx_compile_status!=GL_TRUE){
		char log[2048];int log_length;
		glGetShaderInfoLog(vtx_id,2048,(GLsizei*)&log_length,log);
		std::cerr<<"Error: [OpenGLShaderProgram] vertex shader compile error: "<<log<<std::endl;
		glDeleteShader(vtx_id);
		return false;}

	frg_id=glCreateShader(GL_FRAGMENT_SHADER);
	const char* frg_shader_string=frg_shader.c_str();
	GLint frg_string_length=(GLint)frg_shader.length()+1;
	glShaderSource(frg_id,1,&frg_shader_string,&frg_string_length);
	glCompileShader(frg_id);
	GLint frg_compile_status;
	glGetShaderiv(frg_id,GL_COMPILE_STATUS,&frg_compile_status);
	if(frg_compile_status!=GL_TRUE){
		char log[2048];int log_length;
		glGetShaderInfoLog(frg_id,2048,(GLsizei*)&log_length,log);
		std::cerr<<"Error: [OpenGLShaderProgram] fragment shader compile error: "<<log<<std::endl;
		glDeleteShader(frg_id);
		return false;}

	prg_id=glCreateProgram();
	glAttachShader(prg_id,vtx_id);
	glAttachShader(prg_id,frg_id);

	if(use_geo){
		geo_id=glCreateShader(GL_GEOMETRY_SHADER_EXT);
		const char* geo_shader_string=geo_shader.c_str();
		GLint geo_string_length=(GLint)geo_shader.length()+1;
		glShaderSource(geo_id,1,&geo_shader_string,&geo_string_length);
		glCompileShader(geo_id);
		GLint geo_compile_status;
		glGetShaderiv(geo_id,GL_COMPILE_STATUS,&geo_compile_status);
		if(geo_compile_status!=GL_TRUE){
			char log[2048];int log_length;
			glGetShaderInfoLog(geo_id,2048,(GLsizei*)&log_length,log);
			std::cerr<<"Error: [OpenGLShaderProgram] geometry shader compile error: "<<log<<std::endl;
			glDeleteShader(geo_id);
			return false;}

		glAttachShader(prg_id,geo_id);
		glProgramParameteriEXT(prg_id,GL_GEOMETRY_INPUT_TYPE_EXT,geo_input_type);
		glProgramParameteriEXT(prg_id,GL_GEOMETRY_OUTPUT_TYPE_EXT,geo_output_type);
		glProgramParameteriEXT(prg_id,GL_GEOMETRY_VERTICES_OUT_EXT,max_geo_vtx_output);}

	glLinkProgram(prg_id);
	GLint prg_link_status;
	glGetProgramiv(prg_id,GL_LINK_STATUS,&prg_link_status);
	if(prg_link_status!=GL_TRUE){
		char log[2048];int log_length;
		glGetShaderInfoLog(prg_id,2048,(GLsizei*)&log_length,log);
		std::cerr<<"Error: [OpenGLShaderProgram] program link error: "<<log<<std::endl;
		glDeleteProgram(prg_id);
		return false;}

	glDeleteShader(vtx_id);
	glDeleteShader(frg_id);
	if(use_geo)glDeleteShader(geo_id);
	compiled=true;
	return true;
}

//////////////////////////////////////////////////////////////////////////
////OpenGLShaderLibrary

OpenGLShaderLibrary* OpenGLShaderLibrary::Instance(){static OpenGLShaderLibrary instance;return &instance;}
std::shared_ptr<OpenGLShaderProgram> OpenGLShaderLibrary::Get(const std::string& name)
{
	auto search=shader_hashtable.find(name);
	if(search!=shader_hashtable.end())return search->second;
	else return std::shared_ptr<OpenGLShaderProgram>(nullptr);
}

OpenGLShaderLibrary::OpenGLShaderLibrary(){Initialize_Shaders();}
void OpenGLShaderLibrary::Initialize_Shaders()
{

	Initialize_Headers();
	Add_Shader(vpos_vtx_shader,ucolor_frg_shader,"vpos");
	Add_Shader(vpos_model_vtx_shader,ucolor_frg_shader,"vpos_model");

	Add_Shader(vcolor_vtx_shader,vcolor_frg_shader,"vcolor");
	Add_Shader(psize_vtx_shader,ucolor_frg_shader,"psize_ucolor");
	Add_Shader(vnormal_vfpos_vtx_shader,vnormal_vfpos_lt_frg_shader,"vnormal_lt");
	Add_Shader(vclip_vfpos_vtx_shader,gcolor_frg_shader,"gcolor_bk");
	Add_Shader(vpos_model_vnormal_vfpos_vtx_shader,vnormal_vfpos_dl_fast_frg_shader,"vpos_model_vnormal_dl_fast");
}

void OpenGLShaderLibrary::Initialize_Headers()
{
	shader_header_hashtable.insert(std::make_pair("version",version));
	shader_header_hashtable.insert(std::make_pair("material",material));
	shader_header_hashtable.insert(std::make_pair("phong_dl_func",phong_dl_func));
	shader_header_hashtable.insert(std::make_pair("phong_pl_func",phong_pl_func));
	shader_header_hashtable.insert(std::make_pair("phong_sl_func",phong_sl_func));
	shader_header_hashtable.insert(std::make_pair("phong_dl_fast_func",phong_dl_fast_func));

	OpenGLUbos::Bind_Shader_Ubo_Headers(shader_header_hashtable);
}

std::string OpenGLShaderLibrary::Parse(const std::string& shader) const
{
	std::string s=shader;
	std::replace(s.begin(),s.end(),'~','#');	////replace ~ with #, fix for linux compiling
	std::string name="#include";
	size_type p1=s.find(name);
	while(p1!=std::string::npos){
		size_type p2=s.find(' ',p1);
		size_type p3=s.find(';',p1);
		if(p2==std::string::npos||p3==std::string::npos)break;
		size_type n_var=p3-p2-1;
		std::string var=s.substr(p2+1,n_var);
		auto hash_pair=shader_header_hashtable.find(var);
		if(hash_pair==shader_header_hashtable.end())break;
		const std::string& replace=hash_pair->second;
		size_type n_replace=p3-p1+1;
		s.replace(p1,n_replace,replace);
		p1=s.find(name);
	}
	std::replace(s.begin(),s.end(),'~','#');	////replace ~ with #, fix for linux compiling
	return s;
}

void OpenGLShaderLibrary::Add_Shader(const std::string& vtx_shader,const std::string& frg_shader,const std::string& name)
{
	{std::shared_ptr<OpenGLShaderProgram> shader=std::make_shared<OpenGLShaderProgram>();
	shader->Initialize(Parse(vtx_shader),Parse(frg_shader));shader->name=name;
	shader_hashtable.insert(std::make_pair(shader->name,shader));}	
}

std::shared_ptr<OpenGLShaderProgram> OpenGLShaderLibrary::Get_Shader(const std::string& name)
{return OpenGLShaderLibrary::Instance()->Get(name);}
