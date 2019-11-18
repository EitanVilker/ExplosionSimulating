//#####################################################################
// Opengl viewer
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include "OpenGLViewerDriver.h"

int main(int argc,char* argv[])
{
	std::string output_dir="output";
	for(int i=0;i<argc;i++){
		if(strcmp(argv[i],"-o")==0){
			output_dir=std::string(argv[++i]);}}
	std::cout<<"[OpenGLViewer arguments]: -o ="<<output_dir<<std::endl;

	std::shared_ptr<OpenGLViewer> viewer=nullptr;
	viewer.reset(new OpenGLViewerMassSpring());
	
	viewer->output_dir=output_dir;
	viewer->use_2d_display=false;

	viewer->Initialize();
	viewer->Set_Offscreen_Output_Dir(output_dir+"/_images");
	viewer->Run();
	viewer->Finish();

	return 0;
}
