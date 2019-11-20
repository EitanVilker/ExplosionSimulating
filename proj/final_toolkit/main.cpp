//#####################################################################
// Main
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#include <iostream>
#include "GridFluid.h"
#include "ToolkitDriver.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
	int driver=1;
	string directory_path;

	// if there is an arg, it is the directory path
	if(argc>0){
		directory_path = argv[0];
	}
	// should have an else saying yo need a directory path

	switch(driver){
	case 1:{

		//pass directory path to driver
		ToolkitDriver<3> driver;
		driver.Initialize();
		driver.Run();
	}break;
	}

}

#endif
