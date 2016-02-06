#include "Mesh.h"
#include <iostream>

int main(int argc, char* argv[]){
    std::cout << "Hellooo!" << "\n";
    Mesh mesh;
//    mesh.loadMF("data/Bunny6000.mesh");
    mesh.simplifyMesh("/Users/rickardbergeling/GitHub/CS3242-3D-Modeling-and-Animation/CS3242-Assignment01/CS3242-Assignment01/data/Bunny6000.mesh", "data/output.mesh", 600);
	/* if(argc!=4) std::cout<<"Usage: ./exe input output faceCnt\n";
	else{
		Mesh mesh;
		mesh.simplifyMesh(argv[1],argv[2],atoi(argv[3]));
	}
	return 0;*/
    return 0;
}