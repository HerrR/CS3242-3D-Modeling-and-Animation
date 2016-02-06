#include "Mesh.h"

// * Pekare (Peka på vart objektet ligger)
// & Referens (modifiera objektet direkt)

Mesh::Mesh(const char* filename){	
	loadMF(filename);
}

void Mesh::loadMF(const char* filename){
	if(V.size()>0) V.clear();
	if(F.size()>0) F.clear();
    
	std::ifstream infile;
	infile.open(filename, std::ios::in);
    
	std::string strbuff;
	while(std::getline(infile,strbuff)){
		std::stringstream ss;
		ss<<strbuff;
		char type;
		ss>>type;
		if(type=='v'){
			Vertex v;
			ss>>v.x>>v.y>>v.z;
			V.push_back(v);
		}
		else if(type=='f'){
			Face f;
			ss>>f.a>>f.b>>f.c;
			F.push_back(f);
		}
	}
    std::cout << "V size: " << V.size() << std::endl;
    std::cout << "F size: " << F.size() << std::endl;
	infile.close();
}

void Mesh::writeMF(const char* filename){
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	std::string strbuff;
	for(uint i=0;i<V.size();i++){
		outfile<<"v "<<V[i].x<<" "<<V[i].y<<" "<<V[i].z<<std::endl;
	}
	for(uint i=0;i<F.size();i++){
		outfile<<"f "<<F[i].a<<" "<<F[i].b<<" "<<F[i].c<<std::endl;
	}
	outfile.close();
}

void processMesh(){
    std::cout << "Process mesh initialized!" << std::endl;
//    std::cout << typeof(Mesh);
}

void Mesh::simplifyMesh(const char* input, const char* output, int faceCnt){
	// you may assume inputs are always valid
	loadMF(input);
	processMesh();
	std::cout<<"Original face count: "<<HEF.size()<<std::endl;
	// do mesh simplification here

	//reverseProcessMesh();
	//writeMF();
}



void Mesh::convertMesh()
{
}

void Mesh::revertMesh()
{
}

//std::vector<HEVertex> Mesh::neighborVertices(HEVertex v)
//{
//}
//
//std::vector<HEFace> Mesh::neighborFaces(HEVertex v)
//{
//}
//
//std::vector<HEVertex> Mesh::adjacentVertices(HEFace f)
//{
//
//}

int Mesh::Vcnt(){
    /* Half Edge Vertex */
	return HEV.size();
}

int Mesh::Fcnt(){
    /* Half Edge Face */
	return HEF.size();
}