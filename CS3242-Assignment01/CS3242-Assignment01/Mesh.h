#ifndef _MESH_H_
#define _MESH_H_

#include "Eigen/Dense"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

typedef unsigned int uint;

//data structure of indexed face set
typedef struct {
	//3d coordinates
	float x;
	float y;
	float z;
} Vertex;

typedef struct{
	//three vertex ids
	uint a,b,c;
} Face;

struct HEVertex; struct HEFace; struct HEEdge;

typedef struct HEVertex{
    float x;
    float y;
    float z;
    int id;
    HEEdge* edge; // Outgoing halfedge
} HEVertex;

typedef struct HEEdge {
    HEVertex* vert;   // Vertex it points to
    HEEdge* twin;   // Oppositely oriented adjacent halfedge
    HEFace* face;   // Adjacent face
    HEEdge* next;   // Next halfedge around the face
    int id;
} HEEdge;

typedef struct HEFace {
    HEEdge* edge;  // one of the halfedges bordering the face
    int id;
} HEFace;

class Mesh{
private:
	vector<Vertex> V;
	vector<Face> F;

	vector<HEVertex*> HEV;
	vector<HEEdge*> HEE;
	vector<HEFace*> HEF;
public:
	Mesh() {};
	Mesh(const char*);
	//load a Mesh from .mesh file
	void loadMF(const char*);
	//write a Mesh to .mesh file (no header)
	void writeMF(const char*);
	//simplify a mesh
	void simplifyMesh(const char* input, const char* output, int faceCnt);
	//turn indexed face set to halfedge
	void convertMesh();
    // Collapse edge
    void collapseEdge(HEEdge *e);
	//turn halfedge to indexed face set
	void revertMesh();
	//helper methods
    Eigen::Matrix4f getQ(HEVertex* v);
    vector<HEEdge*> neighborEdges(HEFace* f);
	vector<HEVertex*> neighborVertices(HEVertex* v);
	vector<HEFace*> neighborFaces(HEVertex* v);
	vector<HEVertex*> adjacentVertices(HEFace* f);
    vector<HEEdge*> associatedEdges(HEVertex* v);
    
	//return vertex count
	int Vcnt();
	//return face count
	int Fcnt();
};
#endif