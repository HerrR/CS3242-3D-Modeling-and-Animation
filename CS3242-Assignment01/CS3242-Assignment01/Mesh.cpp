#include "Mesh.h"
#include <typeinfo>
#include <map>
#include <algorithm>

// * Pekare (Peka på vart objektet ligger)
// & Referens (modifiera objektet direkt)

using namespace std;

Mesh::Mesh(const char* filename){	
	loadMF(filename);
}

void Mesh::loadMF(const char* filename){
	if(V.size()>0) V.clear();
	if(F.size()>0) F.clear();
    
	ifstream infile;
	infile.open(filename, ios::in);
    
	string strbuff;
	while(getline(infile,strbuff)){
		stringstream ss;
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
    cout << "V size: " << V.size() << endl;
    cout << "F size: " << F.size() << endl;
	infile.close();
}

void Mesh::writeMF(const char* filename){
	ofstream outfile;
	outfile.open(filename, ios::out);
	string strbuff;
	for(uint i=0;i<V.size();i++){
		outfile<<"v "<<V[i].x<<" "<<V[i].y<<" "<<V[i].z<<endl;
	}
	for(uint i=0;i<F.size();i++){
		outfile<<"f "<<F[i].a<<" "<<F[i].b<<" "<<F[i].c<<endl;
	}
	outfile.close();
}

void Mesh::simplifyMesh(const char* input, const char* output, int faceCnt){
	// you may assume inputs are always valid
    cout << "Loading mesh... \n";
	loadMF(input);
    
    cout << "Converting mesh... \n";
	convertMesh();
    
    cout << HEV.size() << " HE Vertices created" << endl;
    cout << HEF.size() << " HE Faces created " << endl;
    cout << HEE.size() << " HE Edges created" << endl;
    
	cout<< "Original face count: " << HEF.size()<<endl;
	// do mesh simplification here
    cout << "Simplify mesh here..."<<endl;
    
    cout<< "Output face count: " << HEF.size()<<endl;
    
    cout << "Reverting mesh..." << endl;
	revertMesh();
    cout << "Writing output file..." << endl;
    writeMF(output);
}


void Mesh::convertMesh()
{
    // Create HE Vertex for each normal vertex
    for (int i = 0; i < V.size(); i++) {
        HEVertex v;
        v.x = V[i].x;
        v.y = V[i].y;
        v.z = V[i].z;
        HEV.push_back(v);
    }
    
    // Key: Index of a Half Edge in the HEE vector.
    // Value: vertex index of the origin vertex in HEV.
    map<int, int> HasOriginAt;
    
    // Key: index of a HEE
    // Value: index in HEV of HEE.vertex
    map<int, int> IsPointingTo;
    
    // Key: Vertex index in HEV
    // Value: Vector of indexes in HEE of edges pointing to key.
    map<int, vector<int>> PointedAtBy;
    
    for (int i = 0; i < F.size(); i++) {
        // Create 3 HE Edges and one HE Face
        HEFace f;
        HEEdge e1, e2, e3;

        // Push edges to HEE Vector and face to HEF vector.
        HEF.push_back(f);
        HEE.push_back(e1);
        HEE.push_back(e2);
        HEE.push_back(e3);
        
        // Indexes of vertex 1,2 and 3 in HEV.
        int vIndex1 = F[i].a-1;
        int vIndex2 = F[i].b-1;
        int vIndex3 = F[i].c-1;
        
        // Indexes of e1, e2 and e3 in HEE.
        int eIndex1 = (int)HEE.size()-3;
        int eIndex2 = (int)HEE.size()-2;
        int eIndex3 = (int)HEE.size()-1;
        
        cout << "Vertex indexes: " << vIndex1 << ", " << vIndex2 << ", " << vIndex3 << endl;
        cout << "Edge indexes: " << eIndex1 << ", " << eIndex2 << ", " << eIndex3 << endl;
        
        // Add adjacent edge to face
        HEF[HEF.size()-1].edge = &HEE[eIndex1];
        
        // Assign one outgoing Half Edge for each
        HEV[vIndex1].edge = &HEE[eIndex1];
        HEV[vIndex2].edge = &HEE[eIndex2];
        HEV[vIndex3].edge = &HEE[eIndex3];
        
        // Assign vertex that each edge is pointing to
        HEE[eIndex1].vert = &HEV[vIndex2];
        HEE[eIndex2].vert = &HEV[vIndex3];
        HEE[eIndex3].vert = &HEV[vIndex1];
        
        // Assigns face for each HE
        HEE[eIndex1].face = &HEF[HEF.size()-1];
        HEE[eIndex2].face = &HEF[HEF.size()-1];
        HEE[eIndex3].face = &HEF[HEF.size()-1];
        
        // Assign null twins
        HEE[eIndex1].twin = NULL;
        HEE[eIndex2].twin = NULL;
        HEE[eIndex3].twin = NULL;
        
        // Assign next for each HE
        HEE[eIndex1].next = &HEE[eIndex2];
        HEE[eIndex2].next = &HEE[eIndex3];
        HEE[eIndex3].next = &HEE[eIndex1];
        
        // Update origin map for finding twins.
        HasOriginAt[eIndex1] = vIndex1;
        HasOriginAt[eIndex2] = vIndex2;
        HasOriginAt[eIndex3] = vIndex3;

        // Update which index the edges are pointing to
        IsPointingTo[eIndex1] = vIndex2;
        IsPointingTo[eIndex2] = vIndex3;
        IsPointingTo[eIndex3] = vIndex1;

        // Update PointedAtBy map.
        PointedAtBy[vIndex2].push_back(eIndex1);
        PointedAtBy[vIndex3].push_back(eIndex2);
        PointedAtBy[vIndex1].push_back(eIndex3);
        
    }
    
    for (int i = 0; i < HEE.size(); i++) {
        cout << &HEE[i] << endl;
        cout << &HEE[i].next << endl;
        cout << &HEE[i].next->next << endl;
        cout << &HEE[i].next->next->next << endl << endl;
    }
    
////     Loop through the half edge list, given a half edge, find its origin using the map
//    for (int i = 0; i < HEE.size(); i++) {
//        
//        for(int j = 0; j < PointedAtBy[HasOriginAt[i]].size(); j++){
//            if(HasOriginAt[PointedAtBy[HasOriginAt[i]][j]] == IsPointingTo[i]){
//                cout << "Twin found for HEE["<<i<<"] as HEE[" << PointedAtBy[HasOriginAt[i]][j] << "]" <<endl;
//                HEE[i].twin = &HEE[PointedAtBy[HasOriginAt[i]][j]];
////                HEE[PointedAtBy[HasOriginAt[i]][j]].twin = &HEE[i];
//            }
//        }
//    }
    
}

void Mesh::revertMesh()
{
//    V.clear();
//    F.clear();
    // Loop through each HE Face and create faces

}

vector<HEVertex> Mesh::neighborVertices(HEVertex v)
{
    // NeighborVertices of a vertex v, returns a vector of HEVertexes
    // Look at what is in the end of each outgoing half edge vertex
    vector<HEVertex> neighbours;

    return neighbours;
}

vector<HEFace> Mesh::neighborFaces(HEVertex v)
{
    vector<HEFace> adjacentFaces;
    // Get the neighboring faces of a face, look at the half edges adjacent half edges and find their close faces.
    return adjacentFaces;
}

vector<HEVertex> Mesh::adjacentVertices(HEFace f)
{
    // Get the vertices making up a face
    vector<HEVertex> adjacentVerts;
    return adjacentVerts;
}

int Mesh::Vcnt(){
    /* Half Edge Vertex */
	return (int)HEV.size();
}

int Mesh::Fcnt(){
    /* Half Edge Face */
	return (int)HEF.size();
}