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
    
    cout << "\nConverting mesh... \n";
	convertMesh();
    
    cout << HEV.size() << " HE Vertices created" << endl;
    cout << HEF.size() << " HE Faces created " << endl;
    cout << HEE.size() << " HE Edges created \n" << endl;
    
    cout << "Original face count: " << HEF.size() << endl;
	// do mesh simplification here
    
//    collapseEdge(*HEE[6]);
    int nextDeletion;
    srand((int)time(0));
    while(HEF.size() > faceCnt){
        nextDeletion = rand() % HEE.size();
//        cout << "Deleting edge " << nextDeletion << "...\n";
        collapseEdge(*HEE[nextDeletion]);
    }

    
//    cout << endl << "Simplification goes here..."<< endl << endl;
    
    cout<< "\nOutput face count: " << HEF.size()<<endl;
    
    cout << "\nReverting mesh..." << endl;
	revertMesh();
    cout << "\nWriting output file..." << endl;
    writeMF(output);
}

void Mesh::collapseEdge(HEEdge e){
    // TODO
    // 4 Half Edges will be linked through 2 new twin connections CHECK
    // Update vertex pointers for all HEs pointing at vertex to be deleted CHECK
    
    // Delete adjacent faces CHECK
    // Delete one of the vertices (the one being pointed to by HE selected)
    // Delete all HE's in face being deleted
    
    vector<HEEdge> deleteThese = Mesh::neighborEdges(e.face);
    vector<HEEdge> deleteTheseToo = Mesh::neighborEdges(e.twin->face);
    deleteThese.insert(deleteThese.end(), deleteTheseToo.begin(), deleteTheseToo.end());
    
    // Update twin pairs
//    cout << "Setting twin of " <<
//        e.next->twin->id << " pointing at " <<
//        e.next->twin->vert->id << " to " <<
//        e.next->next->twin->id << " pointing at " <<
//        e.next->next->twin->vert->id << endl;
    e.next->twin->twin = e.next->next->twin;
    
//    cout << "Setting twin of " <<
//        e.next->next->twin->id << " pointing at "<<
//        e.next->next->twin->vert->id << " to " <<
//        e.next->twin->id << " pointing at " <<
//        e.next->twin->vert->id << endl;
    e.next->next->twin->twin = e.next->twin;
    
//    cout << "Setting twin of " << e.twin->next->twin->twin->id << " to " <<  e.twin->next->next->twin->id << endl;
//    cout << "Setting twin of " <<
//    e.twin->next->twin->id << " pointing at "<<
//    e.twin->next->twin->vert->id << " to " <<
//    e.twin->next->next->twin->id << " pointing at " <<
//    e.twin->next->next->twin->vert->id << endl;
    e.twin->next->twin->twin = e.twin->next->next->twin;
    
//    cout << "Setting twin of " <<
//    e.twin->next->next->twin->id << " pointing at "<<
//    e.twin->next->next->twin->vert->id << " to " <<
//    e.twin->next->twin->id << " pointing at " <<
//    e.twin->next->twin->vert->id << endl;
//    cout << "Setting twin of " << e.twin->next->next->twin->twin->id << " to " <<  e.twin->next->twin->id << endl << endl;
    e.twin->next->next->twin->twin = e.twin->next->twin;
    
//    cout << "\nNew vertex to be pointed to " << e.twin->vert->id << endl;
//    cout << "Vertex to be deleted: " << e.vert->id << endl;
    
    HEEdge* nextEdge = e.next->twin;
    do{
//        cout << "Repointing edge "<< nextEdge->id <<" currently pointing at " << nextEdge->vert->id << " to " << e.twin->vert->id << endl;
        nextEdge->vert = e.twin->vert;
        nextEdge = nextEdge->next->twin;
    } while(nextEdge->vert->id != e.twin->vert->id);
    
//    cout << endl;

    for(int i = 0; i<HEV.size(); i++){
        if(HEV[i]->id == e.vert->id){
//            cout << "Deleting vertex with id " << HEV[i]->id << endl;
            HEV.erase(HEV.begin()+i);
            break;
        }
    }

//    cout << "Faces to be deleted: " << e.face->id << " and " << e.twin->face->id << endl;
    int deletedFaces = 0;
    for (int i = 0; i<HEF.size(); i++) {
        if( (HEF[i]->id == e.face->id) || (HEF[i]->id == e.twin->face->id) ){
//            cout << "Deleting face with id " << HEF[i]->id << endl;
            HEF.erase(HEF.begin()+i);
            deletedFaces ++;
            if(deletedFaces == 2){
                break;
            }
            i --;
        }
    }

    int deletedEdges = 0;
    bool stepBack = false;
    for(int i = 0; i<HEE.size(); i++){
        if(stepBack){
            i--;
            stepBack = false;
        }

        for (int j = 0; j < deleteThese.size(); j++) {
            if(HEE[i]->id == deleteThese[j].id){
//                cout << "Deleting edge  " << HEE[i]->id << endl;
                HEE.erase(HEE.begin()+i);
                deleteThese.erase(deleteThese.begin()+j);
                
                if(deletedEdges == deleteThese.size()){
                    break;
                } else {
                    deletedEdges ++;
                }
                if(i == 0){
                    stepBack = true;
                } else {
                    i--;
                }
            }
        }
    }
}


void Mesh::convertMesh()
{
    // Create HE Vertex for each normal vertex
    for (int i = 0; i < V.size(); i++) {
        HEVertex* v = new HEVertex();
        v->x = V[i].x;
        v->y = V[i].y;
        v->z = V[i].z;
        v->id = i;
        HEV.push_back(v);
    }
    
    // Key: ID of Half Edge in HEE
    // Value: ID of HEV at origin of key
    map<int, int> HasOriginAt;
    
    // Key: ID of Half Edge in HEE
    // Value: ID of vertex being pointed to
    map<int, int> IsPointingTo;
    
    // Key: Vertex index in HEV
    // Value: Vector of indexes in HEE of edges pointing to key.
    map<int, vector<int>> PointedAtBy;
    
    for (int i = 0; i < F.size(); i++) {
        // Create 3 HE Edges and one HE Face
//        HEFace f;
        HEFace *f = new HEFace();
        HEEdge *e1 = new HEEdge();
        HEEdge *e2 = new HEEdge();
        HEEdge *e3 = new HEEdge();

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
//        
//        cout << "Iteration " << i << endl;
//        cout << "Vertex indexes: " << vIndex1 << ", " << vIndex2 << ", " << vIndex3 << endl;
//        cout << "Edge indexes: " << eIndex1 << ", " << eIndex2 << ", " << eIndex3 << endl << endl;

        
        // Add adjacent edge to face
        HEF[i]->edge = HEE[eIndex1];
        HEF[i]->id = i;
        
        // Assign IDs to edges
        HEE[eIndex1]->id = eIndex1;
        HEE[eIndex2]->id = eIndex2;
        HEE[eIndex3]->id = eIndex3;
        
        // Assign one outgoing Half Edge for each
        HEV[vIndex1]->edge = HEE[eIndex1];
        HEV[vIndex2]->edge = HEE[eIndex2];
        HEV[vIndex3]->edge = HEE[eIndex3];
        
        // Assign vertex that each edge is pointing to
        HEE[eIndex1]->vert = HEV[vIndex2];
        HEE[eIndex2]->vert = HEV[vIndex3];
        HEE[eIndex3]->vert = HEV[vIndex1];
        
        // Assigns face for each HE
        HEE[eIndex1]->face = HEF[i];
        HEE[eIndex2]->face = HEF[i];
        HEE[eIndex3]->face = HEF[i];

        // Assign next for each HE
        HEE[eIndex1]->next = HEE[eIndex2];
        HEE[eIndex2]->next = HEE[eIndex3];
        HEE[eIndex3]->next = HEE[eIndex1];
        
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
    
//     Loop through the half edge list, given a half edge, find its origin using the map
    for (int i = 0; i < HEE.size(); i++) {
        
        for(int j = 0; j < PointedAtBy[HasOriginAt[i]].size(); j++){
            
            if(HasOriginAt[PointedAtBy[HasOriginAt[i]][j]] == IsPointingTo[i]){
//                cout << "Twin found for "<< HEE[i]->id+1 <<" as " << HEE[PointedAtBy[HasOriginAt[i]][j]]->id+1 << "" <<endl;
                HEE[i]->twin = &(*HEE[PointedAtBy[HasOriginAt[i]][j]]);
                break;
            }
        }
    }
}

void Mesh::revertMesh()
{
    V.clear();
    F.clear();
    map<int, int> locatedAt;
    for (int i = 0; i<HEV.size(); i++) {
        Vertex v;
        v.x = HEV[i]->x;
        v.y = HEV[i]->y;
        v.z = HEV[i]->z;
        V.push_back(v);
        locatedAt[HEV[i]->id] = i+1;
    }
    
    for(int i = 0; i<HEF.size(); i++){
        Face f;
        vector<HEVertex> adjacentVerts = adjacentVertices(*HEF[i]);
        f.a = locatedAt[adjacentVerts[0].id];
        f.b = locatedAt[adjacentVerts[1].id];
        f.c = locatedAt[adjacentVerts[2].id];
//        cout << "Created face with vertices " << f.a << ", " << f.b << ", " << f.c << endl;
        F.push_back(f);
    }
}
vector<HEEdge> Mesh::neighborEdges(HEFace* f){
    vector<HEEdge> neighborEdges;
    
    HEEdge* firstEdge = f->edge;
    HEEdge* nextEdge = firstEdge->next;
    neighborEdges.push_back(*nextEdge);
    do{
        nextEdge = nextEdge->next;
        neighborEdges.push_back(*nextEdge);
    } while(nextEdge->id != firstEdge->id);
    
    return neighborEdges;
}

vector<HEVertex> Mesh::neighborVertices(HEVertex v)
{
    // NeighborVertices of a vertex v, returns a vector of HEVertexes
    // Look at what is in the end of each outgoing half edge vertex
    vector<HEVertex> neighbourVerts;
    
    HEVertex *firstVertex = v.edge->vert;
    HEEdge *nextEdge = v.edge;
    
    do{
        nextEdge = nextEdge->twin->next;
        neighbourVerts.push_back(*nextEdge->vert);
    } while(nextEdge->vert != firstVertex);
    
    return neighbourVerts;
}



vector<HEFace> Mesh::neighborFaces(HEVertex v)
{
    vector<HEFace> adjacentFaces;
    
    HEEdge *firstEdge = v.edge;
    HEEdge *nextEdge = v.edge->next;
    adjacentFaces.push_back(*nextEdge->twin->face);
    do{
        adjacentFaces.push_back(*nextEdge->twin->face);
        nextEdge = nextEdge->next;
    } while(nextEdge != firstEdge);
    
    return adjacentFaces;
}

vector<HEVertex> Mesh::adjacentVertices(HEFace f)
{
    // Get the vertices making up a face
    vector<HEVertex> adjacentVerts;
    
    HEEdge *firstEdge = f.edge;
    HEEdge *nextEdge = f.edge->next;
    
    adjacentVerts.push_back(*firstEdge->vert);
    
    do{
        adjacentVerts.push_back(*nextEdge->vert);
        nextEdge = nextEdge->next;
    } while(nextEdge != firstEdge);
    
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