#include "Mesh.h"
#include "Eigen/Dense"
#include <typeinfo>
#include <map>
#include <algorithm>

using namespace std;

#define DEBUG 1

typedef struct {
    float Q;
    HEEdge *e;
} EQTuple;

bool tupleSort(EQTuple* i, EQTuple* j) { return (i->Q<j->Q); }
bool edgeSort(HEEdge* i, HEEdge* j) { return (i->id < j->id); }


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

Eigen::Matrix4f Mesh::getQforEdge(HEEdge* e){
    // Get the neighbouring faces of a vertex v
    vector<HEFace*> neighborFaces = Mesh::neighborFaces(e->vert);
    
    // Place to store planeErrorQuadrics
    vector<Eigen::Matrix4f> planeErrorQuadrics;
    
    // Calculate the plane error quadric for each face
    for(int j = 0; j<neighborFaces.size(); j++){
        // Get the three corner vertexes of the plane
        HEVertex firstVertex = *neighborFaces[j]->edge->vert; // First point
        HEVertex secondVertex = *neighborFaces[j]->edge->next->vert; // Second point
        HEVertex thirdVertex = *neighborFaces[j]->edge->next->next->vert; // Third point
        
        // Calculate two vectors from the three points
        Eigen::Vector3f v1(secondVertex.x - firstVertex.x, secondVertex.y - firstVertex.y, secondVertex.z - firstVertex.z);
        Eigen::Vector3f v2(thirdVertex.x - firstVertex.x, thirdVertex.y - firstVertex.y, thirdVertex.z - firstVertex.z);
        
        // Calculate plane normal
        Eigen::Vector3f faceNormal = v1.cross(v2);
        float d = -(faceNormal[0] * e->vert->x + faceNormal[1] * e->vert->y + faceNormal[2] * e->vert->z);
        
        // Calculate plane error quadric
        Eigen::Vector4f facePlaneEquation(faceNormal[0], faceNormal[1], faceNormal[2], d);
        Eigen::Matrix4f Kp = facePlaneEquation*facePlaneEquation.transpose();
        planeErrorQuadrics.push_back(Kp);
    }
    
    // Sum the plane error quadric for each face
    Eigen::Matrix4f Q = Eigen::MatrixXf::Zero(4, 4);
    for(int j = 0; j < planeErrorQuadrics.size(); j++){
        Q += planeErrorQuadrics[j];
    }
    
    return Q;
}


Eigen::Matrix4f Mesh::getQ(HEVertex* v){
    
    // Get the neighbouring faces of a vertex v
    vector<HEFace*> neighborFaces = Mesh::neighborFaces(v);
    
    // Place to store planeErrorQuadrics
    vector<Eigen::Matrix4f> planeErrorQuadrics;
    
    // Calculate the plane error quadric for each face
    for(int j = 0; j<neighborFaces.size(); j++){
        // Get the three corner vertexes of the plane
        HEVertex firstVertex = *neighborFaces[j]->edge->vert; // First point
        HEVertex secondVertex = *neighborFaces[j]->edge->next->vert; // Second point
        HEVertex thirdVertex = *neighborFaces[j]->edge->next->next->vert; // Third point
        
        // Calculate two vectors from the three points
        Eigen::Vector3f v1(secondVertex.x - firstVertex.x, secondVertex.y - firstVertex.y, secondVertex.z - firstVertex.z);
        Eigen::Vector3f v2(thirdVertex.x - firstVertex.x, thirdVertex.y - firstVertex.y, thirdVertex.z - firstVertex.z);
        
        // Calculate plane normal
        Eigen::Vector3f faceNormal = v1.cross(v2);
        float d = -(faceNormal[0] * v->x + faceNormal[1] * v->y + faceNormal[2] * v->z);
        
        // Calculate plane error quadric
        Eigen::Vector4f facePlaneEquation(faceNormal[0], faceNormal[1], faceNormal[2], d);
        Eigen::Matrix4f Kp = facePlaneEquation*facePlaneEquation.transpose();
        planeErrorQuadrics.push_back(Kp);
    }
    
    // Sum the plane error quadric for each face
    Eigen::Matrix4f Q = Eigen::MatrixXf::Zero(4, 4);
    for(int j = 0; j < planeErrorQuadrics.size(); j++){
        Q += planeErrorQuadrics[j];
    }
    
    return Q;
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
    
    // Map with all the vertex IDs and their error quadric
    map<int, Eigen::Matrix4f> errorQuadricsMap;
    
    // Get the Q matrices for each vertex
    if(DEBUG) cout << "Getting Q matrixes for all vertexes... " << endl;
    for(int i = 0; i < HEV.size() ; i++){
        Eigen::Matrix4f Q = getQ(HEV[i]);
        errorQuadricsMap[HEV[i]->id] = Q;
    }
    
    // QList is a vector of tuples containing a reference to a half edge and an associated Q value which is the cost of deleting that edge
    vector<EQTuple*> QList;
    
    // Get Q values for all the edges
    if(DEBUG) cout << "Getting Q values for all edges..." << endl;
    for(int i = 0; i<HEE.size(); i++){
        Eigen::Vector4f v(HEE[i]->twin->vert->x, HEE[i]->twin->vert->y, HEE[i]->twin->vert->z, 1);
        EQTuple *e = new EQTuple{
            v.transpose() * errorQuadricsMap[HEE[i]->vert->id] * v, // Q
            HEE[i] // e
        };
        QList.push_back(e);
    }
    
    cout << "Wanting to reduce faces from " << HEF.size() << " to " << faceCnt << endl;
    cout << "\nReducing mesh..." << endl;
    
    while((HEF.size()) > faceCnt && (HEF.size()>2)){
        // Sort the QList so that the edge generating the least error gets to first place
        sort(QList.begin(), QList.end(), tupleSort);
        while(true){
            if(QList.size() == 1)break;
            if(edgeDeletionMap[QList[0]->e] || vertexDeletionMap[QList[0]->e->vert]){
                QList.erase(QList.begin());
                cout << "Disaster avoided!" << endl;
            } else {
                break;
            }
        }
        
        if(DEBUG) cout << "\nFirst in QList: edge " << QList[0]->e->id << " with twin " << QList[0]->e->twin->id << endl;
        if(DEBUG){
            // Look for EQ Tuples in QList which contains an edge that has been deleted.
            for(int i = 0; i<QList.size(); i++){
                if(edgeDeletionMap[QList[i]->e]){
                    cout << "----- Edge in QList not valid "<< QList[i]->e->id << endl;
                }
            }
        }
        
        // Find the edges that are going to be deleted (those within the adjacent faces of e and e->twin
        vector<HEEdge*> f1Edges = neighborEdges(QList[0]->e->face);
        vector<HEEdge*> f2Edges = neighborEdges(QList[0]->e->twin->face);
        vector<HEEdge*> edgesToBeDeleted;
        
        edgesToBeDeleted.insert(edgesToBeDeleted.end(), f1Edges.begin(), f1Edges.end());
        edgesToBeDeleted.insert(edgesToBeDeleted.end(), f2Edges.begin(), f2Edges.end());
        
        if(DEBUG){
            cout << "Edges to be deleted: " << endl;
            for(int i = 0; i < edgesToBeDeleted.size(); i++){
                cout << edgesToBeDeleted[i]->id << " ";
            }
            cout << endl;
        }
        
        // Find the edges that are going to be updated (whose Q value will be affected by deletion
        vector<HEEdge*> edgesToBeUpdated;
        //        vector<HEEdge*> e1Edges = associatedEdges(QList[0]->e->vert);
        //        vector<HEEdge*> e2Edges = associatedEdges(QList[0]->e->twin->vert);
        vector<HEEdge*> e1Edges = edgesPointingToSameVert(QList[0]->e);
        vector<HEEdge*> e2Edges = edgesPointingToSameVert(QList[0]->e->twin);
        
        edgesToBeUpdated.insert(edgesToBeUpdated.end(), e1Edges.begin(), e1Edges.end());
        edgesToBeUpdated.insert(edgesToBeUpdated.end(), e2Edges.begin(), e2Edges.end());
        
        // Don't need to update the edges that are going to be deleted
        for(int i = 0; i<edgesToBeUpdated.size(); i++){
            for(int j = 0; j < edgesToBeDeleted.size(); j++){
                if(edgesToBeUpdated[i] == edgesToBeDeleted[j]){
                    edgesToBeUpdated.erase(edgesToBeUpdated.begin()+i);
                    i--;
                }
            }
        }
        
        if(DEBUG) cout << "\nCollapsing edge " << QList[0]->e->id << " with twin " << QList[0]->e->twin->id << " \n";
        collapseEdge(QList[0]->e);
        
        // Reposition vertex
        QList[0]->e->twin->vert->x = (QList[0]->e->vert->x + QList[0]->e->twin->vert->x)/2;
        QList[0]->e->twin->vert->y = (QList[0]->e->vert->y + QList[0]->e->twin->vert->y)/2;
        QList[0]->e->twin->vert->z = (QList[0]->e->vert->z + QList[0]->e->twin->vert->z)/2;
        
        // Delete EQ Tuples for edges that have been deleted from the mesh
        vector<EQTuple*>::iterator it;
        int numDeletions = 0;
        for(it=QList.begin(); it != QList.end();){
            if(find(edgesToBeDeleted.begin(), edgesToBeDeleted.end(), (*it)->e) != edgesToBeDeleted.end()){
                if(DEBUG) cout << "Deleting EQ tuple with ref to edge " << (*it)->e->id << " \n";
                QList.erase(it);
                numDeletions ++;
                --it;
            }
            ++it;
        }
        
        // Update Q value for edges affected by the change
        int numUpdates = 0;
        for(int i = 0; i<QList.size(); i++){
            for(int j = 0; j < edgesToBeUpdated.size(); j++){
                if(QList[i]->e == edgesToBeUpdated[j]){
                    Eigen::Matrix4f Qmat = getQforEdge(QList[i]->e);
                    Eigen::Vector4f v(QList[i]->e->twin->vert->x, QList[i]->e->twin->vert->y, QList[i]->e->twin->vert->z, 1);
                    QList[i]->Q = v.transpose()*Qmat*v;
                    numUpdates++;
                    break;
                }
            }
            if(numUpdates == edgesToBeUpdated.size())break;
        }
        //        for(it=QList.begin(); it != QList.end();){
        //            if(find(edgesToBeUpdated.begin(), edgesToBeUpdated.end(), (*it)->e) != edgesToBeUpdated.end()){
        //                Eigen::Matrix4f Qmat = getQforEdge((*it)->e);
        //                Eigen::Vector4f v((*it)->e->twin->vert->x, (*it)->e->twin->vert->y, (*it)->e->twin->vert->z, 1);
        //                (*it)->Q = v.transpose()*Qmat*v;
        //            }
        //            ++it;
        //        }
    }
    
    cout<< "\nOutput face count: " << HEF.size()<<endl;
    
    cout << "\nReverting mesh..." << endl;
    revertMesh();
    cout << "\nWriting output file..." << endl;
    writeMF(output);
}

void Mesh::collapseEdge(HEEdge *e){
    e->twin->vert->edge = e->twin->next->twin;
    e->vert->edge = e->next->twin->next;
    
    // Find which edges are going to get deleted
    vector<HEEdge*> edgesToBeDeleted;
    vector<HEEdge*> edgesToBeDeleted1 = Mesh::neighborEdges(e->face);
    vector<HEEdge*> edgesToBeDeleted2 = Mesh::neighborEdges(e->twin->face);
    edgesToBeDeleted.insert(edgesToBeDeleted.end(), edgesToBeDeleted1.begin(), edgesToBeDeleted1.end());
    edgesToBeDeleted.insert(edgesToBeDeleted.end(), edgesToBeDeleted2.begin(), edgesToBeDeleted2.end());
    
    // Find which vertexes will have to be updated
    vector<HEVertex*> vertsToBeUpdated;
    vector<HEVertex*> neighborVerts1 = neighbourVertsFromEdge(e);
    vector<HEVertex*> neighborVerts2 = neighbourVertsFromEdge(e->twin);
    vertsToBeUpdated.insert(vertsToBeUpdated.end(), neighborVerts1.begin(), neighborVerts1.end());
    vertsToBeUpdated.insert(vertsToBeUpdated.end(), neighborVerts2.begin(), neighborVerts2.end());
    //    vector<HEVertex*> neighborVerts1 = neighborVertices(e->vert);
    //    vector<HEVertex*> neighborVerts2 = neighborVertices(e->twin->vert);
    
    // Update twin pairs
    if(DEBUG) cout << "Twin of " << e->next->twin->id;
    e->next->twin->twin = e->next->next->twin;
    if(DEBUG) cout << " set to " << e->next->twin->twin->id << " \n";
    
    if(edgeDeletionMap[e->next->twin]){
        cout << "Twin set to deleted edge " << endl;
    }
    
    if(edgeDeletionMap[e->next->twin->twin]){
        cout << "Twin set to deleted edge " << endl;
    }
    
    if(DEBUG) cout << "Twin of " << e->next->next->twin->id;
    e->next->next->twin->twin = e->next->twin;
    if(DEBUG) cout << " set to " << e->next->next->twin->twin->id << " \n";
    
    if(DEBUG) cout << "Twin of " << e->twin->next->twin->id;
    e->twin->next->twin->twin = e->twin->next->next->twin;
    if(DEBUG) cout << " set to " << e->twin->next->twin->twin->id << " \n";
    
    if(DEBUG) cout << "Twin of " << e->twin->next->next->twin->id;
    e->twin->next->next->twin->twin = e->twin->next->twin;
    if(DEBUG) cout << " set to " << e->twin->next->next->twin->twin->id << " \n";
    
    // Delete Half Edges within the adjacent faces
    vector<HEEdge*>::iterator HEEdgeIterator;
    for(HEEdgeIterator=HEE.begin(); HEEdgeIterator != HEE.end();){
        for(int i = 0; i<edgesToBeDeleted.size(); i++){
            if(find(edgesToBeDeleted.begin(), edgesToBeDeleted.end(), (*HEEdgeIterator)) != edgesToBeDeleted.end()){
                if(DEBUG) cout << "Deleting edge " << (*HEEdgeIterator)->id << " \n";
                edgeDeletionMap[(*HEEdgeIterator)] = true;
                HEE.erase(HEEdgeIterator);
                --HEEdgeIterator;
            }
        }
        ++HEEdgeIterator;
    }
    
    // Delete faces
    int deletedFaces = 0;
    for (int i = 0; i<HEF.size(); i++) {
        if( (HEF[i]->id == e->face->id) || (HEF[i]->id == e->twin->face->id) ){
            faceDeletionMap[HEF[i]] = true;
            HEF.erase(HEF.begin()+i);
            if(DEBUG) cout << "Deleting face " << HEF[i]->id << " \n";
            deletedFaces ++;
            if(deletedFaces == 2){break;}
            i --;
        }
    }
    
    // Remove vertex
    for(int i = 0; i<HEV.size(); i++){
        if(HEV[i]->id == e->vert->id){
            if(DEBUG) cout << "Deleting vertex " << HEV[i]->id << " \n";
            vertexDeletionMap[HEV[i]] = true;
            HEV.erase(HEV.begin()+i);
            break;
        }
    }
    
    // Repoint all edges currently pointing at vertex to be deleted
    HEEdge* nextEdge = e->next->twin;
    do{
        if(DEBUG) cout << "Repointing the outgoing edge of vertex " << nextEdge->twin->vert->id << " to " << nextEdge->id << " \n";
        nextEdge->twin->vert->edge = nextEdge;
        if(DEBUG) cout << "Setting the vertex at the end of edge " << nextEdge->id << " to " << e->twin->vert->id << " \n";
        nextEdge->vert = e->twin->vert;
        nextEdge = nextEdge->next->twin;
    } while(nextEdge->vert->id != e->twin->vert->id);
    
    // Reassign outgoing half edge for neighbourhood vertexes
    for( int i = 0; i < vertsToBeUpdated.size(); i++ ){
        
        // Outgoing edge of vertex has been deleted
        if(edgeDeletionMap[vertsToBeUpdated[i]->edge]){
            
            if(DEBUG)cout << "Ougoing edge of vertex "<< vertsToBeUpdated[i]->id << " has been deleted" << endl;
            
            if(!vertexDeletionMap[vertsToBeUpdated[i]]){
                
                if(DEBUG)cout << "Repointed from " << vertsToBeUpdated[i]->edge->id;
                
                do{
                    vertsToBeUpdated[i]->edge = vertsToBeUpdated[i]->edge->twin->next;
                } while(edgeDeletionMap[vertsToBeUpdated[i]->edge]);
                if(DEBUG)cout << " to " << vertsToBeUpdated[i]->edge->id <<" "<< endl;
                
                if(DEBUG){
                    if(vertsToBeUpdated[i]->edge->twin->vert != vertsToBeUpdated[i]){
                        cout <<"------- RAISE THE ALARMS! ------- \nI am a valid edge who is not the twin of my sibling, meaning the twin of my outgoing edge is not pointing at me! My id is "<< vertsToBeUpdated[i]->id << "\n\n";
                    }
                }
            } else {
                // Vertex has also been deleted
                for(int j = 0; j<HEE.size(); j++){
                    if(HEE[j]->vert == vertsToBeUpdated[i]){
                        if(!edgeDeletionMap[HEE[j]->twin]){
                            vertsToBeUpdated[i]->edge = HEE[j]->twin;
                            break;
                        } else {
                            if(DEBUG){
                                cout << "Twin of valid edge has been deleted... " << endl;
                            }
                        }
                    }
                }
            }
        }
        //            if(vertsToBeUpdated[i] == e->vert){
        //                // Vertex will be deleted
        //                do{
        //                    vertsToBeUpdated[i]->edge = vertsToBeUpdated[i]->edge->twin->next;
        //                } while(edgeDeletionMap[vertsToBeUpdated[i]->edge]);
        //
        //            } else if (vertsToBeUpdated[i] == e->twin->vert){
        //                //
        //                for(int j = 0; j<HEE.size(); j++){
        //                    if(HEE[j]->vert == vertsToBeUpdated[i]){
        //                        if(!edgeDeletionMap[HEE[j]->twin]){
        //                            vertsToBeUpdated[i]->edge = HEE[j]->twin;
        //                            break;
        //                        }
        //                    }
        //                }
        //            } else {
        //
        //                do{
        //                    vertsToBeUpdated[i]->edge = vertsToBeUpdated[i]->edge->next->next->twin;
        //
        //                } while(edgeDeletionMap[vertsToBeUpdated[i]->edge]);
        //            }
        //
        //        }
        //
        //
        ////            do{
        ////                vertsToBeUpdated[i]->edge = vertsToBeUpdated[i]->edge->next->next->twin;
        ////            }while(edgeDeletionMap[vertsToBeUpdated[i]->edge]);
        //
        //        }
        //
        //        if(find(HEE.begin(), HEE.end(), vertsToBeUpdated[i]->edge) != HEE.end()){
        //
        //            // Vertex is not the one which is going to be deleted
        //            if(!(vertsToBeUpdated[i]->id == e->vert->id)){
        //                // Repoint the edge to a valid edge
        //                do{
        //                    vertsToBeUpdated[i]->edge = vertsToBeUpdated[i]->edge->next->next->twin;
        //                }while(edgeDeletionMap[vertsToBeUpdated[i]->edge]);
        ////                } while(!(find(HEE.begin(), HEE.end(), vertsToBeUpdated[i]->edge) != HEE.end()));
        //                // Vertex to be updated is the same vertex that is going to be removed.
        //            } else {
        //
        //                // This vertex is going to be deleted, but for safety-reasons we repoint the edge anyway
        //                for(int j = 0; j<HEE.size(); j++){
        //                    if(HEE[j]->vert == vertsToBeUpdated[i]){
        //                        vertsToBeUpdated[i]->edge = HEE[j]->twin;
        //                        break;
        //                    }
        //                }
        //            }
        //        } else {
        //            for(int j = 0; j<HEE.size(); j++){
        //                if(HEE[j]->vert == vertsToBeUpdated[i]){
        //                    if(find(HEE.begin(), HEE.end(), HEE[j]->twin) != HEE.end()){
        //                        vertsToBeUpdated[i]->edge = HEE[j]->twin;
        //                        break;
        //                    }
        //                }
        //            }
        //        }
        if(DEBUG) {
            if(edgeDeletionMap[vertsToBeUpdated[i]->edge]){
                cout << "--- Reassigned edge not in HEE vector, id "<< vertsToBeUpdated[i]->edge->id << " for vertex "<< vertsToBeUpdated[i]->id <<" \n";
            }
        }
    }
    
    if(DEBUG){
        // Look for broken twin pairs
        for(int i = 0; i<HEE.size(); i++){
            if(HEE[i]->twin->twin->id != HEE[i]->id){
                cout << "--- Broken twin pair between " << HEE[i]->id << " pointing at "<< HEE[i]->vert->id << " and " << HEE[i]->twin->id<< " pointing at "<< HEE[i]->twin->vert->id << " \n";
            }
        }
    }
    
    if(DEBUG){
        // Check if any of the valid edges are pointing to a deleted edge
        for(int i = 0; i<HEV.size(); i++){
            if(edgeDeletionMap[HEV[i]->edge]){
                cout << "--- Outgoing edge of " << HEV[i]->id << " pointing at deleted edge " << HEV[i]->edge->id << " this means that the repointing algorithm isn't working as intended. \n" << endl;
            }
        }
    }
    
    if(DEBUG){
        // Check for illegal twins
        // An illegal twin is an edge who has a twin not present in the HEE vector.
        for(int i = 0; i<HEE.size(); i++){
            if(edgeDeletionMap[HEE[i]->twin]){
                cout << "----- Active edge " << HEE[i]->id << " has a dead twin " << HEE[i]->twin->id << " \n";
            }
        }
    }
    
    if(DEBUG){
        for(int i = 0; i<HEE.size(); i++){
            if(HEE[i]->next->next->next != HEE[i]){
                cout << "---- Broken face with id " << HEE[i]->face->id << endl;
            }
        }
    }
}


void Mesh::convertMesh()
{
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
    
    //    if(DEBUG){
    for(int i = 0; i<HEE.size(); i++){
        if(!(find(HEV.begin(), HEV.end(), HEE[i]->vert) != HEV.end())){
            cout << "Can't find vert of edge " << HEE[i]->id << " which is " << HEE[i]->vert->id << " \n";
        }
    }
    //    }
    
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
        vector<HEVertex*> adjacentVerts = adjacentVertices(HEF[i]);
        
        f.a = locatedAt[adjacentVerts[0]->id];
        f.b = locatedAt[adjacentVerts[1]->id];
        f.c = locatedAt[adjacentVerts[2]->id];
        F.push_back(f);
        
        if(DEBUG){
            if(locatedAt[adjacentVerts[0]->id] == 0){
                cout << "Face " <<HEF[i]->id << " contains vertex " << adjacentVerts[0]->id << " not present in vertex list" << endl;
            }
            
            if(locatedAt[adjacentVerts[1]->id] == 0){
                cout << "Face " <<HEF[i]->id << " contains vertex " << adjacentVerts[1]->id << " not present in vertex list" << endl;
            }
            
            if(locatedAt[adjacentVerts[2]->id] == 0){
                cout << "Face " <<HEF[i]->id << " contains vertex " << adjacentVerts[2]->id << " not present in vertex list" << endl;
            }
        }
    }
}
vector<HEEdge*> Mesh::neighborEdges(HEFace *f){
    vector<HEEdge*> neighborEdges;
    
    if(DEBUG){
        // Check if adjacent edge is a valid half edge, should always be the case
        if(!(find(HEE.begin(), HEE.end(), f->edge) != HEE.end())){
            cout << "------- Edge not present in HEE "<< f->edge->id << " \n";
        }
    }
    
    HEEdge *firstEdge = f->edge;
    HEEdge *nextEdge = firstEdge->next;
    neighborEdges.push_back(nextEdge);
    
    do{
        nextEdge = nextEdge->next;
        neighborEdges.push_back(nextEdge);
    } while(nextEdge->id != firstEdge->id);
    
    return neighborEdges;
}

vector<HEVertex*> Mesh::neighborVertices(HEVertex *v)
{
    vector<HEVertex*> neighbourVerts;
    
    HEVertex *firstVertex = v->edge->vert;
    HEEdge *nextEdge = v->edge;
    
    do{
        nextEdge = nextEdge->twin->next;
        neighbourVerts.push_back(nextEdge->vert);
    } while(nextEdge->vert != firstVertex);
    
    return neighbourVerts;
}

vector<HEEdge*> Mesh::associatedEdges(HEVertex *v){
    vector<HEEdge*> edges;
    
    HEEdge* firstEdge = v->edge;
    HEEdge *nextEdge = firstEdge->twin->next;
    edges.push_back(nextEdge);
    
    do{
        nextEdge = nextEdge->twin->next;
        edges.push_back(nextEdge);
    } while(nextEdge->id != firstEdge->id);
    
    return edges;
}


vector<HEEdge*> Mesh::edgesPointingToSameVert(HEEdge *e){
    vector<HEEdge*> edges;
    
    HEEdge* firstEdge = e;
    HEEdge *nextEdge = firstEdge->next->twin;
    edges.push_back(nextEdge);
    
    do{
        nextEdge = nextEdge->next->twin;
        edges.push_back(nextEdge);
    } while(nextEdge->id != firstEdge->id);
    
    return edges;
}

vector<HEFace*> Mesh::neighborFaces(HEVertex *v)
{
    vector<HEFace*> adjacentFaces;
    
    HEEdge *firstEdge = v->edge;
    adjacentFaces.push_back(firstEdge->face);
    
    HEEdge *nextEdge = firstEdge->twin->next;
    do{
        adjacentFaces.push_back(nextEdge->face);
        nextEdge = nextEdge->twin->next;
    } while(nextEdge != firstEdge);
    
    return adjacentFaces;
}

vector<HEVertex*> Mesh::adjacentVertices(HEFace *f)
{
    vector<HEVertex*> adjacentVerts;
    
    HEEdge *firstEdge = f->edge;
    HEEdge *nextEdge = f->edge->next;
    
    adjacentVerts.push_back(firstEdge->vert);
    
    do{
        adjacentVerts.push_back(nextEdge->vert);
        nextEdge = nextEdge->next;
    } while(nextEdge != firstEdge);
    
    return adjacentVerts;
}

vector<HEVertex*> Mesh::neighbourVertsFromEdge(HEEdge* e){
    vector<HEVertex*> neighbourVerts;
    
    HEEdge *nextEdge = e;
    do{
        nextEdge = nextEdge->twin->next;
        neighbourVerts.push_back(nextEdge->vert);
    } while(nextEdge != e);
    
    return neighbourVerts;
}

int Mesh::Vcnt(){
    /* Half Edge Vertex */
    return (int)HEV.size();
}

int Mesh::Fcnt(){
    /* Half Edge Face */
    return (int)HEF.size();
}