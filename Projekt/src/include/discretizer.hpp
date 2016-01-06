#ifndef __DISCRETIZER_HPP_
#define __DISCRETIZER_HPP_

#include<string>
#include<sstream>
#include<fstream>
#include<vector>

namespace Icarus
{

//*** Vertex ***//

struct Vertex {
    float x, y, z;
};

//*** Face ***//

class Face
{
public:
    Face(int num_vertices, std::vector<int>& id_vertices,
        std::vector<Vertex>& vertices, Vertex normal);
    
    bool pointInsideYz(Vertex point);
    
    Vertex get_vertex(int local_id) { return _vertices[local_id]; }
    
    Vertex get_normal() { return _normal; }
    
private:
    int _num_vertices;
    std::vector<Vertex> _vertices;
    Vertex _normal;
    Vertex _last_point;
    bool _last_check;
};

//*** Object ***//

class Object
{
public:
    Object(std::string name);
    
    void set_vertex(float x, float y, float z);
    
    void set_normal(float x, float y, float z);
    
    void set_face(int num_vertices, std::vector<int>& id_vertices, int local_id_normal);
    
	/**
	* \returns Gibt 'o' zurück, wenn sich der Punkt im Objekt befindet, sonst 'a'.
	*/
    char pointInside(Vertex point);
    
private:
    std::string _name;
    int _num_vertices;
    int _num_faces;
    std::vector<Vertex> _vertices;
    std::vector<Vertex> _normals;
    std::vector<Face> _faces;
};

/***
 * \returns Gibt für jeden Gitterpunkt 'a' (Frei), 'o' (innerhalb von Objekt)
 *			oder 'b' (Rand des Gesamtgebietes) zurück. Die Indizes werden von
 *			innen nach außen als x,y,z linear abgerollt.
 */
std::vector<char> discretizer(std::string inputFile,
    float h, int nx, int ny, int nz);

void save_discretizer(std::vector<char> discretized_points,
                      std::string outputFile,
                      int nx, int ny, int nz);

}//namespace Icarus

#endif//__DISCRETIZER_HPP_
