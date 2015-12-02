// -std=c++11 noetig

#ifndef __DISCRETIZER_HPP_
#define __DISCRETIZER_HPP_

#include<string>
#include<sstream>
#include<fstream>
#include<vector>

namespace Icarus
{

struct Vertex {
    float x, y, z;
};

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

class Object
{
public:
    Object(std::string name);
    
    void set_vertex(float x, float y, float z);
    
    void set_normal(float x, float y, float z);
    
    void set_face(int num_vertices, std::vector<int>& id_vertices, int local_id_normal);
    
    char pointInside(Vertex point);
    
private:
    std::string _name;
    int _num_vertices;
    int _num_faces;
    std::vector<Vertex> _vertices;
    std::vector<Vertex> _normals;
    std::vector<Face> _faces;
};

}//namespace Icarus

#endif//__DISCRETIZER_HPP_
