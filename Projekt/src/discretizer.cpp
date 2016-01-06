#ifndef __DISCRETIZER_CPP_
#define __DISCRETIZER_CPP_

#include "include/discretizer.hpp"

namespace Icarus
{

//****************  Face  ****************//

Face::Face(int num_vertices, std::vector<int>& id_vertices,
        std::vector<Vertex>& vertices, Vertex normal):
    _num_vertices(num_vertices),
    _normal(normal),
    _last_point(Vertex{1000.0f, 1000.0f, 1000.0f}), // evtl gibt es ja noch eine schlauere Loesung...
    _last_check(false)
{
    for (int v(0); v<_num_vertices; ++v)
    {
        _vertices.push_back(vertices[id_vertices[v]]);
    }
    _vertices.push_back(_vertices[0]);    // dann ist ***[_num_vertices]=***[0], macht einige Zugriffe einfacher
}


bool Face::pointInsideYz(Vertex point)
{
    if (_last_point.y == point.y && _last_point.z == point.z)
    {
        return _last_check;
    }
    // wuerde hier ein 'else' das ganze uebersichtlicher machen?
    // ich mag es nach dem Testen so viel wie moeglich zu kuerzen :)
    bool inside(false);
    _last_point = point;
    for (int l(0); l<_num_vertices; ++l)
    {
        if (((_vertices[l].z > point.z) != (_vertices[l+1].z > point.z))
            || ((_vertices[l].z >= point.z) != (_vertices[l+1].z >= point.z)))
        {
            float temp((point.z - _vertices[l].z) / (_vertices[l+1].z - _vertices[l].z)
                * (_vertices[l+1].y - _vertices[l].y) + _vertices[l].y);
            if (temp > point.y-1e-6f && temp < point.y+1e-6f)
                return _last_check = true;
            else if (temp > point.y)
                inside = !inside;
        }
    }
    return _last_check = inside;
}

//****************  Object  ****************//

Object::Object(std::string name):
    _name(name),
    _num_vertices(0),
    _num_faces(0)
{
}

void Object::set_vertex(float x, float y, float z)
{
    _vertices.push_back(Vertex{x, y, z});
    ++_num_vertices;
}

void Object::set_normal(float x, float y, float z)
{
    _normals.push_back(Vertex{x, y, z});
}

void Object::set_face(int num_vertices, std::vector<int>& id_vertices, int local_id_normal)
{
    _faces.push_back(Face(num_vertices, id_vertices, _vertices, _normals[local_id_normal]));
    ++_num_faces;
}

char Object::pointInside(Vertex point)
{
    int in_and_out(0);
    for (int f(0); f<_num_faces; ++f)
    {
        if (_faces[f].pointInsideYz(point))
        {
            /*
            Vertex vector{_faces[f].get_vertex(0).x - point.x,
                _faces[f].get_vertex(0).y - point.y,
                _faces[f].get_vertex(0).z - point.z};
            float scalprod(vector.x * _faces[f].get_normal().x +
                vector.y * _faces[f].get_normal().y + vector.z * _faces[f].get_normal().z);
            */
            float scalprod((_faces[f].get_vertex(0).x - point.x) * _faces[f].get_normal().x
                + (_faces[f].get_vertex(0).y - point.y) * _faces[f].get_normal().y
                + (_faces[f].get_vertex(0).z - point.z) * _faces[f].get_normal().z);
            if (scalprod > 1e-6f)
            {
                ++in_and_out;
            }
            else if (scalprod < -1e-6f)
            {
                --in_and_out;
            }
        }
    }
    if (in_and_out != 0)
    {
        return 'o';
    }
    return 'a';
}

//****************  discretizer  ****************//

std::vector<char>discretizer(std::string inputFile,
    float h, int nx, int ny, int nz)
// TODO: "von wo bis wo" einstellen koennen
{
    //*** Ein- und Ausgabedatei oeffnen ***//
    
    std::ifstream fin;
    fin.open(inputFile);
    // TODISCUSS: Dateien lieber nur so lange wie noetig oeffnen?
    // dh: fin oeffnen, lesen, schliessen - rechnen - fout oeffnen, schreiben, schliessen
    // dadurch waeren waehrend der Rechnung alle Dateien geschlossen.


    //*** Daten einlesen ***//

    std::vector<Object> objects; // Vektor fuer alle Objekte
    
    // diverse Variablen, die zum Einlesen benoetigt werden
    std::string line;
    std::string type;
    std::stringstream stream;
    std::string name;
    float x, y, z;
    int num_objects(0), num_vertices(0), num_faces(0), num_normals(0);
    int num_vertices_per_face(0), num_slashes(0), id(0), diff_vertices_id(-1), diff_normals_id(-1), id_normal(0);
    std::vector<int> vertices_ids_of_face;
    
    // das ganze Dokument wird Zeile fuer Zeile durchgelesen
    // Format/ Aufbau einer .obj Datei von Blender, siehe: http://www.martinreddy.net/gfx/3d/OBJ.spec
    while (getline(fin, line))
    {
        stream.str(line);   // ein Stringstream vereinfacht das Einlesen der Daten einer Zeile wesentlich
        stream >> type;
        
        if (type == "o")    // neues Objekt
        {
            stream >> name;
            objects.push_back(Object(name));    // add new object to the list
            ++num_objects;
            diff_vertices_id = num_vertices + 1;    // needed to calculate the local vertex id
            diff_normals_id = num_normals + 1;    // needed to calculate the local normal id
        }//== "o"
        
        else if (type == "v")   // neuer Eckpunkt von aktuellem Objekt
        {
            stream >> x >> y >> z;
            objects[num_objects-1].set_vertex(x, y, z);
            ++num_vertices;
        }//== "v"
        
        else if (type == "vn")  // neuer Normalenvektor von aktuellem Objekt
        {
            stream >> x >> y >> z;
            objects[num_objects-1].set_normal(x, y, z);
            ++num_normals;
        }//== "vn"
        
        else if (type == "f")   // neue Flaeche von aktuellem Objekt (nur Polygone werden sinnvoll eingelesen)
        {
            // erst wird der Aufbau der Zeile untersucht
            stream >> line;
            if (line.find_first_of('/') != line.find_last_of('/') &&
                line.find_last_of('/') != std::string::npos)
            {
                num_slashes = 2;
            }
            else if (line.find_first_of('/') != std::string::npos)
            {
                num_slashes = 1;
            }
            else
            {
                num_slashes = 0;
            }
            
            // Wenn noch kein Normalenvektor existiert, muss er noch berechnet werden (TODO)
            if (line.find_last_of('/') + 1 < line.length())
            {
                stream.seekg(2);
                getline(stream, line, '/');
                if (num_slashes == 2)
                    getline(stream, line, '/');
                stream >> id_normal;
                id_normal -= diff_normals_id;
            }
            
            // zugehoerige Eckpunkte werden eingelesen
            vertices_ids_of_face.erase(vertices_ids_of_face.begin(), vertices_ids_of_face.end());
            stream.seekg(2);    // first char is f, then space, then first vertex-id
            num_vertices_per_face = 0;
            while (stream >> id)
            {
                vertices_ids_of_face.push_back(id - diff_vertices_id);  // Vertex-ID einlesen
                getline(stream, line, ' '); // der Rest an Information des Punktes wird nicht benoetigt
                ++num_vertices_per_face;
            }
            
            // Flache wird dem aktuellen Objekt hinzugefuegt
            objects[num_objects-1].set_face(num_vertices_per_face, vertices_ids_of_face, id_normal);
            ++num_faces;
        }//== "f"
        
        /*else if (type == "usemtl")
        * {
        *     // TODO: spaeter, wenn auch Waermeleitung in Festkoerpern betrachtet wird
        * }//== "usemtl"
        */
        
        stream.clear();
    }//while(getline)

    //*** Raum diskretisieren ***//

    std::vector<char> discretized_points(nx*ny*nz);
    // Preufe fuer jeden Punkt, ob Luft oder Gegenstand
    
    //{
    for (int z(0); z<nz; ++z)
    {
        for (int y(0); y<ny; ++y)
        {
            for (int x(0); x<nx; ++x)
            {
                Vertex point{(float)x*h, (float)y*h, (float)z*h};
                char what('a');
				// prüfe, ob punkt in irgendeinem objekt liegt
                for (int o(0); o<num_objects && what=='a'; ++o)
                {
                    what = objects[o].pointInside(point);
                }
                discretized_points.push_back(what);
            }//x-loop
        }//y-loop
    }//z-loop
    
    // Pruefe fuer alle Objekt-Punkte, ob es Randpunkte sind
    for (int z(0); z<nz; ++z)
    {
        for (int y(0); y<ny; ++y)
        {
            for (int x(0); x<nx; ++x)
            {
                if (discretized_points[z*ny*nx + y*nx + x] != 'a')
                {
					// globaler rand des gebiets
                    if (x==0 || y==0 || z==0 || x==nx-1 || y==ny-1 || z==nz-1)
                    {
                        // TODO
                        discretized_points[z*ny*nx + y*nx + x] = 'b';
                    }
					// wenn min. ein nachbar frei und ich nicht, bin ich rand
                    else if ((discretized_points[z*ny*nx + y*nx + x + 1] == 'a')
                    || (discretized_points[z*ny*nx + y*nx + x - 1] == 'a')
                    || (discretized_points[z*ny*nx + (y+1)*nx + x] == 'a')
                    || (discretized_points[z*ny*nx + (y-1)*nx + x] == 'a')
                    || (discretized_points[(z+1)*ny*nx + y*nx + x] == 'a')
                    || (discretized_points[(z-1)*ny*nx + y*nx + x] == 'a'))
                    {
                        discretized_points[z*ny*nx + y*nx + x] = 'b';
                    }
                }//!= 'a'
            }//x-loop
        }//y-loop
    }//z-loop
    //*/
    
    // hier ist discretized_points fertig
    return discretized_points;
}

void save_discretizer(std::vector<char> discretized_points,
                      std::string outputFile,
                      float h,
                      int nx, int ny, int nz)
{
    std::ofstream fout(outputFile);

    //*** Diskretisierung abspeichern ***//
    fout << nx << " " << ny << " " << nz << std::endl;
    
    /*
    for (int i(0); i<nx*ny*nz; ++i)
    {
        //fout << discretized_points[i] << " ";
        fout << discretized_points[i] << std::endl;
        // TODISCUSS: direkt hintereinander, Leerzeichen oder neue Zeile?
    }
    */
    
    ///*
    for (int z(0); z<nz; ++z)
    {   fout << "z = " << z*h << std::endl;
        for (int y(0); y<ny; ++y)
        {
            for (int x(0); x<nx; ++x)
            {
                fout << ((discretized_points[z*ny*nx + y*nx + x] == 'o')? 'o' : ((discretized_points[z*ny*nx + y*nx + x] == 'b')? 'b' : '.'));
                //fout << discretized_points[z*ny*nx + y*nx + x];
            }//z-loop
            fout << std::endl;
        }//y-loop
        fout << std::endl << std::endl;
    }//x-loop
    //*/
 }

}//namespace Icarus

#endif//__DISCRETIZER_CPP_
