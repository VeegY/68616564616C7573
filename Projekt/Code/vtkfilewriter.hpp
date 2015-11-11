// #################################################################################################
//			Studienprojekt Modellbildung & Simulation - 2015/16
// #################################################################################################
// 					Header: vtkFileWriter.hpp
// ------------------------------------Doxygen-Dokumentation----------------------------------------
///  \file vtkFileWriter.hpp
///  \brief
///  Stellt Klasse fuer die Ausgabe der Ergebnisse in eine .vtk Datei bereit.
///
//#################################################################################################


#include "vector.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

// ========================================DOKUMENTATION============================================
///\fn		template<typename type> string toString(type val)
///\brief	Konvertiert verschiedene Datentypen zu einem String.
///         Alternative, falls to_String(...) nicht benutzt werden kann.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp des zu konvertierenden Objekts.
///@param	val		Objekt das Konvertiert werden soll.
// =================================================================================================
template<typename type>
std::string toString(const type val)
{
    stringstream ss;
    ss << val;
    string str = ss.str();
    return str;
}

// ========================================DOKUMENTATION============================================
///\class 	vtkFileWriter
///\brief 	Erstellt .vtk Dateien und fuellt diese mit Daten.
///         Die Dateien sind vom Typ structuerd point
// =================================================================================================
class vtkFileWriter
{
        bool* cellDataWrittenLast;
        bool* pointDataWrittenLast;
        string filename;
        string title;
        unsigned int tsteps;
        int xdim, ydim, zdim;
        int cells, points;

    public:


// ========================================DOKUMENTATION============================================
///\brief	Konstruktor
// -------------------------------------------------------------------------------------------------
///@param 	_filename       string mit dem Dateinamen OHNE DATEIENDUNG
///@param   _title          Titel der Datei
///@param   _xdim           Anzahl der Punkte in x Richtung
///@param   _ydim           Anzahl der Punkte in y Richtung
///@param   _zdim           Anzahl der Punkte in z Richtung
///@param   timesteps       Anzahl der Zeitschritte
// -------------------------------------------------------------------------------------------------
        vtkFileWriter(string _filename, string _title, int _xdim, int _ydim, int _zdim, unsigned int timesteps):
            filename(_filename),
            xdim(_xdim),
            ydim(_ydim),
            zdim(_zdim),
            tsteps(timesteps),
            title(_title)
            {
                cells=((xdim-1)*(ydim-1)*(zdim-1));
                points=(xdim*ydim*zdim);
                cellDataWrittenLast=new bool[tsteps]();
                pointDataWrittenLast=new bool[tsteps]();
                ofstream file;
                string hfname(filename);
                if (tsteps <= 1){
                    hfname.append(".vtk");
                    file.open(hfname.c_str(), ios::out | ios::trunc );
                    if (file.is_open())
                    {
                        file << "# vtk DataFile Version 3.0\n";
                        file << title << endl;
                        file << "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS ";
                        file << xdim << " " << ydim << " " << zdim << endl;
                        file << "ASPECT_RATIO 1 1 1\nORIGIN 0 0 0\n\n";
                        file.close();
                    }
                    else cout << "Unable to open file";
                }else
                {
                  for (int i=1; i<=tsteps; ++i){
                    hfname=filename;
                    hfname.append(".vtk.");
                    hfname.append(toString(i));                                //falls to_string nicht läuft mit tostring(i)
                    file.open(hfname.c_str(), ios::out | ios::trunc );
                    if (file.is_open())
                    {
                        file << "# vtk DataFile Version 3.0\n";
                        file << title << endl;
                        file << "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS ";
                        file << xdim << " " << ydim << " " << zdim << endl;
                        file << "ASPECT_RATIO 1 1 1\nORIGIN 0 0 0\n\n";
                        file.close();
                    }
                    else cout << "Unable to open file";
                  }
                }
            };
// =================================================================================================


// ========================================DOKUMENTATION============================================
//					    Destruktor
///\brief	Destruktor
// -------------------------------------------------------------------------------------------------
        ~vtkFileWriter()
        {
            delete cellDataWrittenLast;
            delete pointDataWrittenLast;
        };
// =================================================================================================

// ========================================DOKUMENTATION============================================
///\brief	Fuegt skalarwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	data[]	Array der Quelldaten
///@param   length  Länge des Arrays
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addPointDataToTimestep(const type data[], const int length, const int timestep, string name)
        {
            if (length==points){
                ofstream file;
                string hfname(filename);
                hfname=filename;
                hfname.append(".vtk.");
                hfname.append(toString(timestep+1));                                //falls to_string nicht läuft mit tostring(i)
                file.open(hfname.c_str(), ios::out | ios::app);
                if (file.is_open())
                {
                    if (pointDataWrittenLast[timestep]==false)
                    {
                        file << "POINT_DATA " << points <<endl;
                    }
                    file << "SCALARS " << name <<" float"<<endl;
                    file << "LOOKUP_TABLE default" << endl;
                    for (int i=0; i<points; ++i)
                    {
                        file << static_cast<float>(data[i])<<endl;
                    }
                    file << endl;
                    file.close();
                    pointDataWrittenLast[timestep]=true;
                    cellDataWrittenLast[timestep]=false;
                }

                else cout << "Unable to open file";
            }else
            {
                cout<<"false dimension: no data written!"<<endl;
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt skalarwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	data    Vektor mit den zu schreibenden Daten
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addPointDataToTimestep(const Vector<type>& data,const int timestep, const string name)
        {
            if(data.dim()==points)
            {
                ofstream file;
                string hfname(filename);
                hfname=filename;
                hfname.append(".vtk.");
                hfname.append(toString(timestep+1));                                //falls to_string nicht läuft mit tostring(i)
                file.open(hfname.c_str(), ios::out | ios::app);
                if (file.is_open())
                {
                    if (pointDataWrittenLast[timestep]==false)
                    {
                        file << "POINT_DATA " << points <<endl;
                    }
                    file << "SCALARS " << name <<" float"<<endl;
                    file << "LOOKUP_TABLE default" << endl;
                    for (int i=0; i<points; ++i)
                    {
                        file << static_cast<float>(data[i])<<endl;
                    }
                    file << endl;
                    file.close();
                    pointDataWrittenLast[timestep]=true;
                    cellDataWrittenLast[timestep]=false;
                }else cout << "Unable to open file";
            }else
            {
                cout<<"false dimension: no data written!"<<endl;
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt skalarwertige Zelldaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	data[]	Array der Quelldaten
///@param   length  Länge des Arrays
///@param   timestep Zeitschritt, dem die Daten hinzugefügt werden
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addCellDataToTimestep(const type* data[], const int length, const int timestep, const string name)
        {
            if (length==cells){
                ofstream file;
                string hfname(filename);
                hfname=filename;
                hfname.append(".vtk.");
                hfname.append(toString(timestep+1));                                //falls to_string nicht läuft mit tostring(i)
                file.open(hfname.c_str(), ios::out | ios::app);
                if (file.is_open())
                {
                    if (cellDataWrittenLast[timestep]==false)
                    {
                        file << "CELL_DATA " << cells <<endl;
                    }
                    file << "SCALARS " << name <<" float"<<endl;
                    file << "LOOKUP_TABLE default" << endl;
                    for (int i=0; i<cells; ++i)
                    {
                        file << static_cast<float>(data[i])<<endl;
                    }
                    file << endl;
                    file.close();
                    pointDataWrittenLast[timestep]=false;
                    cellDataWrittenLast[timestep]=true;
                }

                else cout << "Unable to open file";
            }else
            {
                cout<<"false dimension: no data written!"<<endl;
            }

        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt skalarwertige Zelldaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	data    Vektor mit den zu schreibenden Daten
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addCellDataToTimestep(const Vector<type>& data, const int timestep, const string name)
        {
            if(data.dim()==cells)
            {
                ofstream file;
                string hfname(filename);
                hfname=filename;
                hfname.append(".vtk.");
                hfname.append(toString(timestep+1));                                //falls to_string nicht läuft mit tostring(i)
                file.open(hfname.c_str(), ios::out | ios::app);
                if (file.is_open())
                {
                    if (cellDataWrittenLast[timestep]==false)
                    {
                        file << "CELL_DATA " << cells <<endl;
                    }
                    file << "SCALARS " << name <<" float"<<endl;
                    file << "LOOKUP_TABLE default" << endl;
                    for (int i=0; i<cells; ++i)
                    {
                        file << static_cast<float>(data[i])<<endl;
                    }
                    file << endl;
                    file.close();
                    pointDataWrittenLast[timestep]=false;
                    cellDataWrittenLast[timestep]=true;
                }else cout << "Unable to open file";
            }else
            {
                cout<<"false dimension: no data written!"<<endl;
            }

        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt skalarwertige Punktdaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	data[]	Array der Quelldaten
///@param   length  Länge des Arrays
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addPointDataToAll(const type data[], int length, string name)
        {
            for (int i=0; i<tsteps; ++i){
                addPointDataToTimestep(data, length, i, name);
            }
        }


// ========================================DOKUMENTATION============================================
///\brief	Fuegt skalarwertige Punktdaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	data    Vektor mit den zu schreibenden Daten
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addPointDataToAll(const Vector<type>& data,const string name)
        {
            for (int i=0; i<tsteps; ++i){
                addPointDataToTimestep(data, i, name);
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt skalarwertige Zelldaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	data[]	Array der Quelldaten
///@param   length  Länge des Arrays
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addCellDataToAll(const type data[], int length, string name)
        {
            for (int i=0; i<tsteps; ++i){
                addCellDataToTimestep(data, length, i, name);
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt skalarwertige Zelldaten zu allen Zeitschritten hinzu
///         Die Daten werden in der Datei als float abgespeichert
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	data    Vektor mit den zu schreibenden Daten
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addCellDataToAll(const Vector<type>& data, const string name)
        {
            for (int i=0; i<tsteps; ++i){
                addCellDataToTimestep(data, i, name);
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt vektorwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	datax[]	Array der x Komponenten
///@param	datay[]	Array der y Komponenten
///@param	dataz[]	Array der z Komponenten
///@param   length  Länge der Arrays
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addPointVecToTimestep(const type datax[], const type datay[], const type dataz[], const int length, const int timestep, string name)
        {
            if (length==points){
                ofstream file;
                string hfname(filename);
                hfname=filename;
                hfname.append(".vtk.");
                hfname.append(toString(timestep+1));                                //falls to_string nicht läuft mit tostring(i)
                file.open(hfname.c_str(), ios::out | ios::app);
                if (file.is_open())
                {
                    if (pointDataWrittenLast[timestep]==false)
                    {
                        file << "POINT_DATA " << points <<endl;
                    }
                    file << "VECTORS " << name <<" float"<<endl;

                    for (int i=0; i<points; ++i)
                    {
                        file << static_cast<float>(datax[i])<<" "<< static_cast<float>(datay[i])<<" "<< static_cast<float>(dataz[i])<<endl;
                    }
                    file << endl;
                    file.close();
                    pointDataWrittenLast[timestep]=true;
                    cellDataWrittenLast[timestep]=false;
                }

                else cout << "Unable to open file";
            }else
            {
                cout<<"false dimension: no data written!"<<endl;
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt vektorwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	    Datentyp der Quelldaten
///@param	datax       Vektor mit den x Komponenten der zu schreibenden  Daten
///@param	datay       Vektor mit den y Komponenten der zu schreibenden  Daten
///@param	dataz       Vektor mit den z Komponenten der zu schreibenden  Daten
///@param   timestep    Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name        Name der Daten
// =================================================================================================
        template<typename type>
        void addPointVecToTimestep(const Vector<type>& datax, const Vector<type>& datay, const Vector<type>& dataz, const int timestep, const string name)
        {
            if(datax.dim()==points && datay.dim()==points && dataz.dim()==points)
            {
                ofstream file;
                string hfname(filename);
                hfname=filename;
                hfname.append(".vtk.");
                hfname.append(toString(timestep+1));                                //falls to_string nicht läuft mit tostring(i)
                file.open(hfname.c_str(), ios::out | ios::app);
                if (file.is_open())
                {
                    if (pointDataWrittenLast[timestep]==false)
                    {
                        file << "POINT_DATA " << points <<endl;
                    }
                    file << "VECTORS " << name <<" float"<<endl;
                    for (int i=0; i<points; ++i)
                    {
                        file << static_cast<float>(datax[i]) << " " << static_cast<float>(datay[i]) << " " << static_cast<float>(dataz[i])<<endl;
                    }
                    file << endl;
                    file.close();
                    pointDataWrittenLast[timestep]=true;
                    cellDataWrittenLast[timestep]=false;
                }else cout << "Unable to open file";
            }else
            {
                cout<<"false dimension: no data written!"<<endl;
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt vektorwertige Zelldaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	datax[]	Array der x Komponenten
///@param	datay[]	Array der y Komponenten
///@param	dataz[]	Array der z Komponenten
///@param   length  Länge der Arrays
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addCellVecToTimestep(const type datax[], const type datay[], const type dataz[], const int length, const int timestep, string name)
        {
            if (length==cells){
                ofstream file;
                string hfname(filename);
                hfname=filename;
                hfname.append(".vtk.");
                hfname.append(toString(timestep+1));                                //falls to_string nicht läuft mit tostring(i)
                file.open(hfname.c_str(), ios::out | ios::app);
                if (file.is_open())
                {
                    if (cellDataWrittenLast[timestep]==false)
                    {
                        file << "CELL_DATA " << cells <<endl;
                    }
                    file << "VECTORS " << name <<" float"<<endl;
                    for (int i=0; i<cells; ++i)
                    {
                        file << static_cast<float>(datax[i])<<" "<< static_cast<float>(datay[i])<<" "<< static_cast<float>(dataz[i])<<endl;
                    }
                    file << endl;
                    file.close();
                    cellDataWrittenLast[timestep]=true;
                    pointDataWrittenLast[timestep]=false;
                }

                else cout << "Unable to open file";
            }else
            {
                cout<<"false dimension: no data written!"<<endl;
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt vektorwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	    Datentyp der Quelldaten
///@param	datax       Vektor mit den x Komponenten der zu schreibenden  Daten
///@param	datay       Vektor mit den y Komponenten der zu schreibenden  Daten
///@param	dataz       Vektor mit den z Komponenten der zu schreibenden  Daten
///@param   timestep    Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name        Name der Daten
// =================================================================================================
        template<typename type>
        void addCellVecToTimestep(const Vector<type>& datax, const Vector<type>& datay, const Vector<type>& dataz, const int timestep, const string name)
        {
            if(datax.dim()==cells && datay.dim()==cells && dataz.dim()==cells)
            {
                ofstream file;
                string hfname(filename);
                hfname=filename;
                hfname.append(".vtk.");
                hfname.append(toString(timestep+1));                                //falls to_string nicht läuft mit tostring(i)
                file.open(hfname.c_str(), ios::out | ios::app);
                if (file.is_open())
                {
                    if (cellDataWrittenLast[timestep]==false)
                    {
                        file << "CELL_DATA " << cells <<endl;
                    }
                    file << "VECTORS " << name <<" float"<<endl;
                    for (int i=0; i<cells; ++i)
                    {
                        file << static_cast<float>(datax[i]) << " " << static_cast<float>(datay[i]) << " " << static_cast<float>(dataz[i])<<endl;
                    }
                    file << endl;
                    file.close();
                    cellDataWrittenLast[timestep]=true;
                    pointDataWrittenLast[timestep]=false;
                }else cout << "Unable to open file";
            }else
            {
                cout<<"false dimension: no data written!"<<endl;
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt vektorwertige Punktdaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	datax[]	Array der x Komponenten
///@param	datay[]	Array der y Komponenten
///@param	dataz[]	Array der z Komponenten
///@param   length  Länge der Arrays
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addPointVecToAll(const type datax[], const type datay[], const type dataz[], int length, string name)
        {
            for (int i=0; i<tsteps; ++i){
                addPointVecToTimestep(datax, datay, dataz, length, i, name);
            }
        }

// ========================================DOKUMENTATION============================================
///\brief	Fuegt vektorwertige Punktdaten zu allen Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	    Datentyp der Quelldaten
///@param	datax       Vektor mit den x Komponenten der zu schreibenden  Daten
///@param	datay       Vektor mit den y Komponenten der zu schreibenden  Daten
///@param	dataz       Vektor mit den z Komponenten der zu schreibenden  Daten
///@param   name        Name der Daten
// =================================================================================================
        template<typename type>
        void addPointVecToAll(const Vector<type>& datax, const Vector<type>& datay, const Vector<type>& dataz, const string name)
        {
            for (int i=0; i<tsteps; ++i){
                addPointCellToTimestep(datax, datay, dataz, i, name);
            }
        }


// ========================================DOKUMENTATION============================================
///\brief	Fuegt vektorwertige Zelldaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	Datentyp der Quelldaten
///@param	datax[]	Array der x Komponenten
///@param	datay[]	Array der y Komponenten
///@param	dataz[]	Array der z Komponenten
///@param   length  Länge der Arrays
///@param   name    Name der Daten
// =================================================================================================
        template<typename type>
        void addCellVecToAll(const type datax[], const type datay[], const type dataz[], int length, string name)
        {
            for (int i=0; i<tsteps; ++i){
                addCellVecToTimestep(datax, datay, dataz, length, i, name);
            }
        }


// ========================================DOKUMENTATION============================================
///\brief	Fuegt vektorwertige Zelldaten zu allen Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam	type 	    Datentyp der Quelldaten
///@param	datax       Vektor mit den x Komponenten der zu schreibenden  Daten
///@param	datay       Vektor mit den y Komponenten der zu schreibenden  Daten
///@param	dataz       Vektor mit den z Komponenten der zu schreibenden  Daten
///@param   name        Name der Daten
// =================================================================================================
        template<typename type>
        void addCellVecToAll(const Vector<type>& datax, const Vector<type>& datay, const Vector<type>& dataz, const string name)
        {
            for (int i=0; i<tsteps; ++i){
                addCellVecToTimestep(datax, datay, dataz, i, name);
            }
        }

};


