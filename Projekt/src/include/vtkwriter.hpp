// #################################################################################################
//     Studienprojekt Modellbildung & Simulation - 2015/16
// #################################################################################################
//     Header: vtkwriter.hpp
// ------------------------------------Doxygen-Dokumentation----------------------------------------
///  \file VtkFileWriter.hpp
///  \brief
///  Stellt Klasse fuer die Ausgabe der Ergebnisse in eine .vtk Datei bereit.
///
//#################################################################################################

#ifndef __VTKWRITER_HPP_
#define __VTKWRITER_HPP_

#include "logger.hpp"
#include "fullvector.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

namespace Icarus
{
class vtkWriter
{
    bool* _cell_data_written_last;
    bool* _point_data_written_last;
    std::string _filename;
    std::string _title;
    size_t _tsteps;
    size_t _xdim, _ydim, _zdim;
    size_t _cells, _points;

public:

///\brief   Konstruktor
// -------------------------------------------------------------------------------------------------
///@param   filename       string mit dem Dateinamen OHNE DATEIENDUNG
///@param   title          Titel der Datei
///@param   xdim           Anzahl der Punkte in x Richtung
///@param   ydim           Anzahl der Punkte in y Richtung
///@param   zdim           Anzahl der Punkte in z Richtung
///@param   timesteps       Anzahl der Zeitschritte
    vtkWriter(std::string filename, std::string title, size_t xdim, size_t ydim, size_t zdim, size_t timesteps);

///\brief   Konstruktor
// -------------------------------------------------------------------------------------------------
///@param   filename       string mit dem Dateinamen OHNE DATEIENDUNG
///@param   title          Titel der Datei
///@param   xdim           Anzahl der Punkte in x Richtung
///@param   ydim           Anzahl der Punkte in y Richtung
///@param   zdim           Anzahl der Punkte in z Richtung
///@param   h              Ortschrittweite in alle Richtungen
///@param   timesteps       Anzahl der Zeitschritte
    vtkWriter(std::string filename, std::string title, size_t xdim, size_t ydim, size_t zdim, double h, size_t timesteps);

///\brief   Destruktor
    ~vtkWriter();

///\brief   Fuegt skalarwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   data[]  Array der Quelldaten
///@param   length  Länge des Arrays
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten

    template<typename type>
    void addPointDataToTimestep(const type data[], const size_t length, const size_t timestep, std::string name);

///\brief   Fuegt skalarwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   data    Vektor mit den zu schreibenden Daten
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten

    template<typename type>
    void addPointDataToTimestep(const FullVector<type>& data, const size_t timestep, const std::string name);


///\brief   Fuegt skalarwertige Zelldaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   data[]  Array der Quelldaten
///@param   length  Länge des Arrays
///@param   timestep Zeitschritt, dem die Daten hinzugefügt werden
///@param   name    Name der Daten

    template<typename type>
    void addCellDataToTimestep(const type data[], const size_t length, const size_t timestep, const std::string name);


///\brief   Fuegt skalarwertige Zelldaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   data    Vektor mit den zu schreibenden Daten
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten

    template<typename type>
    void addCellDataToTimestep(const FullVector<type>& data, const size_t timestep, const std::string name);


///\brief   Fuegt skalarwertige Punktdaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   data[]  Array der Quelldaten
///@param   length  Länge des Arrays
///@param   name    Name der Daten

    template<typename type>
    void addPointDataToAll(const type data[], size_t length, std::string name);



///\brief   Fuegt skalarwertige Punktdaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   data    Vektor mit den zu schreibenden Daten
///@param   name    Name der Daten

   template<typename type>
   void addPointDataToAll(const FullVector<type>& data,const std::string name);


///\brief   Fuegt skalarwertige Zelldaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   data[]  Array der Quelldaten
///@param   length  Länge des Arrays
///@param   name    Name der Daten

    template<typename type>
    void addCellDataToAll(const type data[], size_t length, std::string name);


///\brief   Fuegt skalarwertige Zelldaten zu allen Zeitschritten hinzu
///         Die Daten werden in der Datei als float abgespeichert
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   data    Vektor mit den zu schreibenden Daten
///@param   name    Name der Daten

    template<typename type>
    void addCellDataToAll(const FullVector<type>& data, const std::string name);

///\brief   Fuegt vektorwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   datax[] Array der x Komponenten
///@param   datay[] Array der y Komponenten
///@param   dataz[] Array der z Komponenten
///@param   length  Länge der Arrays
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten

    template<typename type>
    void addPointVecToTimestep(const type datax[], const type datay[], const type dataz[], const size_t length, const size_t timestep, std::string name);


///\brief   Fuegt vektorwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type        Datentyp der Quelldaten
///@param   datax       Vektor mit den x Komponenten der zu schreibenden  Daten
///@param   datay       Vektor mit den y Komponenten der zu schreibenden  Daten
///@param   dataz       Vektor mit den z Komponenten der zu schreibenden  Daten
///@param   timestep    Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name        Name der Daten

    template<typename type>
    void addPointVecToTimestep(const FullVector<type>& datax, const FullVector<type>& datay, const FullVector<type>& dataz, const size_t timestep, const std::string name);



///\brief   Fuegt vektorwertige Zelldaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   datax[] Array der x Komponenten
///@param   datay[] Array der y Komponenten
///@param   dataz[] Array der z Komponenten
///@param   length  Länge der Arrays
///@param   timestep Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name    Name der Daten

    template<typename type>
    void addCellVecToTimestep(const type datax[], const type datay[], const type dataz[], const size_t length, const size_t timestep, std::string name);

///\brief   Fuegt vektorwertige Punktdaten zu einem Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type        Datentyp der Quelldaten
///@param   datax       Vektor mit den x Komponenten der zu schreibenden  Daten
///@param   datay       Vektor mit den y Komponenten der zu schreibenden  Daten
///@param   dataz       Vektor mit den z Komponenten der zu schreibenden  Daten
///@param   timestep    Zeitschritt, dem die Daten hinzugefuegt werden
///@param   name        Name der Daten

    template<typename type>
    void addCellVecToTimestep(const FullVector<type>& datax, const FullVector<type>& datay, const FullVector<type>& dataz, const size_t timestep, const std::string name);


///\brief   Fuegt vektorwertige Punktdaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   datax[] Array der x Komponenten
///@param   datay[] Array der y Komponenten
///@param   dataz[] Array der z Komponenten
///@param   length  Länge der Arrays
///@param   name    Name der Daten

    template<typename type>
    void addPointVecToAll(const type datax[], const type datay[], const type dataz[], size_t length, std::string name);

///\brief   Fuegt vektorwertige Punktdaten zu allen Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type        Datentyp der Quelldaten
///@param   datax       Vektor mit den x Komponenten der zu schreibenden  Daten
///@param   datay       Vektor mit den y Komponenten der zu schreibenden  Daten
///@param   dataz       Vektor mit den z Komponenten der zu schreibenden  Daten
///@param   name        Name der Daten

    template<typename type>
    void addPointVecToAll(const FullVector<type>& datax, const FullVector<type>& datay, const FullVector<type>& dataz, const std::string name);



///\brief   Fuegt vektorwertige Zelldaten zu allen Zeitschritten hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type    Datentyp der Quelldaten
///@param   datax[] Array der x Komponenten
///@param   datay[] Array der y Komponenten
///@param   dataz[] Array der z Komponenten
///@param   length  Länge der Arrays
///@param   name    Name der Daten

    template<typename type>
    void addCellVecToAll(const type datax[], const type datay[], const type dataz[], size_t length, std::string name);


///\brief   Fuegt vektorwertige Zelldaten zu allen Zeitschritt hinzu.
///         Die Daten werden in der Datei als float abgespeichert.
// ------------------------------------------------------------------------------------------------
///@tparam  type        Datentyp der Quelldaten
///@param   datax       Vektor mit den x Komponenten der zu schreibenden  Daten
///@param   datay       Vektor mit den y Komponenten der zu schreibenden  Daten
///@param   dataz       Vektor mit den z Komponenten der zu schreibenden  Daten
///@param   name        Name der Daten

    template<typename type>
    void addCellVecToAll(const FullVector<type>& datax, const FullVector<type>& datay, const FullVector<type>& dataz, const std::string name);

};

}//namespace Icarus

#include "vtkwriter.tpp"

#endif // __VTKWRITER_HPP_
