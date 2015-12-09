/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                proto.hpp
* Erstellt:                 23.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Prototypen für häufige Klassenarten.
*/

#ifndef __PROTO_HPP_
#define __PROTO_HPP_

namespace Icarus
{

class NonCopyable
{
protected:
    NonCopyable() { }
private:
    NonCopyable(const NonCopyable&) = delete;
    NonCopyable& operator=(const NonCopyable&) = delete;
};

class Interface
{
protected:
    Interface();
};

}

#endif // __PROTO_HPP_
