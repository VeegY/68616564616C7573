// das hier ist mehr eine Demo, als ein Test

#include"include/discretizer.hpp"   // Klassendeklarationen
#include"include/discretizer.cpp"   // Funktionen der Klassen und Funktion discretizer()

#include<iostream>  // zu testzwecken

using namespace Icarus;

int main()
{
    discretizer("threeblocks.obj", "threeblocks.list", 0.1f, 100, 100, 100);
    return 0;
}
