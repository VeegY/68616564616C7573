// das hier ist mehr eine Demo, als ein Test

#include"include/discretizer.hpp"   // Klassendeklarationen
#include"include/discretizer.cpp"   // Funktionen der Klassen und Funktion discretizer()

using namespace Icarus;

int main()
{
    discretizer("board.obj", "board.list", 0.0005f, 200, 200, 25);
    return 0;
}
