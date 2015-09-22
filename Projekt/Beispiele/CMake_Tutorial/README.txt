===============================================================================
===				  CMake Readme				   ====
===============================================================================

Mithilfe CMake laesst sich das Build-System eines Projekts verwalten. Dazu wird
zunaechst ein sogenanntes build-Verzeichnis benoetigt. In unserem Fall ist dies
'build'. Hier sind alle Verwaltungsdateien des Build-Systems untergebracht. Ich
habe CMake gerade so eingerichtet, dass die Verwaltungsdateien von CMake und 
Make in einer Verzeichnisstruktur gespeichert wird, die innerhalb von 'build'
erzeugt wird. Insbesondere weist sie dieselbe Verzeichnisstruktur wie das
'Code'-Verezichnis auf.

Zunaechst moechte Ich kurz auf die Nutzung der verschiedenen Makros eingehen,
die Ich implementiert habe. Anschliessend werde Ich die wesentlichen Elemente
der 'CMakeLists.txt'-Dateien erlaeutern, auf denen CMake basiert.

===============================================================================
===				     Nutzung				   ====
===============================================================================

build/cmake ../Code
Im 'build'-Ordner wird CMake mit einem Verweis auf das 'Code' ausgefuehrt.
Damit sucht CMake nach der 'CMakeLists.txt' in 'Code'. Zunaechst wird die 
Verzeichnisstruktur erstellt, indem die gelinkten Subverzeichnisse nach eigenen
CMakeLists.txt durchsucht werden.
Anschliessend erzeugt der Compiler die in der 'CMakeLists.txt' definierten
Makros. Die entsprechenden Verwaltungsarchitekturen werden im 'build'-
Verzeichnis erstellt.

Buildoption CMAKE_BUILD_TYPE:
CMake laesst in verschiedenen Buildvarianten kompilieren. Dabei unterscheiden
wir zwischen den Build Typen 'None', 'Debug' und 'Release'. Anleitung:
   build/cmake -DCMAKE_BUILD_TYPE=[Build_Typ] ../Code

'None':
Standardeinstellung die ohne Spezifikation (wie oben) verwendet wird
'Debug':
Die Compilerflags koennen/werden gerade so gewaehlt, dass das der Code
hinsichtlich des Debugvorgangs optimiert wird.
'Release':
Das Build System wird dahingehend ausgelegt, dass der Compiler den Code
hinsichtlich der Leistungsfaehigkeit optimiert.

(Optional lassen sich auch weitere und eigene Build Typen aufsetzen, fuer
uns reicht das so jedoch völlig aus)

MAKRO: make all / make
Mithilfe dieses Makros wird das komplette build (bzw. alle dem Makro ALL 
zugewiesenen Befehle) kompiliert. Das heisst, dass alle Compiler aufgesetzt
werden und die Verzeichnisstruktur nach Programm-Makros durchsucht wird.
Dabei werden die einzelnen Programme mit den angegebenen Compilerflags
kompiliert. Dabei habe Ich CMake so aufgesetzt, dass im 'Build'-Verzeichnis ein
'EXE'-Verzeichnis erzeugt wird, in das automatisch saemtliche
ausfuehrbaren Programme gespeichert werden.

MAKRO: make doxy
Dieser Befehl erzeugt die Doxygen-Dokumentation im Verzeichnis 'Dokumentation'.

- - - WIP - - - : Muss mich damit noch eingehender auseinandersetzen.
ctest
Mithilfe CMake lassen sich ausserdem die einzelnen Unit-Tests ausfuehren.
Einige Beispiele dazu sind innerhalb des Verzeichnisses 
'Beispiele/Doxygen_Beispiel/Code/'
zu finden.
- - - WIP - - -

===============================================================================
===				  CMakeLists.txt			   ====
===============================================================================

Die 'CMakeLists.txt'-Dateien sind die primaeren CMake Veraltungsdateien. In
jedem Verzeichnis, das dem Build-System hinzugefuegt werden soll muss sich eine
'CMakeLists.txt'-Datei befinden. Hier werden dann entsprechend der unten
erlaeuterten Befehle Compiler aufgesetzt, Verzeichnisse verknuepft und
schliesslich die Progamme kompiliert.

-------------------------------------------------------------------------------

find_package(Beispiel_Paket)
Verschiedene Pakete lassen sich mithilfe des 'find_package'-Befehls laden.

set(CMAKE_CXX_FLAGS "-[Compilerflag]")
Hier werden die Compilerflags des GCC Compilers definiert. Gleiches gilt fuer
CUDA_NVCC_FLAGS
Um den Compilerflags weitere Flags hinzuzufuegen wird schlicht
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -[Compilerflag]")
verwendet.

include_directories(${PROJECT_SOURCE_DIR} [Ordner1] [Ordner2])
Mit diesem Befehl laessen sich beliebig viele Verzeichnis dem Build-System
hinzufuegen. Es ist empfehlenswert jeweils nur die Ordner zu linken, die sich
in dem selben Verzeichnis der jeweiligen 'CMakeLists.txt'-Datei befinden.

add_subdirectory([Ordner])
Mit diesem Befehl wird dem Build-System mitgeteilt, dass sich in dem angegebenem
Subverzeichnis eine 'CMakeLists.txt'-Datei befindet. Dieses Subverzeichnis muss
vorher mit dem 'include_directories'-Befehl eingefuehrt worden sein.

add_executable([Makroname] [Dateiname].cpp)
Dieser Befehl kompiliert die angegebene Datei mit den angebenen Compilerflags
ueber den GCC-Compiler. Dabei wird im 'build/EXE'-Verzeichnis das zugehoerige 
Programm erzeugt.

add_test([name] [command])
Mithilfe dieses Befehls laesst sich fuer ctest ein Test hinzufuegen. Damit 
lassen sich Unit-Tests umsetzen. Anschliessend wird der Test mit 'make test' 
oder 'ctest (-V)' ausgefuehrt. Die Tests sollten keine Eingaben erfordern.

cuda_add_executable([Makroname] [Dateiname].cu)
Das entsprechende CUDA-Aequivalent zu 'add_executable'.

message(${CUDA_NVCC_FLAGS})
Mit diesem Befehl lassen sich ausgaben innerhalb des build-Vorgangs einbinden.
In diesem Fall werden die verwendeten Compiler-Befehle ausgegeben.

