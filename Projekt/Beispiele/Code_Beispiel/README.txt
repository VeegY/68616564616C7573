===============================================================================
===				  CMake Readme				   ====
===============================================================================

In diesem kleinen Tutorial moechte Ich kurz erklaeren wie hier die das Built-
System aufzusetzen ist. Dabei werde Ich darauf eingehen, wie man mithilfe CMake
das Build-System aufsetzt, das Build anschliessend kompiliert und schliesslich
wie man die einzelnen Programme ausfuehren kann. Abschliessend werde Ich noch
darauf eingehen wie sich die Doxygen Dokumentation kompilieren laesst.

-------------------------------------------------------------------------------
			   Aufsetzen des Build-Systems
-------------------------------------------------------------------------------
Zunaechst sollte man einen Blick in das 'build'-Verzeichnis werfen. Dieses
Verzeichnis ist zunaechst leer. Schaut man nun in das 'Code'-Verzeichnis
findet man dort den Code den wir gerne kompilieren wuerden. Dafuer muessen wir
aber erstmal das Build aufsetzen. Anleitung:

1. Oeffne die Konsole und navigiere ins 'build'-Verzeichnis
2. Lade das CMake-Modul (module load cmake)
3. Lade das Build System (cmake ../)

Wenn wir uns nun wieder das 'build'-Verzeichnis anschauen sieht man, dass hier
einige Verwaltungsdateien erzeugt wurden. Diese sind 'Out of Source', da wir
nicht wollen, dass unser 'Code'-Verzeichnis unnoetigigerweise zu
unuebersichtlich wird. Des weiteren wurde ein 'EXE' Ordner erstellt
der zunaechst leer ist. Bitte beachtet dass sich die 'CMakeLists.txt'-Hauptdatei 
fuer unser Projekt im 'Code'-Verzeichnis befindet. (also cmake ../Code)

-------------------------------------------------------------------------------
			     Kompillieren des Codes
-------------------------------------------------------------------------------
Nun wollen wir den Code im 'Code'-Verzeichnis kompilieren. Dazu muessen wir den
Code ueber das Buildsystem mit den Compilern verknuepfen und anschliessend
ausfuehren. Anleitung:

1. Kompilieren des Codes (make all)
2. Fertig

Das war es schon! Mit nur einer Zeile wurde der komplette Code kompiliert und
wenn man nun in das 'EXE'-Verzeichnis schaut sieht man, dass hier der Code die 
einzelnen Programme gespeichert wurden. Um diese nun auszufuehren navigiert man
in der Konsole ins eben dieses Verzeichnis und fueht einfach wie gewohnt das 
gewuenschte Programm aus. (EXE/[Programm] in der Kommandozeile)

-------------------------------------------------------------------------------
			   Hinzufuegen von eigenem Code
-------------------------------------------------------------------------------
Was ist aber nun wenn ihr euren eigenen Code dem Build System hinzufuegen und
kompilieren wollt? Dazu habe Ich einen Ordner 'Neuer_Code' im Hauptverzeichnis
erstellt, der eine CSR-Implementierung enthaelt. Um diesen Code dem Build-
System hinzuzufuegen bedarf es lediglich weniger Zeilen Code. Anleitung:

1. Verschiebe den CSR Ordner in das 'Code'-Verzeichnis
2. Oeffne die 'CMakeLists.txt'-Datei im 'Code'-Verzeichnis
3. Ergaenze das Build-Systen um das 'CSR'-Verzeichnis
	include_directories(${CMAKE_CURRENT_SOURCE_DIR} CSR)
4. Teile CMake mit, dass sich im 'CSR'-Verzeichnis eine weitere 'CMakeLists.txt'-Datei befindet
	add_subdirectory(CSR)
5. Erzeuge eine 'CMakeLists.txt'-Datei im 'CSR'-Verzeichnis
6. Uebergebe CMake die zu kompilierende Quelltextdatei
	add_executable(CSR_CPU_Test CSR_CPU_TEST.cpp)
7. Starte CMake erneut (cmake ../)
8. Kompiliere den Code (make all)

Wie wir nun in der Konsole erkennen koennen wurde dem Code das 'CSR_CPU_TEST'-
Programm hinzugefuegt. Dieses laesst sich nun wie gewohnt im 'EXE'- Verzeichnis 
ausfuehren. (EXE/[Programm] in der Kommandozeile)

-------------------------------------------------------------------------------
			 Kompilieren der Doxygen Dokumentation
-------------------------------------------------------------------------------
Zunaechst werfen wir einen Blick in das 'Dokumentation'-Verzeichnis. Darin
befindet sich die Konfigurationsdatei von Doxygen 'doxygen_config'. Ausserdem
habe Ich im Projektverzeichnis zusaetzlich das Doxygen Manual beigelegt. Das
kompilieren der Dokumentation ist sehr leicht. Anleitung:

1. Kompilieren der Doxygen Dokumentation (make doxy)
2. Fertig

Wir sehen nun, dass viele Dateien im 'Dokumentation'-Verzeichnis dazugekommen
sind. Wenn wir nun eine html Datei oeffnen, wird im Browser die Dokumentation
erzeugt.
