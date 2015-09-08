===============================================================================
===				Doxygen Readme				   ====
===============================================================================

In diesem README moechte Ich euch nahelegen wie Doxygen dazu verwendet werden
kann, um die Struktur unseres Codes in einer übersichtlichen Dokumentation zu
veranschlaulichen. Doxygen verwendet eine einheitliche Syntax die zur 
Kommentierung und Beschreibung des Codes angewendet wird. In diesem Ordner
habe Ich einige Doxygen-Templates angefuegt, die ihr in eueren Code kopieren
koennt. Viel Spass.

Im ersten Schritt dieses README's werde Ich darauf eingehen wie Ihr die 
Dokumentationkompilieren koennt. Ausserdem werde Ich im zweiten Teil naeher
auf die einzelnen Kommentierungen eingehen.

===============================================================================
===		   Kompilierung und Nutzung der Dokumentation		   ====
===============================================================================

Die Doxygen Kommentierung laesst sich mithile des Build System's ueber den
Compiler erzeugen. Bevor wir damit starten koennen, benoetigt Doxygen eine 
Setup Datei, die als 'doxygen_config' im Dokumentations-Verzeichnis gespeichert 
ist.

Ich moechte nun auf die einzelnen Konfigurationsoptionen eingehen die Ich 
(bisher) fuer unsere Dokumentation gewaehlt habe. Ausfuehrliche Erlauertungen
zu den findet ihr in der Konfigurationsdatei. Der gesamte Dokumenationsordner 
wird fortan nicht (regelmaessig) im Versionskontrollsystem commitet. Ihr koennt
also gerne eigene Konfigurationen erstellen.

-------------------------------------------------------------------------------

PROJECT_NAME           = "Studienprojekt 2015_16"
Mit dieser Variable laesst sich der Titel der Dokumentation bearbeiten.

OUTPUT_LANGUAGE        = German
Diese Variable steuert die in der Dokumentation verwendete Sprache.

- - - WIP - - - : Muss noch getestet werden
OPTIMIZE_OUTPUT_FOR_C  = YES
Optimiert die Ausgabe auf die Verwendendung von C-Code.
- - - WIP - - - : Muss noch getestet werden

QUIET                  = YES
WARNINGS               = YES
Mit dieser Konfiguration beschraenkt sich die Compilerausgabe im Kompilierungs-
prozess nur noch auf Warnungen (hervorgerufen durch falsche Kommentierung).

INPUT                  = ../Code/
Mithilfe dieser Variable wird das Code-Verzeichnis gewaehlt.

RECURSIVE              = YES
FILE_PATTERNS          = *.c *.cpp *.h *.hpp *.cu
Aktivieren der rekursiven Suche des Code-Verzeichnisses nach den mithilfe
FILE_PATTERNS ausgewaehlten Dateitypen.

GENERATE_HTML          = YES
GENERATE_LATEX         = NO
GENERATE_RTF           = NO
GENERATE_XML           = NO
Gewaehlte Formate in denen die Dokumentation erzeugt werden soll. Hier: HTML.

HTML_OUTPUT            = ../Dokumentation
Die Dokumentation wird in dem hier angegbenem Verzeichnis erzeugt.

- - - WIP - - - : Muss noch getestet werden
GENERATE_TREEVIEW      = YES
- - - WIP - - - : Muss noch getestet werden

-------------------------------------------------------------------------------

Mithilfe 'make doxy' im 'Build'-Verzeichnis laesst sich die Dokumentation
erzuegen. Des weiteren wurde auch wurde die Doxygen Dokumentation in 'ALL' 
aufgenommen. Das bedeutet mit 'make all' laesst sich sowohl der Code, als auch
die Doxygen Dokumentation kompilieren.

Um die Dokumentation zu oeffnen reicht es aus die 'annotated.html' -Datei im
Browser auszufuehren. Des Weiteren befindet sich das offizielle Doxygen-Manual
im Dokumentationsverzeichnis.

===============================================================================
===			Erlaeuterung der Kommentierung		  	   ====
===============================================================================

Ein Doxygen-Kommentar laesst sich in C-Code mithilfe '///' einleiten und der
makierte Code wird von Doxygen ausgewertet. Die wichtigsten Befehle moechte
Ich nun naeher erlauern. Eingehendere Beispiele finden sich auch im Manual.

-------------------------------------------------------------------------------

\brief
Der 'brief' nachfolgende Text generiert eine Kurzbeschreibung des kommentierten
Codes. Die Kommentierung wird mit einem Punkt abgeschlossen und jeglich
weiterer Text wird als detaillierte Beschreibung des Codes verwendet.
Dieser Befehl wird jeweils dem zu beschreibenden Code vorangestellt.
Beispiel:
///\brief
///Diese Klasse implementiert einen Kreis.
///Diese Klasse implementiert einen Kreis mit dem Radius 'r'.
class Kreis
{
double _r;

Kreis(double r) : _r(r) {}
}

\param
Mit 'param' kann die Beschreibung eines Parametes einer Funktionen erzeugt
werden. Fuer Templates ist der entsprechende Befehl 'tparam'.
Dieser Befehl wird jeweils dem zu beschreibenden Code vorangestellt.
Beispiel:
///\param r Radius des Kreises
Kreis(double r) : _r(r) {}

<
Mit diesem Vergleichszeichen lassen sich einzelne Variablen beschreiben.
Dieser Befehl wird jeweils der zu beschreibenden Variable nachgestellt.
Beispie:
double d(2*r); ///< d dim
