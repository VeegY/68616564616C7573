===============================================================================
===		 	   Ueber Git-Versionskontrolle		   	   ====
===============================================================================

In Softwareprojekten, in denen verschiedene Leute an demselben Code arbeiten
wollen wird mit Git ein Versionskontrollsystem bereitgestellt. Dafuer kann der 
temporaere Status der Dateien eines Verzeichnisses als 'commit' unter Verwendung 
einer speziellen 'Branch-Struktur' gespeichert werden. Dazu spaeter mehr.

Außerdem laesst sich damit die Arbeit verschiedener Nutzer am Code fuer das 
Projekt synchronisieren. Im Rahmen unseres Studienprojekts haben wir dafuer 
ueber die Github Platform ein sogenanntes Repository bereitsgestellt, mithilfe 
dessen sich der Arbeitsaustausch verwalten laesst. Alle weiteren Informationen 
zu Github, werden in einem eigenen (kurzen) Tutorial zusammengefasst.

Zunaechst moechte Ich die Nutzung des Git-Repos erklaeren. Dabei werde Ich
ausserdem naeher erlaeutern wie Git im Rahmen unseres Projekts genutzt werden 
sollte um die Arbeit fuer alle Beteiligten moeglichst angenehm zu gestalten.
Abschließend werde Ich noch auf die Struktur des Github-Repositorys eingehen.

===============================================================================
===				   Nutzung				   ====
===============================================================================

Zunaechst sei angemerkt, dass sich fuer alle git-Befehle eine Dokumentation 
aufrufen laesst. Mit 'Enter' lassen sich die weiteren Zeilen der Hilfe scrollen.
Um die Hilfe wieder zu verlassen reicht 'q' zum 'quitten'.
   git --help [Befehl] -> Bsp: git --help clone
Eine weitere gute Referenz ist
   https://git-scm.com/book/de/v1

------------------------Herunterladen des Repositories-------------------------
git 
Bevor mit der Arbeit am Code begonnen wird, sollte sich jeder User eine lokale
Kopie des Github-Repositorys auf seinem Mathe-Account erstellen. 
   git clone https://github.com/VeegY/stuprotm1516.git
Der Server wird nun mithilfe des Befehl
   git remote -v
als sogenanntes Remote Repository 'origin' der lokalen Kopie aufgefuehrt.
Alternativ laesst sich dem Server auch manuell ein Remote Name zuordnen.
   git remote add [remote-name] [server]
    -> Bsp: git remote add stuprotm https://github.com/VeegY/stuprotm1516.git
Diese Variante laedt das Repository jedoch nicht herunter.

Auf dieser lokalen Kopie koennt ihr nun euren Code erstellen und mit dem Code 
der Anderen kombinieren und testen.

---------------------------Konfiguration von Github----------------------------

Damit deine Beitraege im Github-Repository zugeordnet werden koennen, muessen 
zunaechst die wichtigsten Daten in deiner Git-Konfiguration registriert werden. 
   git config --global user.name "[username]"							
   git config --global user.email "[email]"
Mithilfe '--global' wird git dabei gerade so auf aufgesetzt, dass alle weiteren 
git-Projekte mit denselben Angaben konfiguriert werden.
Desweiteren laesst sich so auch der Standardeditor einstellen, der fuer Commit-
Nachrichten verwendet werden soll (etwa gedit oder emacs). Dazu spaeter mehr.
   git config --global core.editor [Editor]
Um deine Angaben zu ueberpruefen verwende
   git config --list

Des Weiteren ist im 'stuprotm1516'-Verzeichnis die Konfiguration '.gitignore' 
gespeichert. Hier koennen Dateitags und Verzeichnisse angegeben werden, die von 
der Versionskontrolle ausgeschlossen werden sollen. Diese Datei wird als 
Verwaltungsdatei nicht im Verzeichnis angezeigt, kann aber ueber 
   gedit .gitignore
bearbeitet werden.

---------------------------Verwendung von Branches-----------------------------

Git arbeitet mit einer sogenannten Branch-Struktur. Dabei handelt es sich um 
spezielle Zeiger, die auf den Commit-Zustand des Verzeichnisses referenzieren. 
Das heisst, dass verschiedene Versionen des Repositories (unter gewissen 
Einschraenkungen) gleichzeitig existieren koennen. In unserem Fall waere die 
erste Anlaufstelle der Branch 'develop'. Die Liste aller lokalen Branches 
laesst sich mit folgendem Befehl einsehen
   git branch 
Mit 'git branch -r' lassen sich alle Branches des Github-Repositorys anzeigen.
Anschliessend waehlen wir den Branch aus auf dem wir arbeiten wollen
   git checkout [local-branch]
Wir sind damit auf einen anderen Branch gewechselt. Dabei wird der Zustand des 
lokalen Repositorys mit dem Zustand des Branches synchronisiert auf den wir 
gewechselt sind. Der aktive Branch ist stets mit dem Zeiger 'HEAD' referenziert. 
Alle Dateien, werden auf den Zustand gesetzt zu dem der Branch commitet wurde.

Mithilfe des 'branch'-Befehls laesst sich ein neuer Branch erzeugen. Dabei wird 
der aktuelle Zustand des lokalen Repositorys als Zeiger des Branches verwendet
   git branch [local-branch]
Eine weitere Moeglichkeit einen Branch zu erstellen und insbesondere direkt 
darauf zu wechseln ist
   git checkout -b [local-branch]
Moechte man einen Branch loeschen, laesst sich dies mithilfe
   git branch -d [local-branch]
erzielen.

--------------------------Versionskontrolle mit Git----------------------------

Die Dateien in deinem (lokalen) Repository basieren auf einem Statussystem mit 
verschiedenen Dateizustaenden. Diese Dateizustaende ergeben in der Summe den 
'commit'-Zustand des Verzeichnisses. Der aktuelle Zustand der Dateien, bzw des 
Repos laesst sich wie folgt pruefen
   git status

Eine Datei die sich im Grundzustand wird als 'unmodified' bezeichnet. Dieser 
Zustand wird erzielt, falls eine Datei seit des letzten commits nicht veraendert 
wurde. Eine bearbeitete Datei nimmt entsprechend den Zustand 'modified' an. Ist
eine Datei innerhalb des aktuellen branches nicht bekannt wird sie als 
'untracked' markiert.

Ein Zentrales Element der Versionskontrolle mit Git stellt die sogenannte 
'Staging Area' dar. Eine Datei kann dieser mit den Befehlen 'add' und 'reset'
hinzugefuegt, bzw. entfernt werden.
   git add [filename]
   git reset [filename]
Dabei sind die Dateien genau in dem Zustand erfasst, wie sie ge'add'ed wurden. 
Das heisst wird eine Datei spaeter geaendert, muss sie erst wieder geadded 
werden, damit die Aenderungen erfasst werden. Mehrere Dateien koennen der 
Staging Area unter Verwendung der Stern-Schreibweise
   git add [Verzeichnis]/\*.[tag] -> Beispiel: git add Projekt/Code/\*.cpp
simultan und strukturiert hinzugefuegt werden. (resp. reset)

Der Finale Schritt um die Aenderungen an den Dateien geltend zu machen ist es
die Staging Area zu commiten. Dabei wird der Zeiger des aktiven Branches auf ein
Abbild des Repositories verschoben, das die Aenderungen der Dateien erfasst.
   git commit -m "[message]"		
Dabei sei darauf hingewiesen, dass ein commit stets mit einer Commit-Nachricht 
versehen werden sollte. Andernfalls kommt es nicht zu einem Commit. Alternativ 
wird die Commit-Nachricht ohne '-m' mithilfe des gewaehlten Editors erzeugt 
Um die Aenderungen eines Commits einzusehen wird
   git log
verwendet. Einzelne Dateien lassen sich aus dem Commit-Zustand des Repositories
mithilfe des Befehls 'rm' entfernen. Dabei sei darauf hingewiesen, dass diese
Datei nicht im Verzeichnis sondern nur aus der Versionskontrolle geloescht wird.
   git rm --cached [filename]
Ohne '--cached' wird die Datei auch innerhalb des Verzeichnis geloescht.

-----------------Synchronisierung mit dem Github Repository--------------------

Bevor wir aber beginnen Dateien in unserem Github-Repository herumzuschieben
und damit alles durcheinander bringen, sollten wir zunaechst die aktuelle Kopie 
des Github Repositories herunterladen. Eine Moeglichkeit ist
   git pull [remote-name]
Dieser Befehl laedt alle Daten herunter, die noch nicht in deinem lokalen
Repository sind und updatet gleichzeitig alle Dateien des Repositorys. Dazu ist
es jedoch erforderlich, dass der gewuenschte Branch des Remote Repository als
Ziel von git clone verwendet wurde.

Eine alternative und gezieltere  Vorgehensweise kann mit einer Kombination von
fetch und merge erzielt werden. Dazu verwendet man
   git fetch [remote-name] -> Bsp: git fetch origin
Dieser Befehl laedt die Daten des Github Branches [remote-name] als lokale Kopie.
Um einen Remote Branch mit dem aktuell aktiven lokalen Branch zu verbinden
   git merge [remote-branch] -> Bsp: git merge origin/develop
Falls ihr es entdeckt, lasst erstmal besser die Finger von 'rebase'!

Andersherum laesst sich der aktive Branch des lokalen Repository mit dem
gewaehlten Branch des Remote Repositorys mergen.
   git push [remote-name] [remote-branch] -> Bsp: git push origin origin/develop

------------------------------Github Repository--------------------------------

Github:
Wenn ihr einen Branch im Github-Repository loeschen wollt, waehlt im
Projektverzeichnis 'Branches' und loescht ueber das Muelleimer-Icon. Loeschen
koennen nur der Ersteller des Branches sowie der Eigentuemer des Repositories.
Das ist in diesem Fall Kevin.

===============================================================================
==				Projekt-Struktur			   ====
===============================================================================

Speziell fuer unser Projekt haben wir uns folgende Branch-Struktur ueberlegt:

Master-Branch:
Dieser Branch stellt das entgueltige Produkt dar. Ihr solltet (und koennt) nicht 
auf diesen Branch pushen. Die Idee ist es hier, dass der Code der auf den
Master-Branch gepushed wird fuer alle Architekturen compilierbar und ausfuehrbar
sein soll. Damit wird erzielt, dass stets eine Version existiert, die fuer alle
verwendbar ist.

Develop-Branch:
In diesem Bereich landen eure Code-Teile. Arbeitet ihr lokal an einer Funktion 
duerft ihr diese nach "Develop" pushen sobald sie lokal auf eurem Rechner laeuft.

Debug-Branches:
Diese Branches sollten von euch erstellt werden, wenn ihr ein Problem oder Bug
in eurem Code habt und fuer die Loesung Hilfe benoetigt. 
Ein Kuerzel im Branch-Namen waere beim erstellen wuenschenswert um gleichnamige 
Debug-Branches zu vermeiden (zb. Debug-Kevin).
Zusaetzlich gibt es einen Ordner im Repository fuer Bug-Reports (.txt-Dateien)
auf denen ihr euren Bug naeher erlaeutern koennt.

