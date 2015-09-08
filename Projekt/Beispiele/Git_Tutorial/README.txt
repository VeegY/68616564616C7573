===============================================================================
===				  GIT Readme				   							   ====
===============================================================================

Git ist eine Versionenkontrolle und hilft beim parallelen Arbeiten am 
Code für unser Projekt. Git bietet zusätzlich die Möglichkeit Zweige(branch) zu 
bilden und diese später wieder zu vereinen (mergen). Auf einem (zum diesem 
Zeitpunkt  noch unbekannten) Server gibt es ein Repository, dass die Grundlage 
bildet. Dieses Repository kann man "clonen" und dann lokal am eigenen Rechner,
auch offline, programmieren. Alle Änderungen müssen manuell "commited" werden. 
Damit Git weiß welche Dateien "commited" werden sollen, muss man sie mit "add" 
in eine Zwischenzustand bringen. Durch dieses "commit" erzeugt man lokal eine neue 
Version. Soll diese Version für alle zugänglich sein, kann man seine Version 
auf den Server "pushen". Unten stehen viele Befehle nochmal detaillierter 
beschrieben.



===============================================================================
===				     Nutzung				   							   ====
===============================================================================

Speziell für unser Projekt haben wir uns folgende Branch-Struktur überlegt:

Master-Branch:
Dieser Branch stellt das entgültige Produkt dar. Ihr sollt und könnt nicht 
auf diesen Branch pushen. Alles was im Master-Branch landet ist laufender und
getesteter Code.
Ihr sollt aber von diesem Master-Branch alle Neuigkeiten erhalten die ihr braucht,
könnt also einen "pull" auf dem Master-Branch durchführen.

Develop-Branch:
In diesem Bereich landen eure Code-Teile. Arbeitet ihr lokal an einer Funktion 
dürft ihr diese nach "Develop" pushen sobald sie lokal auf eurem Rechner läuft.

Debug-Branch:
Dieser Branch wird von euch erstellt und ist für alle zugänglich sobald ihr ein
Problem oder Bug in eurem Code habt. Ein Kürzel im Branch-Namen wäre beim erstellen
wünschenswert um elf gleichnamige Debug-Branches zu vermeiden (zb. Debug-Kevin).
Zusätzlich gibt es einen Ordner im Repository für .txt-Dateien auf denen ihr
euren Bug für andere näher beschreiben könnt.

===============================================================================
===				     git-commands										   ====
===============================================================================

****	globale Infos				****
git config --global user.name "Kevin"							
git config --global user.email "kevin.hollmann@tu-dortmund.de"	


****	Arbeiten mit dem Server		****
git clone username@host:/path/to/repository (Username und Server stehen noch nicht fest)
git remote add origin <server> 				(legt das push-Ziel fest)
git remote -v 								(listet alle eingestellten remote-repositorys auf)
git push origin <branchname>  				(für euch: branchname=develop)
git pull 									(merged Serveraktualisierungen mit euren)


**** 	Lokales Arbeiten 			****
git add <filename> 							(Alle Dateien adden die ihr commiten wollt)
git commit -m "Message"						(Alle "add"-Datein werden lokal commited, mit Notiz bitte)


**** 	branches					****
git checkout -b <branchname> 				(Erstellt einen neuen Branch namens "branchname")
git checkout <branchname>					(Wechselt auf den branch "branchname")
git branch									(Listet alle branches auf)