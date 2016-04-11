Dieser Ordner enthält eine (stand-alone) Demo, die die Skalierbarkeit
einer dünnbesetzten Bandmatrix-Vektormultiplikation sowie eines
darauf aufbauenden CG-Lösers (z.B. 5-Punkte-Stern) auf Icarus demonstriert.

Aufruf z.B. mit:

mpirun --hostfile my_hosts -np 8 --map-by node ./mvb 4 d 3000

D.S. 9.4.16
