Aufbau:
"
127.0.0.1 slots=1
icarus-xx slots=y
...
"
Gerade das was zwischen den Gänsefüßen steht, die Punkte sind aber nur Platzhalter
und sollen zeigen, dass noch andere icarus Nodes folgen können!

erste Zeile: lokaler Host (Gateway Knoten)(nicht notwendig)
zweite bis beliebige Zeile: Nodes als Host
slots enstrpicht der Anzahl der Prozesse pro Host

Aufruf:
mpirun -np X -hostfile path/to/hostfile.txt icarus


