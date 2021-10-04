#! /usr/bin/python3

import sys
import csv


# Read-Klasse: Speichert die Sequenz und den Namen eines Reads
class Read:
    # Initialisierung aus zwei FASTA-Zeilen: Header in der ersten, Basen in der zweiten
    def __init__(self, lines):
        self.name = lines[0].strip()
        self.bases = lines[1].strip()

    # Fügt die k-mere der Länge kmersize des Reads zu der übergebenen kmer-Liste hinzu
    def add_kmers(self, kmersize, kmers):
        for pos in range(0, len(self.bases)-kmersize+1):
            kmer = self.bases[pos:(pos+kmersize)]
            if kmer not in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1
    
    # Ausgabemethoden zum Debuggen
    def __str__(self):
        return self.name + ": " + self.bases[:20] + "..."
    
    def __repr__(self):
        return self.__str__()


# Klasse zum Speichern eines Knotens im DBG
class DBGnode:
    # Sequenz speichern und zwei Dictionaries für die ausgehenden (eto)
    # und eingehenden (efrom) Kanten. Diese Dictionaries halten als Schlüssel
    # jeweils DBGnodes, zu/von denen die Kante geht und als Wert die Information,
    # wie häufig diese Kante hinzugefügt wurde.
    def __init__(self, seq):
        self.seq = seq
        self.eto = dict()
        self.efrom = dict()
    
    # Methoden, um eine Kante von/zu dem Knoten hinzuzufügen
    def add_edge_to(self, eto):
        if eto not in self.eto:
            self.eto[eto] = 0
        self.eto[eto] += 1
    
    def add_edge_from(self, efrom):
        if efrom not in self.efrom:
            self.efrom[efrom] = 0
        self.efrom[efrom] += 1
   
    # Methoden, um die Sequenzen der anderen Knoten zu generieren, zu/von denen
    # eine Kante theoretisch möglich wäre. Da es sich um einen DBG handelt, muss
    # die Überlappung für eine Kante zwingendermaßen k-1 Basen lang sein, und da
    # nur 4 Basen (A, G, T, C) zugelassen sind, gibt es für jeden Knoten nur exakt
    # 4 andere Knoten, zu denen eine Kante möglich ist, und 4, von denen eine Kante
    # kommen kann.
    def get_potential_from(self):
        res = []
        for x in "AGTC":
            res += [x + self.seq[:-1]]
        return res

    def get_potential_to(self):
        res = []
        for x in "AGTC":
            res += [self.seq[1:] + x]
        return res

    # Methoden, um zu überprüfen, ob der Knoten mit einem nächsten/vorhergehenden
    # Knoten verschmolzen werden kann. Dies ist nur möglich, falls die Verbindung
    # in beide Richtungen eindeutig ist (beispielsweise für eine Erweiterung Richtung
    # to: self.eto darf nur eine Kante enthalten, und der Knoten, zu dem diese Kante
    # geht, darf nur einen Eintrag in self.efrom haben: Dann ist die Verbindung eindeutig
    # und die Knoten dürfen miteinander verbunden werden)
    def can_join_to(self):
        if len(self.eto) != 1:
            return False
        if len(next(iter(self.eto.keys())).efrom) != 1:
            return False
        return True
    
    def can_join_from(self):
        if len(self.efrom) != 1:
            return False
        if len(next(iter(self.efrom.keys())).eto) != 1:
            return False
        return True

    # Methoden, um den Knoten mit einem nächsten/vorhergehenden Knoten zu verbinden.
    # Falls die Verbindung möglich ist, muss die Sequenz dieses Knotens entsprechend
    # um die erste/letzte Base des zu verschmelzenden Knotens erweitert werden, und 
    # eine der Kantenlisten des zu verschmelzenden Knotens muss übernommen werden 
    # (wenn beispielsweise mit eto verschmolzen wird, muss die eto-Liste des
    # zu verschmelzenden Knotens übernommen werden - die efrom des zu verschmelzenden
    # Knotens zeigt zu self, das ist Bedingung in can_join, und diese Kante entfällt,
    # da die beiden Knoten ja verschmolzen werden).
    def join_to(self):
        if not self.can_join_to():
            return None
        nextnode = next(iter(self.eto.keys()))
        self.eto = nextnode.eto
        self.seq += nextnode.seq[-1]
        return nextnode

    def join_from(self):
        if not self.can_join_from():
            return None
        prevnode = next(iter(self.efrom.keys()))
        self.efrom = prevnode.efrom
        self.seq = prevnode.seq[0] + self.seq
        return prevnode        
    
    
# FASTA-Datei einlesen und eine Liste von Read-Objekten zurückgeben    
def read_fasta(readfile):
    f = open(readfile, "r")
    readlines = []
    reads = []
    for line in f:
        readlines += [line]
        if len(readlines) == 2:
            reads += [Read(readlines)]
            readlines = []
    return reads

# Aus einer Liste von Reads eine Dictionary von k-meren (key: k-mer-Sequenz, 
# value: Häufigkeit des Vorkommens des k-mers) erstellen. Dabei wird, um es 
# bequemer zu machen, einfach die kmer-dictionary an jeden Read übergeben und
# der Read angewiesen, seine k-mere zu der dictionary hinzuzufügen.
def build_kmerlist(reads, kmersize):
    kmers = dict()
    for read in reads:
        read.add_kmers(kmersize, kmers)
    return kmers

# Graph erstellen.
def build_graph(kmers):
    # Der Graph ist eine dictionary (key: Sequenz des Knotens, 
    # value: entsprechendes DBGnode-Objekt)
    dbg = dict()
    for kmer_s in kmers.keys():
        # Für jedes k-mer aus der Liste wird dieses zunächst, falls noch nicht
        # im DBG, dem DBG hinzugefügt
        if kmer_s not in dbg.keys():
            dbg[kmer_s] = DBGnode(kmer_s)
        kmer = dbg[kmer_s]
        # Dann wird für das k-mer überprüft, mit welchen anderen k-meren es 
        # überlappen könnte. Es gibt nur 4 mögliche ausgehende und eingehende
        # Kanten, die jeweils von get_potential_to() und get_potential_from()
        # zurückgegeben werden. Falls eins dieser k-mere tatsächlich im DBG
        # vorhanden ist, dann werden die nötigen Kanten (jeweils von und zu)
        # erstellt. Somit erfordert das Hinzufügen eines Knotens und aller
        # seiner Kanten nur 8 Abfragen, anstatt es mit allen Knoten des
        # Graphen zu vergleichen.
        for pto in kmer.get_potential_to():
            if pto in dbg.keys():
                dbg[pto].add_edge_from(kmer)
                kmer.add_edge_to(dbg[pto])
        for pfrom in kmer.get_potential_from():
            if pfrom in dbg.keys():
                dbg[pfrom].add_edge_to(kmer)
                kmer.add_edge_from(dbg[pfrom])
    return dbg

# Graph vereinfachen
def simplify_dbg(dbg):
    # Zunächst wird der vereinfachte Graph erstellt - wieder eine dictionary,
    # bei der der key die Sequenz des Knotens ist und der Value das entsprechende
    # DBGnode-Objekt.
    sdbg = dict()
    # In dieser Herangehensweise werden die Knoten der Reihe nach bearbeitet und
    # danach aus dem dbg gelöscht. Also wird so lange weitergemacht, wie im dbg noch
    # (zu verarbeitende) Knoten sind
    while len(dbg) != 0:
        # Zunächst wird irgendein Knoten aus dem dbg genommen
        currentnode = next(iter(dbg.values()))
        # Dieser wird der Liste der in diesem Schritt bearbeiteten Knoten hinzugefügt
        todelete = [currentnode.seq]
        # Solange der Knoten entweder nach to oder nach from verschmolzen werden
        # kann, wird dies getan. join_to und join_from gibt jeweils den Knoten,
        # mit dem verschmolzen wurde, zurück. Dieser wird dann jeweils auch der Liste
        # der in diesem Schritt bearbeiteten Knoten hinzugefügt.
        while currentnode.can_join_to():
            todelete += [currentnode.join_to().seq]
        while currentnode.can_join_from():
            todelete += [currentnode.join_from().seq]
        # Wenn keine weiteren Verschmelzungen mehr möglich sind, wird der (nun
        # maximal in beide Richtungen erweitere) Knoten dem vereinfachten Graph
        # hinzugefügt.
        sdbg[currentnode.seq] = currentnode
        # Danach werden alle Knoten, die in diesem Schritt bearbeitet wurden,
        # aus dem dbg gelöscht.
        for node in todelete:
            dbg.pop(node)
    return sdbg


def write_edges_file(edges):
    h = ['source', 'target']
    f = open('../Visualisierung/edges.csv', 'w')
    w = csv.writer(f)
    w.writerow(h)
    st = set(edges)
    for i in list(st):        
        w.writerow(i)    
    f.close()
    
    
def write_nodes_file(nodes):
    h = ['id']
    f = open('../Visualisierung/nodes.csv', 'w')
    w = csv.writer(f)
    w.writerow(h)
    lis = list(nodes)
    for n in lis:
        w.writerow([n])    
    f.close()
    
    
def prepare_sdbg_for_csv(sdbg):
    nodes = []
    edges = []
    for nodeobj in sdbg.values():
        nodes.append(nodeobj.seq)
        for etoseq in nodeobj.eto.keys():
            edges.append((nodeobj.seq, etoseq.seq))
    return edges, nodes


kmersize=22

# test.fasta
# virus_perfectreads.fasta
# virus_errorreads.fasta
# virus_errorreads2.fasta
# virus_perfectreads.fasta
# virus2_errorreads.fasta
# virus2_errorreads2.fasta

reads = read_fasta("virus2_errorreads.fasta")
print("Found " + str(len(reads)) + " reads.")

kmers = build_kmerlist(reads, kmersize)
print("Found " + str(len(kmers)) + " kmers.")

dbg = build_graph(kmers)
edges = 0
for x in dbg.values():
    edges += len(x.eto)
print("Graph has " + str(len(dbg)) + " nodes with " + str(edges) + " edges.")

dbglen = len(dbg)

sys.stdout.write("\n")

sdbg = simplify_dbg(dbg)
print("Simplified graph has " + str(len(sdbg)) + " nodes")

lens = [len(x) for x in sdbg.keys()]
lens.sort()
print("Contig lengths: " + ", ".join([str(x) for x in lens[::-1]]))

sdbgedges = 0
for x in sdbg.values():
    sdbgedges += len(x.eto)

print(" Bei k=" + str(kmersize) + ": Anzahl k-mere=" + str(len(kmers)) + ", Knoten im DBG=" + str(dbglen) + ", Kanten im DBG=" + str(edges) + ", Knoten im vereinfachten Graph=" + str(len(sdbg)) + ", Kanten im vereinfachten Graph=" + str(sdbgedges) + ", Längen der Contigs=" + ", ".join([str(x) for x in lens[::-1]]))

contigs = open("contigs.fasta", "w")
numcontig = 1
for contig in sdbg.keys():
    contigs.write(">contig" + str(numcontig) + "\n")
    contigs.write(contig)
    contigs.write("\n")
contigs.close()


edges, nodes = prepare_sdbg_for_csv(sdbg)

write_edges_file(edges)
write_nodes_file(nodes)