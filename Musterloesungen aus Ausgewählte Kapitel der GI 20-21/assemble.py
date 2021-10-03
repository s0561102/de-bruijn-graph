#! /usr/bin/python3

import sys

class Read:
    def __init__(self, lines):
        self.name = lines[0].strip()
        self.bases = lines[1].strip()

    def get_kmers(self, kmersize):
        res = {}
        for pos in range(0, len(self.bases)-kmersize+1):
            kmer = self.bases[pos:(pos+kmersize)]
            if kmer not in res:
                res[kmer] = 0
            res[kmer] += 1
        return res

    def __str__(self):
        return self.name + ": " + self.bases[:20] + "..."
    
    def __repr__(self):
        return self.__str__()


class DBGnode:
    def __init__(self, seq):
        self.seq = seq
        self.eto = dict()
        self.efrom = dict()
    
    def add_edge_to(self, eto):
        if eto not in self.eto:
            self.eto[eto] = 0
        self.eto[eto] += 1
    
    def add_edge_from(self, efrom):
        if efrom not in self.efrom:
            self.efrom[efrom] = 0
        self.efrom[efrom] += 1
   
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


class DBGraph:
    def __init__(self):
        self.nodes = {}
        self.kmerlen = None
    
    def add_kmers(self, kmers):
        if len(kmers) == 0:
            return
        ## Falls kmer-Länge (Dimension) des Graphen noch nicht gesetzt ist,
        ## auf den Wert von irgendeinem k-mer aus der übergebenen dictionary
        ## setzen.
        if self.kmerlen is None:
            self.kmerlen = len(next(iter(kmers.keys())))
        for kmer_s in kmers.keys():
            if len(kmer_s) != self.kmerlen:
                raise ValueError("Incompatible k-mer lengths: " + str(self.kmerlen) + " and " + str(len(kmer_s)))
            if kmer_s not in self.nodes.keys():
                self.nodes[kmer_s] = DBGnode(kmer_s)
            kmer = self.nodes[kmer_s]
            ## Für jedes mögliche k-mer, von/zu dem es eine Kante geben könnte,
            ## die entsprechenden Kanten hinzufügen, falls dieses k-mer auch im
            ## Graphen ist (vorsicht: wenn es eine Kante A -> B gibt, sowohl in
            ## A eine kante nach B als auch in B eine Kante von A hinzufügen)
            for pto in kmer.get_potential_to():
                if pto in self.nodes.keys():
                    self.nodes[pto].add_edge_from(kmer)
                    kmer.add_edge_to(self.nodes[pto])
            for pfrom in kmer.get_potential_from():
                if pfrom in self.nodes.keys():
                    self.nodes[pfrom].add_edge_to(kmer)
                    kmer.add_edge_from(self.nodes[pfrom])
        
    def count_edges(self):
        edges = 0
        for kmer, node in self.nodes.items():
            edges += len(node.eto)
        return edges
        
    def __str__(self):
        return "DBG(" + str(self.kmerlen) + ") with " + str(len(self.nodes)) + " nodes and " + str(self.count_edges()) + " edges"
    
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

def build_graph(filename, kmersize):
    reads = read_fasta(filename)
    graph = DBGraph()
    for read in reads:
        graph.add_kmers(read.get_kmers(kmersize))
    return graph
        

dbg = build_graph("virus_perfectreads.fasta", 10)
print(dbg)
