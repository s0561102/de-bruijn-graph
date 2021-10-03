#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 15:28:50 2021
Source: https://towardsdatascience.com/genome-assembly-using-de-bruijn-graphs-69570efcc270
@author: tm
"""
from Bio import SeqIO
import csv
import hashlib


def debruijnize(k_mers):
    nodes = set()
    edges = []
    for k, w in k_mers.items():
        r1 = k[:-1]
        r2 = k[1:]
        #r1_hash = kmer_to_sha(r1)
        #r2_hash = kmer_to_sha(r2)
        nodes.add(r1)
        nodes.add(r2)
        edges.append((r1,r2))
    
    return (nodes,edges)

    
def kmer_to_sha(kmer):
    return hashlib.sha224(kmer.encode('utf-8')).hexdigest()


def build_k_mer(str,k):
    return [str[i:k+i] for i in range(0,len(str)-k+1)]


def read_fasta_reads_to_kmers(filename, k_mer_size):
    k_mers = []
    weighted = {}
    for seq_record in SeqIO.parse(filename, "fasta"):
        k_mers.extend(build_k_mer(str(seq_record.seq), k_mer_size))
        
    for i in k_mers:
        if i not in weighted:
            weighted[i] = 0
        weighted[i] += 1           
    return weighted


def de_bruijn_graphizer(filename, k):
    k_mers = {}
    k_mers = read_fasta_reads_to_kmers(filename, k)
    G = debruijnize(k_mers)
    write_edges_file(G[1])
    write_nodes_file(G[0])
    
    
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
    
    
def print_list(items):
    for i in items:
        print(i)
        
        
de_bruijn_graphizer("test.fasta", 10)