# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 16:06:07 2023

@author: alexandre
"""

import numpy             as np
import pickle            as pk
import networkx          as nx
import matplotlib.pyplot as plt

from optparse                       import OptionParser
from matplotlib                     import cm
from networkx.algorithms.components import connected_components


parser = OptionParser(usage       = 'python graph_clusters.py',
                      prog        = 'GraphClusters',
                      description = 'This is a program used to graph clustering of molecules.') 

parser.add_option('-i', '--inp', 
                  action = 'store', type = 'string',
                  default = 'clusters_PMB1_PMB1.pk',
                  help = 'Input File')

parser.add_option('-o', '--out', 
                  action = 'store', type = 'string',
                  default = 'clusters',
                  help = 'Master Output Filename')

parser.add_option('-b', '--beg', 
                  action = 'store', type = float,
                  default = 0,
                  help = 'Start Time in ns')

parser.add_option('-e', '--end', 
                  action = 'store', type = float,
                  default = 5000,
                  help = 'End Time in ns')

options, args = parser.parse_args()


class Clusters:
    
    def __init__(self,
                 GroupA : str, GroupB : str,
                 Time : np.array,
                 ContactMatrix : np.array):
        
        self.GroupA        = GroupA
        self.GroupB        = GroupB
        self.Time          = Time
        self.ContactMatrix = ContactMatrix


def BuildGraph(Contacts : np.array,
               Cutoff   : float):
    
    G = nx.Graph()
    
    for i in range(Contacts.shape[0]):
        G.add_node(i)
    
    for i in range(Contacts.shape[0]):
        where_contact = np.where(Contacts[i] != 0)[0]
        for j in where_contact:
            if Contacts[i, j] > Cutoff:
                G.add_edge(i, j, weight = Contacts[i,j])
    
    return G

def CountClusters(Subgraphs : list):
    
    N        = 0
    Nodes    = []
    AvDegree = []
    
    for g in Subgraphs:
        nnodes = len(g.nodes())
        if nnodes > 1:
            N += 1
            dens    = nx.density(g)
            AvDegree.append(dens)
        Nodes.append(nnodes)    
        
    return Nodes, N, AvDegree

def NNodeDistribution(NNodes):
    
    N, bins = np.histogram(NNodes, np.arange(2, 21), density = True)
    
    return N, bins


if __name__ == '__main__':
    
    for replica in ['Rep1', 'Rep2', 'Rep3']:
        
        for fast in ['NoFast', 'HalfFast', 'AllFast']:
    
            with open(f'./{replica}/{fast}/{options.inp}', 'rb') as f:
                CLUSTERS = pk.load(f)
            f.close()
    
            ContactMatrix = CLUSTERS.ContactMatrix
            Time          = CLUSTERS.Time
            GroupA        = CLUSTERS.GroupA
            GroupB        = CLUSTERS.GroupB
            
            assert options.beg*1000 >= Time[0]
            assert options.end*1000 <= Time[-1]
            
            dt = Time[1] - Time[0]
            
            idx_0 = int(np.ceil(1000*options.beg/dt))
            idx_1 = int(np.ceil(1000*options.end/dt))
            
            Time = Time[idx_0:idx_1+1]
            
            n_frames = idx_1-idx_0-1
            
            PNNodes_data   = []
            NClusters_data = []
            
            for f, frame in enumerate(ContactMatrix[idx_0:idx_1+1]):
                
                Graph = BuildGraph(frame, 12)
            
                subgraphs  = [Graph.subgraph(c) for c in connected_components(Graph)]
                NNodes, NClusters, AvDegree  = CountClusters(subgraphs)
                
                PNNodes, bins = NNodeDistribution(NNodes)
                
                PNNodes_data.append(PNNodes)
                
                NClusters_data.append([NClusters, np.mean(AvDegree)])
                
                print(f'Done {int(100 * f/(n_frames))} %', end = '\r')
            
            PNNodes_data   = np.array(PNNodes_data)
            NClusters_data = np.array(NClusters_data)
            
            
            TT, BB = np.meshgrid(Time/1000, bins[:-1])
        
            fig, ax = plt.subplots(figsize=(7, 7))
            im = ax.pcolormesh(TT, BB, PNNodes_data.T,
                               cmap=cm.bwr, shading = 'nearest')
        
            ax.set_xlabel('Time / ns',
                          fontsize=20, fontweight='bold')
            ax.set_ylabel('Cluster Size',
                          fontsize=20, fontweight='bold')
            ax.set_yticks(bins[:-1])
            ax.tick_params(labelsize=18)
        
            ax2 = ax.twinx()
            ax2.set_yticks([])
            ax2.set_yticklabels([])
            v1 = np.linspace(PNNodes_data.min(), PNNodes_data.max(), 4, endpoint=True)
            cb = fig.colorbar(im, ticks=v1, ax=ax2, shrink=1, location='top')
            cb.set_label(r'P(Size)',
                         fontsize=18,
                         fontweight='bold',
                         rotation=0,
                         labelpad=0)
            cb.ax.set_xticklabels(["{:4.2f}".format(i) for i in v1], fontsize='14')
        
            fig.tight_layout()
            
            plt.savefig(f'{replica}/{fast}/{options.out}_sizeprob_{options.beg}ns_{options.end}ns_2.svg',
                        format = 'svg',
                        dpi = 600)
            
            
            
            fig, ax = plt.subplots(figsize=(7, 7))
        
            ax.scatter(Time/1000, NClusters_data[:, 0],
                       c=(NClusters_data[:, 1]/NClusters_data[:, 1].max()),
                       cmap=cm.bwr)
            ax.set_xlabel('Time / ns',
                          fontsize=20, fontweight='bold')
            ax.set_ylabel('# Clusters',
                          fontsize=20, fontweight='bold')
            ax.tick_params(labelsize=18)
        
            fig.tight_layout()
            
            #plt.savefig(f'{options.out}_nclust_{options.beg}ns_{options.end}ns.png')
            

            with open(f'./{replica}/{fast}/{options.out}_nclust_{options.beg}ns_{options.end}ns_2.pk', 'wb') as f:
                ncl = np.zeros((Time.shape[0], 2))
                ncl[:,0] = Time
                ncl[:,1] = NClusters_data[:,0]
                pk.dump(ncl, f,
                        protocol = 4)
            f.close()
            
            
            with open(f'./{replica}/{fast}/{options.out}_avnei_{options.beg}ns_{options.end}ns_2.pk', 'wb') as f:
                
                avn = np.zeros((Time.shape[0], 2))
                avn[:,0] = Time
                avn[:,1] = NClusters_data[:,1]
                
                pk.dump(avn, f,
                        protocol = 4)
            
            f.close()
