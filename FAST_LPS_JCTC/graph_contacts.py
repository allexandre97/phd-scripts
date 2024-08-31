# mypy: ignore-errors
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:12:47 2023

@author: alexandre
"""

import MDAnalysis as mda
import numpy             as np
import pickle            as pk
import matplotlib.pyplot as plt
from scipy.integrate import simps

class Contacts:

    def __init__(self,
                 RESNAME1,
                 RESNAME2,
                 GRP1_Names,
                 GRP2_Names,
                 CONTACTS):

        self.RESNAME1   = RESNAME1
        self.RESNAME2   = RESNAME2
        self.GRP1_Names = GRP1_Names
        self.GRP2_Names = GRP2_Names
        self.CONTACTS   = CONTACTS
        
def GraphContacts(CONTACTS, path):
    
    fig, ax = plt.subplots(figsize = (7, 8))

    im = ax.pcolormesh(XX, YY, CONTACTS.T,
                       cmap = 'bwr')

    ax.set_xticks(X[::2])
    ax.set_yticks(Y[::2])
    ax.set_xticklabels(GRP1_Names[::2], rotation = 45)
    ax.set_yticklabels(GRP2_Names[::2])

    ax.set_xlabel(f'{axnames[RESNAME1]}', fontsize = 20, fontweight = 'bold')
    ax.set_ylabel(f'{axnames[RESNAME2]}', fontsize = 20, fontweight = 'bold')

    ax.tick_params(labelsize = 14)

    ax2 = ax.twinx()
    ax2.set_yticks(Y[1::2])
    ax2.set_ylim(ax.get_ylim())
    ax2.tick_params(right = True, labelright = True, labelsize = 14)
    ax2.set_yticklabels(GRP2_Names[1::2])

    ax3 = ax.twiny()
    ax3.set_xticks(X[1::2])
    ax3.set_xlim(ax.get_xlim())
    ax3.tick_params(right = True, labelright = True, labelsize = 14)
    ax3.set_xticklabels(GRP1_Names[1::2], rotation = 45)

    v1 = np.linspace(CONTACTS.min(), CONTACTS.max(), 4, endpoint=True)
    cb = fig.colorbar(im, ticks=v1, ax = ax3, shrink = 1, location = 'top', pad = 0.08)
    cb.set_label(r'P(C)',
                  fontsize = 18,
                  fontweight = 'bold',
                  rotation = 0,
                  labelpad = 5)
    cb.ax.set_xticklabels(["{:4.2E}".format(i) for i in v1], fontsize='14')

    fig.tight_layout()
    plt.savefig(path,
                format = 'svg',
                dpi = 600)

def GraphCoupling(M_GRP1, M_GRP2, path1, path2):
    
    fig, ax = plt.subplots(figsize = (7, 8))

    im = ax.pcolormesh(*np.meshgrid(X, X), M_GRP1.T,
                       cmap = 'bwr')

    ax.set_xticks(X[::2])
    ax.set_yticks(X[::2])
    ax.set_xticklabels(GRP1_Names[::2], rotation = 45)
    ax.set_yticklabels(GRP1_Names[::2])

    ax.set_xlabel(f'{axnames[RESNAME1]}', fontsize = 20, fontweight = 'bold')
    ax.set_ylabel(f'{axnames[RESNAME1]}', fontsize = 20, fontweight = 'bold')

    ax.tick_params(labelsize = 14)

    ax2 = ax.twinx()
    ax2.set_yticks(X[1::2])
    ax2.set_ylim(ax.get_ylim())
    ax2.tick_params(right = True, labelright = True, labelsize = 14)
    ax2.set_yticklabels(GRP1_Names[1::2])

    ax3 = ax.twiny()
    ax3.set_xticks(X[1::2])
    ax3.set_xlim(ax.get_xlim())
    ax3.tick_params(right = True, labelright = True, labelsize = 14)
    ax3.set_xticklabels(GRP1_Names[1::2], rotation = 45)

    v1 = np.linspace(M_GRP1.min(), M_GRP1.max(), 4, endpoint=True)
    cb = fig.colorbar(im, ticks=v1, ax = ax3, shrink = 1, location = 'top', pad = 0.08)
    cb.set_label(r'$\mathbf{K^{%s}}$' % axnames[RESNAME1],
                  fontsize = 18,
                  fontweight = 'bold',
                  rotation = 0,
                  labelpad = 5)
    cb.ax.set_xticklabels(["{:4.2E}".format(i) for i in v1], fontsize='14')

    fig.tight_layout()

    plt.savefig(path1,
                format = 'svg',
                dpi = 600)
    
    fig, ax = plt.subplots(figsize = (7, 8))

    im = ax.pcolormesh(*np.meshgrid(Y, Y), M_GRP2[::-1,::-1].T,
                       cmap = 'bwr')

    ax.set_xticks(Y[::2])
    ax.set_yticks(Y[::2])
    ax.set_xticklabels(GRP2_Names[::2], rotation = 45)
    ax.set_yticklabels(GRP2_Names[::2])

    ax.set_xlabel(f'{axnames[RESNAME2]}', fontsize = 20, fontweight = 'bold')
    ax.set_ylabel(f'{axnames[RESNAME2]}', fontsize = 20, fontweight = 'bold')

    ax.tick_params(labelsize = 14)

    ax2 = ax.twinx()
    ax2.set_yticks(Y[1::2])
    ax2.set_ylim(ax.get_ylim())
    ax2.tick_params(right = True, labelright = True, labelsize = 14)
    ax2.set_yticklabels(GRP2_Names[1::2])

    ax3 = ax.twiny()
    ax3.set_xticks(Y[1::2])
    ax3.set_xlim(ax.get_xlim())
    ax3.tick_params(right = True, labelright = True, labelsize = 14)
    ax3.set_xticklabels(GRP2_Names[1::2], rotation = 45)

    v1 = np.linspace(M_GRP2.min(), M_GRP2.max(), 4, endpoint=True)
    cb = fig.colorbar(im, ticks=v1, ax = ax3, shrink = 1, location = 'top', pad = 0.08)
    cb.set_label(r'$\mathbf{K^{%s}}$' % axnames[RESNAME2],
                  fontsize = 18,
                  fontweight = 'bold',
                  rotation = 0,
                  labelpad = 5)
    cb.ax.set_xticklabels(["{:4.2E}".format(i) for i in v1], fontsize='14')

    fig.tight_layout()

    # plt.savefig(path2,
    #             format = 'svg',
    #             dpi = 600)

def GraphRelevance(rel_p, rel_r, path1, path2, path3, path4):

    fig, ax = plt.subplots(figsize = (7,7))

    sum_rel_p = np.sum(rel_p, axis = 1)

    sorted_ids_p = np.argsort(sum_rel_p)[::-1]

    for n in range(rel_p.shape[1]):
        ev = rel_p[sorted_ids_p,n]
        ax.plot(range(1, X.shape[0]+1), ev/simps(sum_rel_p[sorted_ids_p], dx = 1),
                c = 'royalblue', alpha = 0.5*n/(X.shape[0]-1))

    ax.plot(range(1, X.shape[0]+1), sum_rel_p[sorted_ids_p]/simps(sum_rel_p[sorted_ids_p], dx = 1),
            c = 'firebrick')

    ax.set_xticks(range(1, X.shape[0]+1, 2))
    ax.set_xticklabels(GRP1_Names[sorted_ids_p][::2], rotation = 45)
    ax.set_xlabel(f'{axnames[RESNAME1]} Bead Name', fontsize = 20, fontweight = 'bold')
    ax.set_ylabel(r'$\mathbf{\beta^{%s}}$' % axnames[RESNAME1], fontsize = 20, fontweight = 'bold')
    ax.tick_params(labelsize = 14)

    dummy = ax.twiny()
    dummy.set_xticks(range(2, X.shape[0]+1, 2))
    dummy.tick_params(top = True, labeltop = True, labelsize = 14)
    dummy.set_xticklabels(GRP1_Names[sorted_ids_p][1::2], rotation = 45)
    dummy.set_xlim(ax.get_xlim())

    fig.tight_layout()

    # plt.savefig(path1,
    #             format = 'svg',
    #             dpi = 600)


    fig, ax = plt.subplots(figsize = (7,7))

    sum_rel_r = np.sum(rel_r, axis = 1)

    sorted_ids_r = np.argsort(sum_rel_r)[::-1]

    for n in range(rel_r.shape[1]):
        ev = rel_r[sorted_ids_r,n]
        ax.plot(range(1, Y.shape[0]+1), ev/simps(sum_rel_r[sorted_ids_r], dx = 1),
                c = 'royalblue', alpha = 0.5*n/(Y.shape[0]-1))

    ax.plot(range(1, Y.shape[0]+1), sum_rel_r[sorted_ids_r]/simps(sum_rel_r[sorted_ids_r], dx = 1),
            c = 'firebrick')

    GRP2_Names_r = GRP2_Names[::-1]

    ax.set_xticks(range(1, Y.shape[0]+1, 2))
    ax.set_xticklabels(GRP2_Names_r[sorted_ids_r][::2], rotation = 45)
    ax.set_xlabel(f'{axnames[RESNAME2]} Bead Name', fontsize = 20, fontweight = 'bold')
    ax.set_ylabel(r'$\mathbf{\beta^{%s}}$' % axnames[RESNAME2], fontsize = 20, fontweight = 'bold')
    ax.tick_params(labelsize = 14)

    dummy = ax.twiny()
    dummy.set_xticks(range(2, Y.shape[0]+1, 2))
    dummy.tick_params(top = True, labeltop = True, labelsize = 14)
    dummy.set_xticklabels(GRP2_Names_r[sorted_ids_r][1::2], rotation = 45)
    dummy.set_xlim(ax.get_xlim())

    fig.tight_layout()

    # plt.savefig(path2,
    #             format = 'svg',
    #             dpi = 600)
    
    U_REMP = mda.Universe('./REMP.pdb', in_memory = True)
    U_PMB1 = mda.Universe('./PMB.pdb', in_memory = True)

    beta_pmb = np.zeros_like(sum_rel_p)
    for name, b in zip(GRP1_Names, sum_rel_p):
        idx = list(U_PMB1.atoms.names).index(name)
        beta_pmb[idx] = b/sum_rel_p.max()

    beta_remp = np.zeros_like(sum_rel_r)
    for name, b in zip(GRP2_Names_r, sum_rel_r):
        idx = list(U_REMP.atoms.names).index(name)
        beta_remp[idx] = b/sum_rel_r.max()

    U_REMP.add_TopologyAttr('tempfactors', beta_remp/beta_remp.max())
    U_PMB1.add_TopologyAttr('tempfactors', beta_pmb/beta_pmb.max())

    # with mda.coordinates.PDB.PDBWriter(path3) as writer:
    #     writer.write(U_PMB1)
    # writer.close()


    # with mda.coordinates.PDB.PDBWriter(path4) as writer:
    #     writer.write(U_REMP)
    # writer.close()
        

infile = 'contacts_PMB1_FEMP.pk'#sys.argv[1]

axnames = {'PMB1':'PMB',
           'REMP':'LPS',
           'FEMP':'F-LPS'}

if __name__ == '__main__':
    
    FULLDATA_CONTACTS  = []
    FULLDATA_COUPLING  = []
    FULLDATA_RELEVANCE = []
    
    for replica in ['Rep1', 'Rep2', 'Rep3']:
        
        REP_CONTACTS  = []
        REP_COUPLING  = []
        REP_RELEVANCE = []
        
        for fast in ['HalfFast', 'AllFast']:
    
            with open(f'./{replica}/{fast}/{infile}', 'rb') as f:
                CONTACTS = pk.load(f)
            f.close()
        
            RESNAME1 = CONTACTS.RESNAME1
            RESNAME2 = CONTACTS.RESNAME2
            GRP1_Names = CONTACTS.GRP1_Names
            GRP2_Names = CONTACTS.GRP2_Names
            
            mask = np.char.find(GRP2_Names.astype(str), 'C') != -1
            rev_elements = GRP2_Names[mask][::-1]
            tmp = GRP2_Names.copy()
            tmp[mask] = rev_elements
            
            new_idx = np.array([np.where(GRP2_Names == _)[0][0] for _ in tmp])
            
            GRP2_Names = tmp
            
            CONTACTS = np.array(CONTACTS.CONTACTS)
        
            dt = 5000000/CONTACTS.shape[0]
            W  = int(np.ceil(100000/dt))
        
            idx_0 = int(np.ceil(4000000/dt))
            idx_1 = int(np.ceil(5000000/dt))
        
            DATA = CONTACTS[idx_0:idx_1]
            for d in DATA:
                d /= d.sum()
        
            CONTACTS_MEAN = np.mean(DATA, axis = 0)
            
            CONTACTS_MEAN = CONTACTS_MEAN[:,new_idx]
            
            REP_CONTACTS.append(CONTACTS_MEAN)
        
            X = np.arange(CONTACTS_MEAN.shape[0])
            Y = np.arange(CONTACTS_MEAN.shape[1])
        
            XX, YY = np.meshgrid(X, Y)
            
            GraphContacts(CONTACTS_MEAN,
                          f'./{replica}/{fast}/contacts_{RESNAME1}_{RESNAME2}.svg')
            
        
            M_GRP1 = np.sqrt(CONTACTS_MEAN.dot(CONTACTS_MEAN.T))
            M_GRP2 = np.sqrt(CONTACTS_MEAN[::-1,::-1].T.dot(CONTACTS_MEAN[::-1,::-1]))
            
            REP_COUPLING.append([M_GRP1, M_GRP2])
        
            GraphCoupling(M_GRP1, M_GRP2,
                          f'{replica}/{fast}/matrix_{RESNAME1}_({RESNAME2}).svg',
                          f'{replica}/{fast}/matrix_{RESNAME2}_({RESNAME1}).svg')
        
        
            eigv_p, eigV_p = np.linalg.eigh(M_GRP1)
            eigv_r, eigV_r = np.linalg.eigh(M_GRP2)
        
            rel_p = (eigV_p*eigv_p)
            rel_r = (eigV_r*eigv_r)
            
            REP_RELEVANCE.append([rel_p, rel_r])
            
            GraphRelevance(rel_p, rel_r,
                           f'{replica}/{fast}/contribution_{RESNAME1}_({RESNAME2}).svg',
                           f'{replica}/{fast}/contribution_{RESNAME2}_({RESNAME1}).svg',
                           f'{replica}/{fast}/bindmode_{RESNAME1}_{RESNAME2}.pdb',
                           f'{replica}/{fast}/bindmode_{RESNAME2}_{RESNAME1}.pdb')
        
            
        FULLDATA_CONTACTS.append(REP_CONTACTS)
        FULLDATA_COUPLING.append(REP_COUPLING)
        FULLDATA_RELEVANCE.append(REP_RELEVANCE)
        
        
    GraphContacts(np.mean([FULLDATA_CONTACTS[0][0],
                           FULLDATA_CONTACTS[1][0],
                           FULLDATA_CONTACTS[2][0]], axis = 0),
                  f'./avg_contacts_halffast_{RESNAME1}_{RESNAME2}.svg')
    
    GraphContacts(np.mean([FULLDATA_CONTACTS[0][1],
                           FULLDATA_CONTACTS[1][1],
                           FULLDATA_CONTACTS[2][1]], axis = 0),
                  f'./avg_contacts_allfast_{RESNAME1}_{RESNAME2}.svg')
    
    GraphCoupling(np.mean([FULLDATA_COUPLING[0][0][0],
                           FULLDATA_COUPLING[1][0][0],
                           FULLDATA_COUPLING[2][0][0]], axis = 0) ,
                  np.mean([FULLDATA_COUPLING[0][0][1],
                           FULLDATA_COUPLING[1][0][1],
                           FULLDATA_COUPLING[2][0][1]], axis = 0),
                  f'./avg_matrix_halffast_{RESNAME1}_({RESNAME2}).svg',
                  f'./avg_matrix_halffast_{RESNAME2}_({RESNAME1}).svg')
    
    GraphCoupling(np.mean([FULLDATA_COUPLING[0][1][0],
                           FULLDATA_COUPLING[1][1][0],
                           FULLDATA_COUPLING[2][1][0]], axis = 0),
                  np.mean([FULLDATA_COUPLING[0][1][1],
                           FULLDATA_COUPLING[1][1][1],
                           FULLDATA_COUPLING[2][1][1]], axis = 0),
                  f'./avg_matrix_allfast_{RESNAME1}_({RESNAME2}).svg',
                  f'./avg_matrix_allfast_{RESNAME2}_({RESNAME1}).svg')
    
    GraphRelevance(np.mean([FULLDATA_RELEVANCE[0][0][0],
                            FULLDATA_RELEVANCE[1][0][0],
                            FULLDATA_RELEVANCE[2][0][0]], axis = 0),
                   np.mean([FULLDATA_RELEVANCE[0][0][1],
                            FULLDATA_RELEVANCE[1][0][1],
                            FULLDATA_RELEVANCE[2][0][1]], axis = 0),
                   f'./avg_contribution_halffast_{RESNAME1}_({RESNAME2}).svg',
                   f'./avg_contribution_halffast_{RESNAME2}_({RESNAME1}).svg',
                   f'./avg_bindmode_halffast_{RESNAME1}_{RESNAME2}.pdb',
                   f'./avg_bindmode_halffast_{RESNAME2}_{RESNAME1}.pdb')
    
    GraphRelevance(np.mean([FULLDATA_RELEVANCE[0][1][0],
                            FULLDATA_RELEVANCE[1][1][0],
                            FULLDATA_RELEVANCE[2][1][0]], axis = 0),
                   np.mean([FULLDATA_RELEVANCE[0][1][1],
                            FULLDATA_RELEVANCE[1][1][1],
                            FULLDATA_RELEVANCE[2][1][1]], axis = 0),
                   f'./avg_contribution_allfast_{RESNAME1}_({RESNAME2}).svg',
                   f'./avg_contribution_allfast_{RESNAME2}_({RESNAME1}).svg',
                   f'./avg_bindmode_allfast_{RESNAME1}_{RESNAME2}.pdb',
                   f'./avg_bindmode_allfast_{RESNAME2}_{RESNAME1}.pdb')