#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 18:51:07 2020

@author: tamilieberman
"""

'''
Create an input for use with straincraker.py

'''

# %%
import numpy as np
import argparse
# %%
def loadtaxa(nodesdmp,max_taxa):
         

    f=open(nodesdmp,"r")
    temptree=np.zeros((max_taxa,2),dtype=int)
    

    i=0
    line = f.readline().split("\t")
    print('here')
    while len(line) > 1:
        temptree[i,:]=[int(line[0]), int(line[2])]
        line = f.readline().split("\t")
        i=i+1
        
         
    temptree = temptree[:i,:]     
    f.close()

    return temptree
    
    
# %%
def buildtree(infile, outnpy, max_taxa):
    
    #load tree as child, parent list
    taxa_list = loadtaxa(infile,max_taxa)
  
    #make initial matrix with parent in 2nd row  
    max_levels=40
    max_taxas=np.max(taxa_list)
    tree=-1*(np.ones(((max_taxas+1),max_levels),dtype=int))
    tree[1,0]=1
    
    print(tree.shape)
    
    # would like to build tree from top down,
    # but sometimes parents have higher numbers than children, so need to
    # build this matrix by iteratively searching
    parents_to_populate=np.array([1]);
    i=0

    while len(parents_to_populate)>0:
        i=i+1

        parent_taxa=parents_to_populate[0]
        parent_row=tree[parent_taxa,:]
        tp=(parent_row > 0).nonzero()
        parent_level=tp[-1][-1]
        
        if (i % 10000)==0:
            print(i)
        
        children = taxa_list[taxa_list[:,1]==parent_taxa,0]
        children = children[children != parent_taxa]
        
        parents_to_populate=np.concatenate((parents_to_populate[1:],children))
        #could add in a command to delete these rows from taxa_list
        
        
        for child in children.T:
            #add children to tree
            tree[child,:]=parent_row
            tree[child,(parent_level+1)]=child
    
    np.save(outnpy, tree)
    return tree

# %%     

parser = argparse.ArgumentParser(description='Generates a treematrix for use with straincracker.py')
parser.add_argument('--i', metavar='infile', nargs='?', help='a nodes.dmp file')
parser.add_argument('--o', metavar='outfile', nargs='?', help='destination tree npy file')
parser.add_argument('--n', metavar='ntaxa', type=int, help='max number of taxa')
args = parser.parse_args()
buildtree(args.i,args.o,args.n,)

#python straincracker_loadtree.py --i nodes.dmp --n 3000000 --o treematrix.npy 



