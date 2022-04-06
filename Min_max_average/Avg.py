import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
from scipy.special import comb as choose
from scipy.special import logsumexp as fix
import math
from math import log
from scipy.special import gammaln



def GetAvgPD(tree, k):
    post_order_nodes = list(tree.find_clades(order='postorder'))
    node_lookup = dict((j, i) for (i, j) in enumerate(post_order_nodes))# indexes all nodes
    v = len(post_order_nodes)
    k = k

    #create array of (|v|k) size with all values = inf
    arr = np.full((v, k + 1), 0.0)
    arr[:, 0] = 0 #base case 1- for k=0, d(v,k)=0 regardless of vertex

    #calculate D(v,k) for all nodes
    for node in post_order_nodes:
        #base case 2- D(v,k) = 0 for all leaf nodes
        if node.is_terminal():
            arr[node_lookup[node]] = 0
            # D(v.k) for all k = 0
        else:
            #recursive case
            x = node.clades[0] #right tree
            y = node.clades[1] #left tree
            v_x = tree.distance(node, x) #edge (v,x)
            v_y = tree.distance(node, y) #edge (v,y)
            for p in range(min(k, node.count_terminals()) + 1): # cant select more than |v| nodes at that vertex , thus min(k,|v|)
                #for average the extreme cases always exist that is right and left subtree are always chosen
                sum_x = choose(x.count_terminals(), p) * (v_x + arr[node_lookup[x]][p])
                sum_y = choose(y.count_terminals(), p) * (v_y + arr[node_lookup[y]][p])
                sum_avg = sum_x + sum_y
                mindist = -1
                for r in range(max(0, p - y.count_terminals()), min(x.count_terminals(), p) + 1): #goes over range for possible r values such that r+l =k
                    l = p - r #l = k-r
                    sum_avg += choose(x.count_terminals(),r) * choose(y.count_terminals(),l) *(arr[node_lookup[x]][r]+ arr[node_lookup[y]][l]+v_x+v_y)#+= since avg formula requires addition of all

                arr[node_lookup[node]][p] = sum_avg*(min(1,p)) / (choose(node.count_terminals(),p)*1.0) #normalize over set size

  #  print(arr[node_lookup[tree.root]][p])
    return (arr[v-1][k])
