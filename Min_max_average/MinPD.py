import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade


def GetMinPD(tree, k):
    post_order_nodes = list(tree.find_clades(order='postorder'))
    node_lookup = dict((j, i) for (i, j) in enumerate(post_order_nodes))# indexes all nodes
    v = len(post_order_nodes)
    k = k

    #create array of (|v|k) size with all values = inf
    arr = np.full((v, k + 1), np.inf)
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
                mindist = np.inf
                for r in range(max(0, p - y.count_terminals()), min(x.count_terminals(), p) + 1): #goes over range for possible r values such that r+l =k
                    l = p - r #l = k-r
                    dist = arr[node_lookup[x]][r] + arr[node_lookup[y]][l] + v_x * (min(r, 1)) + v_y * (min(l, 1))#if,r is 0, wont chosose (v,x), same for l=0
                    if dist < mindist:
                        mindist = dist
                arr[node_lookup[node]][p] = mindist # fill array with value

  #  print(arr[node_lookup[tree.root]][p])
    return (arr[v-1][k])




