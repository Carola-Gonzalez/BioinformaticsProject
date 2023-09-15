import treeswift

def get_leaf_labels(tree):
    lset = []
    for leaf in tree.traverse_leaves():
        lset.append(str(leaf.label))
    return lset

def get_leaf_set(tree):
    lset = []
    for leaf in tree.traverse_leaves():
        lset.append(leaf)
    return lset

# Sort leaves of tree in ascending order 
def sort_leaves(lset): 
    lset.sort(key=lambda x: x.edge_length, reverse=False)
    return lset

def print_leaves(lset):
    for node in lset:
        print(node.label) 

def get_degree(node): 
    degree = len(node.child_nodes()) + 1
    return degree

def print_taxa(taxa_set):
    for node in taxa_set: 
        print(node.label)
    
def sum_edges(tree):
    esum = 0
    for node in tree.traverse_postorder():
        if node.is_root():
            pass
        else:
            elen = node.get_edge_length()
            if elen is not None: 
                esum += elen
    return esum

def fast_pd(tree, k): 
    print("HELLOOOO")
    #get lset 
    lset = get_leaf_set(tree)
    #print_taxa(lset)
    
    if len(lset)==k: 
        return tree

    #Sort lset by edge weigth 
    sort_leaves(lset)

    while len(lset) > k: 
        parent = lset[0].get_parent()
        parent.remove_child(lset[0]) 
        parent_degree = get_degree(parent) 
        lset.pop(0)
       
        
        if parent_degree < 3: 
            parent.contract()
            sort_leaves(lset)

    pd = sum_edges(tree)
    return lset, pd

# Labels all internal nodes V1, V2, V2 etc. Root node is designed as V0 and its tw as 0
def label_internal_nodes(tree):
    counter = 0
    for node in tree.traverse_internal(): #traverses the root node, we dont want that 
        node.label = "V" + str(counter)
        counter += 1

# Gets root of tree by traversing in preorder. 
# First node traversed will the root so it should not be computationally taxing
def get_root(tree):
    for node in tree.traverse_preorder():
        if node.is_root():
            node.edge_length = 0 
            node.tw = 0
            return node

def leaf_num(tree):
    leaves = 0
    for node in tree.traverse_leaves():
        leaves+=1
        
    return leaves