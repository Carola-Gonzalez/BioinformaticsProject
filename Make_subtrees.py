import argparse
import numpy
import treeswift
import sys
import pandas as pd
import os 
import statistics
from Bio import SeqIO
import re


# Give Distance matrix of leaves and range of diameter prefered 
# Goes through each cell and check is distance is in diameter range, appends it to pairs list 
# Will go though hald the matrix to not have duplicate pairs in list
# Takes in distance matrix from treeswift (made into dataframe) and diameter  
def get_pairs(df, diam): 
    
    diameter_lower = diam/2
    diameter_upper = diam

    taxa_arr = df.columns.to_numpy()
    pairs = []
    for i in taxa_arr: 
        for j in taxa_arr: 

            if i<j: 
                diameter = df.loc[i][j]
                if diameter > diameter_lower and diameter < diameter_upper:
                    pairs.append([i,j,])
    return pairs


# Gets all leaves of the clade in which the pair is in 
# Get MRCA for a and b, build subtree rooted at ancestor
# Traverse tree subtree and get leaves
# Takes in list of pairs from get pairs method, original tree and mrca matrix from treeswift(made into dataframe) 
def get_clade_leaves(pair, mdf, tree):

    f = pair[0]
    s = pair[1]
    ancestor = mdf.loc[s][f]

    subtree = tree.extract_subtree(ancestor)
    leaves = get_leaf_labels(subtree)

    return leaves


# Switch leaf labels to node object reference
# Takes in a list of labels and dictonary mappin label:reference 
def switch_leaves(leaves, dict):
    
    new_leaf_set = []

    for leaf in leaves: 
        x = dict.get(leaf)
        new_leaf_set.append(x)

    return new_leaf_set


# Switch from object reference to node label
# Takes in list of taxa 
def node_to_label(taxa):
    l = []
    for node in taxa: 
        l.append(node.label)

    return l


# Makes dictonary mapping node label to its object reference
# Takes in distance matrix dataframe
def make_dict(df):
    dict = {}
    df_names = list(df)

    for name in df_names:
        dict[name.label] = name

    return dict


# Input pairs from get pairs method, distance matrix dataframe, k taxa wanted, oringal tree...
# MRCA matrix dataframe and dictornary mapping node.labels:node reference 
# Goes through every pair and get the clade that it is in 
# To the orignal pair, it adds leaves from clade that maintain the maximum diameter wanted (diam from pair)
# If new taxa list is < k, dont do anythign with it, go to next pair 
# Once all possibe leaves have been added, extract a tree from the new taxa list 
# From the new tree, maximize for pd using fast pd algorithm
# Get leaves that maximize pd and extract tree from it 
# Turn tree into newick string and check in tree duplicate is already in subtree list 
# If not append 
def subtrees(pairs, df, k, tree, mdf, dict, subtree_num):

    tree_list = []

    for taxa in pairs:

        temp = get_clade_leaves(taxa,mdf,tree)
        leaves = switch_leaves(temp,dict)
        
        add_leaves(taxa,df,leaves)
      
        if len(taxa) == k: 
           tree_leaves = node_to_label(taxa)
           subtree = tree.extract_tree_with(tree_leaves)
           nw = subtree.newick()

           if check_if_duplicate(nw,tree_list) == True: 
                tree_list.append(nw)

           if len(tree_list)==subtree_num: 
                return tree_list

        elif len(taxa)>k: 
            tree_leaves = node_to_label(taxa)
            temp_subtree = tree.extract_tree_with(tree_leaves)
            kleaves, pd = fast_pd(temp_subtree,k)
            kleaves = node_to_label(kleaves)
            subtree = tree.extract_tree_with(kleaves)
            nw = subtree.newick()

            if  check_if_duplicate(nw,tree_list) == True: 
                tree_list.append(nw)

            if len(tree_list)==subtree_num: 
                return tree_list
  
    return tree_list


# Check if newick string in already in list 
def check_if_duplicate(nw,nw_list):
    for tree in nw_list:
        if nw == tree: 
            return False
    return True


# Adds leaves from clade to taxa list (taxa list being the two leaves that make the pair at first)
# Input taxa (pair), distance matrix, leaves(from subtree made from mrca of pair)
# For each leaf in leaves, only add it to taxa list if its distance to any node in taxa < max_diam (pair diameter)  
def add_leaves(taxa,df,leaves):

    max_diam = df.loc[taxa[0]][taxa[1]]

    for node in leaves: 

            comp = compare_with(node,taxa,df,max_diam)

            if comp == True: 
                taxa.append(node)


# Compare specific leaf from clade, to all the leaves in taxa list 
# Node is leaf from clade, taxa is the leaves with max diameter being that if the first pair, 
# df is distance matrix dataframe and diam is the diameter of the first pair 
def compare_with(node,taxa,df,max_diam):
    
    for leaf in taxa:         
        if df.loc[node][leaf]> max_diam or node==leaf:
            return False 
        
    return True


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


# Sum all the branch weights in this tree
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


# Sort leaves of tree in ascending order 
def sort_leaves(lset): 
    lset.sort(key=lambda x: x.edge_length, reverse=False)
    return lset


# Get degree of specific node 
def get_degree(node): 
    degree = len(node.child_nodes()) + 1
    return degree


# Computes k amount of leaves that maximize pf 
def fast_pd(tree, k): 

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


# Returns a list of subtrees with given diameter and maximize for k taxa 
def subtrees_with_DK(newick,diam,k, subtree_num):
    
    tree = treeswift.read_tree_newick(newick)

    disM = tree.distance_matrix()
    df = pd.DataFrame(disM)
    df = df.fillna(0)

    mdf = tree.mrca_matrix()
    mdf = pd.DataFrame(mdf)
    mdf = mdf.fillna(0)

    dict = make_dict(df)

    pairs = get_pairs(df,diam)

    trees = subtrees(pairs,df,k,tree,mdf,dict,subtree_num)

    return trees

# This function takes a list of newicks strings produced by the subtree function
# And a the a save file and makes a file for each tree
# Make sure to put subtrees with different taxa or diameter in different files
def save_subtrees(subtrees,save_path):
    for tree in subtrees: 
        name_of_file = "Subtree" + str(counter)
        completeName = os.path.join(save_path, name_of_file+".tree")
        file1 = open(completeName, "w+")
        toFile = tree    
        file1.write(toFile)   
        file1.close()
        counter+=1


# Gets the biggest diameter of the tree 
def get_diameter(tree):
    lset, diam = fast_pd(tree,2)
    return diam


# Get smallest branch of the internal branches within a tree
def smallest_branch(tree):
    smallest = 1000

    for node in tree.branch_lengths(terminal=False, internal = True):
        if node != 0 and node != None: 
            if node < smallest: 
                smallest = node 
               
    return smallest
           

# Get biggest branch of the internal branches within a tree
def biggest_branch(tree):
    biggest = 0

    for node in tree.branch_lengths(terminal=False, internal = True):
        if node != 0 and node != None: 
            if node > biggest: 
                biggest = node 
               
    return biggest 


# Get Std of the internal branches of a tree
def get_branch_std(tree):
    sample = []

    for node in tree.branch_lengths(terminal=False, internal = True):
        if node!=0:
            sample.append(node)
    
    std = statistics.stdev(sample)
    return std 


# Function to get file of the stdev of each tree. Indexed by line. 
# Input is a list of file paths with trees of different taxa num and diameter 
def std_of_trees(filepath_list):
    
    for path in filepath_list:
    
        datapath = path + "\\Stdev_of_Branches.txt"
        
        file = open(datapath, "w+")

        for lp in range(100):
            filename = path + "\\Subtree" + str(lp+1)+".tree"
            tree = treeswift.read_tree_newick(filename)
            
            data = get_branch_std(tree)
            
            file.write(str(data)+"\n")
            
        file.close()


# Function to get file of the smallest branch of each tree. Indexed by line. 
def smallest_branches(filepath_list):
    
    for path in filepath_list:
    
        datapath = path + "\\Smallest_Branches.txt"
        
        file = open(datapath, "w+")

        for lp in range(100):
            filename = path + "\\Subtree" + str(lp+1)+".tree"
            tree = treeswift.read_tree_newick(filename)
            
            data = smallest_branch(tree)
            
            file.write(str(data)+"\n")
            
        file.close()


# Function to get file of the diameter of each tree. Indexed by line. 
def diameter_of_trees(filepath_list):
    
    for path in filepath_list:
    
        datapath = path + "\\DiamCalc.txt"
        
        file = open(datapath, "w+")

        for lp in range(100):
            filename = path + "\\Subtree" + str(lp+1)+".tree"
            tree = treeswift.read_tree_newick(filename)
            
            data = get_diameter(tree)
            
            file.write(str(data)+"\n")
            
        file.close()


# Function to get file of the biggest branch of each tree. Indexed by line. 
def biggest_branches(filepath_list):
    
    for path in filepath_list:
    
        datapath = path + "\\Biggest_Branches.txt"
        
        file = open(datapath, "w+")

        for lp in range(100):
            filename = path + "\\Subtree" + str(lp+1)+".tree"
            tree = treeswift.read_tree_newick(filename)
            
            data = biggest_branch(tree)
            
            file.write(str(data)+"\n")
            
        file.close()


# Function to get file of the avg branch of each tree. Indexed by line. 
def avg_branches(filepath_list):
    
    for path in filepath_list:
    
        datapath = path + "\\Avg_Branches.txt"
        
        file = open(datapath, "w+")

        for lp in range(100):
            filename = path + "\\Subtree" + str(lp+1)+".tree"
            tree = treeswift.read_tree_newick(filename)
            
            data = tree.avg_branch_length()
            
            file.write(str(data)+"\n")
            
        file.close()


# Takes in a list of file paths (can not be just 1) and makes fasta file for each tree in the files
# Subtrees must be names the same in each file (Subtree1.tree , Subtree2.tree, Subtree3.tree...etc)
# Two subtree with different taxa (100vs200) need to be named the same for the function to work 
# What seperates them is the parent file they belong to 
def make_fasta_files(original_fasta, filepath_list):

    record_dict = SeqIO.to_dict(SeqIO.parse(original_fasta, "fasta"))

    for file in filepath_list:

        for i in range(100):
            path = file + "\\Subtree" + str(i+1) +".tree"
            tree = treeswift.read_tree_newick(path)
            alingment_name = file + "\\Subtree" + str(i+1) + ".Alignment.fasta"
            
            for node in tree.traverse_leaves():
                taxa_seq = record_dict[node.label].seq
                header = ">" + node.label
                align_file = open(alingment_name, "a")
                align_file.write(header + "\n")
                align_file.write(str(taxa_seq)+"\n")
                align_file.write("\n")
                align_file.close


# https://stackoverflow.com/questions/15166083/remove-periods-from-fasta-file-using-python
# This function makes unaligned fasta files from the aligment ones, therefore, 
# the aligment fasta files must exist and be within the same folder as the trees 
def make_Ufasta_files(filepath_list):

    for file in filepath_list: 

        for i in range(100):
            alingment_file_name = file + "\\Subtree" + str(i+1) + ".Alignment.fasta"
            alingment_outfile_name = file + "\\Subtree" + str(i+1) + ".Unalinged.Ufasta"
            regex = re.compile("[.-]+")  

            with open(alingment_file_name, 'r') as infile, open(alingment_outfile_name,'w') as outfile:
                for line in infile:
                    outfile.write(regex.sub("", line))


