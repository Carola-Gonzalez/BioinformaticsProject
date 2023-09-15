import treeswift
import os.path
import Utils
import statistics

"""
For each folder (4 of them)....
    Creates a file where you Will put the summary stat 

    For each tree in folder (100)
        Go throuth each tree and calculate the summary stat (100 trees)
        Put in the summary state txt
    
    close summary stat file 

"""

"""
What summary stats are we looking at? 
 - Number of taxa 
 - Diameter (fast_pd(tree,2))
 - Mean of branches
 - Biggest Branch
 - Smallest Branch 
 - Standard deviation 

"""

# Gets the biggest diameter of the tree 
def get_diameter(tree):
    lset, diam = Utils.fast_pd(tree,2)
    return diam



# Get smallest branch 
def smallest_branch(tree):
    smallest = 1000

    for node in tree.branch_lengths(terminal=False, internal = True):
        if node != 0 and node != None: 
            if node < smallest: 
                smallest = node 
               
    return smallest
           

# Get biggest branch 
def biggest_branch(tree):
    biggest = 0

    for node in tree.branch_lengths(terminal=False, internal = True):
        if node != 0 and node != None: 
            if node > biggest: 
                biggest = node 
               
    return biggest 

# Get Std 
def get_branch_std(tree):
    sample = []

    for node in tree.branch_lengths(terminal=False, internal = True):
        if node!=0:
            sample.append(node)
    
    std = statistics.stdev(sample)
    return std 

###########################################
# Testing 1

test_tree1 = treeswift.read_tree_newick("C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam\\Subtree1.tree")
test_tree2 = treeswift.read_tree_newick("C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam\\Subtree2.tree")
test_tree3 = treeswift.read_tree_newick("C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam\\Subtree3.tree")
taxa_200 = treeswift.read_tree_newick("C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\200TaxaDiam%2\\Subtree10.tree")

#diam = get_diameter(test_tree1)
#print(diam)

#diam2 = get_diameter(test_tree2)
#print(diam2)

#branch_avg = test_tree.avg_branch_length(terminal=False, internal=True)
#print(branch_avg)

#smallest = smallest_branch(test_tree)
#print(smallest)

#biggest = biggest_branch(test_tree)
#print(biggest)

#std = get_branch_std(test_tree)
#print(std)

for node in taxa_200.branch_lengths(terminal=False, internal=True):
    print(node)

small = smallest_branch(taxa_200)
print("SMALLEST BRANCH")
print(small)
#####################################
#Testing 2


filepath_list = ("C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam", 
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam%2",
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\200Taxa.Diam",
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\200TaxaDiam%2",
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\300TaxaDiam",
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\300TaxaDiam%2")  




"""


for path in filepath_list:
    counter = 1
    
    datapath = path + "\\Stdev_of_Branches.txt"
    
    file = open(datapath, "w+")

    for lp in range(100):
        filename = path + "\\Subtree" + str(counter)+".tree"
        tree = treeswift.read_tree_newick(filename)
        
        data = get_branch_std(tree)
        
        file.write(str(data)+"\n")
        
        counter+=1
    
    file.close()
    
"""





"""
    for lp in range(100):
        filename = "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam\\Subtree" + str(counter)+".tree"
        tree = treeswift.read_tree_newick(filename)
        diam = get_diameter(tree)
        data.append(diam)
        counter+=1
print(data)
"""




"""
for lp in len(filepath_list):

    # Create file for stat 

    for lp in range(100):
        # calculate stat for each tree 
        # write to file 

"""

"""
savepath = "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam"
nameoffile = "DiamCalc"

completeName = os.path.join(savepath,nameoffile +".txt")
file1 = open(completeName, "w+")
toFile = "Diam for Subtree 1: " + str(diam)
file1.write(toFile)
file1.close()
"""




