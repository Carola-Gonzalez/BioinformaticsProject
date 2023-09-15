from Bio import SeqIO
import treeswift
import os.path
import re

fasta_filename = "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\1000,1000,.0000043,.005,long_gap_pdf,GTR+second,35,2.0,1.0\\R0\\rose.aln.true.fasta"
record_dict = SeqIO.to_dict(SeqIO.parse(fasta_filename, "fasta"))


#tree_file = "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam\\Subtree1.tree"
#alingment_file = "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam\\Subtree1Alignment"

#tree = treeswift.read_tree_newick(tree_file)

"""
for node in tree.traverse_leaves():
    gene = record_dict[node.label].seq
    header = ">" + node.label
    file = open(alingment_file, "a")
    file.write(header + "\n")
    file.write(str(gene)+"\n")
    file.write("\n")
    file.close
"""

##################################################


filepath = ("C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam", 
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam%2",
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\200Taxa.Diam",
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\200TaxaDiam%2",
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\300TaxaDiam",
                 "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\300TaxaDiam%2")  



#filepath_list= ("C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam","C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam%2" )
'''
for file in filepath: 
    for i in range(100):
        path = file + "\\Subtree" + str(i+1) +".tree"
        tree = treeswift.read_tree_newick(path)
        alingment_name = file + "\\Subtree" + str(i+1) + ".Alignment.fasta"
        print("Subtree" + str(i+1))
        
        for node in tree.traverse_leaves():
            taxa_seq = record_dict[node.label].seq
            header = ">" + node.label
            align_file = open(alingment_name, "a")
            align_file.write(header + "\n")
            align_file.write(str(taxa_seq)+"\n")
            align_file.write("\n")
            align_file.close
'''

"""
for tree in range (100):
        tree_file = file + "\\Subtree" + str(counter) + ".tree"
        tree = treeswift.read_tree_newick(tree_file)
        alingment_file = file + "\\Subtree" + str(counter) + "Alingment.fasta"

        for node in tree.traverse_leaves():
            gene = record_dict[node.label].seq
            header = ">" + node.label
            file = open(alingment_file, "a")
            file.write(header + "\n")
            file.write(str(gene)+"\n")
            file.write("\n")
            file.close
"""

# https://stackoverflow.com/questions/15166083/remove-periods-from-fasta-file-using-python

'''
test_file = "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam\\Subtree1.Alingment.fasta"
outfilename = "C:\\Users\\Carola\\OneDrive\\BRIDGE2023\\Subtrees\\100Taxa.Diam\\Subtree1.Unalinged.fasta"
regex = re.compile("[.-]+")    

with open(test_file, 'r') as infile, open(outfilename, 'w') as outfile:
    for line in infile:
        outfile.write(regex.sub("", line))
'''


for file in filepath: 
    for i in range(100):
        alingment_file_name = file + "\\Subtree" + str(i+1) + ".Alignment.fasta"
        alingment_outfile_name = file + "\\Subtree" + str(i+1) + ".Unalinged.Ufasta"
        regex = re.compile("[.-]+")  

        with open(alingment_file_name, 'r') as infile, open(alingment_outfile_name,'w') as outfile:
            for line in infile:
                outfile.write(regex.sub("", line))
