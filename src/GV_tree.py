import pandas as pd
import numpy as np
from ete3 import Tree
from graphviz import Source
import re
import os
import argparse


argParser = argparse.ArgumentParser(prog='PROG')
argParser.add_argument('-output','--outputprefix', type=str,default = 'out')  # CHANGE HERE



args = argParser.parse_args()


data = args.outputprefix


with open(data+'inferred_tree.nw','r') as f:
        string = f.read()



string = re.sub('c','',string)
with open('inferred_format.nw','w') as f:
        f.write(string)


def decrease_numbers(text):
    def replace(match):
        num = match.group(0)
        return str(int(num)-1)

    pattern = r'\b\d+\b'
    return re.sub(pattern, replace, text)

# Read the text file
file_path = 'inferred_format.nw'  # Update this with the path to your text file
with open(file_path, 'r') as file:
    original_text = file.read()

# Decrease numbers in the text
modified_text = decrease_numbers(original_text)

# Write the modified text back to the file
with open(file_path, 'w') as file:
    file.write(modified_text)
    file.write(";")

with open('inferred_format.nw','r') as f:
        content = f.read()
        #print(content)


        



# CHANGE THE TREE FROM HERE
# Use cluster information from outside the folder
tree_nw = Tree("inferred_format.nw")

def findPar(tree_nw, split_list):
   for node1 in tree_nw.traverse("preorder"):
     des = list( j for j in node1.iter_descendants() )  # descendants of current node
     ans = list( h for h in node1.get_ancestors() )  # anscestors of current node
     for node2 in tree_nw.traverse("preorder"):
      if len(ans)==0:
       break
      else:
       if node1 != node2 and node2 not in ans:
         if node2 not in des and node1.up!=node2.up:
           if node1.is_leaf():
             mutated_names1 = set([node1.name])
           else:
             mutated_names1 = set(node1.get_leaf_names(is_leaf_fn=None))
           if node2.is_leaf():
             mutated_names2 = set([node2.name])
           else:
             mutated_names2 = set(node2.get_leaf_names(is_leaf_fn=None))
           m = mutated_names1.union( mutated_names2)
           if set((np.array(list(m),dtype='int'))) == set([i for i,val in enumerate(split_list) if val==1]):
            return node1, node2

def findLoss(tree_nw, split_list):
   for node1 in tree_nw.traverse("preorder"):
     des = list( j for j in node1.iter_descendants() )  #list of descendants of current node
     ans = list( h for h in node1.get_ancestors() )  #list of anscestors of current node
     for node2 in tree_nw.traverse("preorder"):
      if len(des)==0:
       break
      else:
       if node1 != node2 and node2 not in ans:
         if node2 in des:
           if node1.is_leaf():
             mutated_names1 = set([node1.name])
           else:
             mutated_names1 = set(node1.get_leaf_names(is_leaf_fn=None))
           if node2.is_leaf():
             non_mutated_names2 = set([node2.name])
           else:
             non_mutated_names2 = set(node2.get_leaf_names(is_leaf_fn=None))
           m = mutated_names1 - non_mutated_names2
           if set((np.array(list(m),dtype='int'))) == set([i for i,val in enumerate(split_list) if val==1]):
            return node1, node2

def findNode(tree_nw, split_list, config):
  if config == 'Het' or config == 'Clonal':
    mutated = [i for i,val in enumerate(split_list) if val==1]
    if len(mutated) == 1:
      return tree_nw & str(mutated[0])
    return tree_nw.get_common_ancestor([str(i) for i in mutated])
  elif config == 'Par':
    mutated = [i for i,val in enumerate(split_list) if val==1]
    node1, node2 = findPar(tree_nw, split_list)
    return node1, node2
  elif config == 'Loss' or config == 'Back':
    mutated = [i for i,val in enumerate(split_list) if val==1]
    node1, node2 = findLoss(tree_nw, split_list)
    return node1, node2
  else:
    print("Problem")
    return 0

   

orig_pos = pd.read_csv(data+'final_df.tsv',sep='\t',header=None) # this is the dataframe which is obtained after post processing and everything

config = pd.read_csv('Genotype configuration.tsv',sep='\t',header=None) # Genotype configuration file
before_post_process_file = pd.read_csv(data+'final_df.tsv',sep='\t',header=None) # final data frame for finding the sites which are mutated
numcells = config.shape[1]-1
print("Number of cells :",numcells)

config.insert(loc = 0,column = 'chr', value = before_post_process_file[0])  # putting chromosome information at first column
config.insert(loc = 1,column = 'site',value = before_post_process_file[1]) # putting the site information at the second column

config = config[config['site'].isin(orig_pos[1])]

print("Shape of configuration file:",config.shape)
cols = [i for i in range(0,numcells)]

config['sum'] = config[cols].sum(axis=1) # sum across all cells
config = config[config['sum'] != 0]  # if sum is more than 0, keep those sites i.e. if more than 0 cell is mutated i.e. at least one cell should be mutated
config.reset_index(inplace=True,drop=True)
config  = config.drop(['sum'],axis=1)

config2 = config.T.drop_duplicates().T # selecting the first occurence of a different genotype across all cells, i.e. finding unique genotypes  which will readily convert the cellular dataframe into clonal one

cluster_labels = list(np.loadtxt(data+'inferred_cluster_labels.txt',dtype='int')) # now finding the labels corresponding to each unique cell
columns = dict()
for i in config2.columns[2:-1]:
  columns[i] = cluster_labels[i]

no_of_clusters = len(set(cluster_labels))
print("Number of clusters:",no_of_clusters)

missing = 0

for cluster_no in range(no_of_clusters):
  if cluster_no not in columns.values():
    cell_instance = cluster_labels.index(cluster_no) # index at which the first instance of that cluster appears in the crowd of all cells
    config2.insert(loc = config2.shape[1]-1,column = cell_instance,value = config[cell_instance])
    columns[cell_instance] = cluster_no
    missing += 1
    # find an instance of that cluster number in the cell configuration and add in the config


config2.rename(columns = columns, inplace = True)  # replace the cell names with the clonal information obtained just above



cols = []
cols.append('chr')
cols.append('site')
for i in range(len(set(cluster_labels))):
  cols.append(i)
cols.append(len(cluster_labels))  # reorder the columns so that there is no confusion 

config2 = config2[cols]
#config2 : clonal level configuration ready

clonal_info = np.array(config2,dtype='object')  # into numpy format for easy element access


i=0
for node in tree_nw.traverse("preorder"): # naming the nodes because by default they are just empty strings
    node.temp = i+len(set(cluster_labels))
    i+=1
for node in tree_nw.traverse("preorder"): # arrow here refers to the arrow we use to tell the relationship between two nodes, so it basically has the 
                                          # information about the children of an internal node
    node.arrow = []
for node in tree_nw.traverse("preorder"): # initialise the labels with an empty list
    node.label = []
for i in range(clonal_info.shape[0]): # every site individually
    chr = clonal_info[i][0][3:]
    site = clonal_info[i][1]
    config = clonal_info[i][-1]
    split_list = clonal_info[i][2:-1]
    if config == 'Clonal' or config == 'Het':
      node = findNode(tree_nw, split_list, config)
      node.label.append(str(chr) + ":" + str(site))
    elif config == 'Par' or config == 'Back' or config == 'Loss':
      node1,node2 = findNode(tree_nw, split_list, config)
      node1.label.append(str(chr) + ":" + str(site))
      if config == 'Back' or config == 'Loss':
        node2.label.append(str(chr) + ":" + str(site) + '-')
      else:
        node2.label.append(str(chr) + ":" + str(site))
    else:
      print("Fault")
#for node in tree_nw.traverse("preorder"):
  #print(node.label)


for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    children = node.get_children()
    node.arrow.extend(children)
string = ''
for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    k=0
    for j in range(len(node.get_children())):
      string_k = str(node.temp-len(set(cluster_labels))) + ' -> ' + str(node.arrow[j].temp-len(set(cluster_labels))) + '\n'
      string = string + string_k
      k += 1

with open(data+'dbsnp_nz_inferred_tree.gv','w') as f:
  f.write('digraph G {\n')
  for node in tree_nw.traverse("preorder"):
    list_of_label = node.label
    delimiter = ','
    string_of_label = delimiter.join(list_of_label)
    f.write(str(node.temp-len(set(cluster_labels)))+' [label="'+string_of_label+'"];\n')
  f.write(string)
  f.write("}")
########################################################################################################################################33

for node in tree_nw.traverse("preorder"):
  if not node.label:
    node.delete(prevent_nondicotomic=False)



for node in tree_nw.traverse("preorder"): # arrow here refers to the arrow we use to tell the relationship between two nodes, so it basically has the 
                                          # information about the children of an internal node
    node.arrow = []


for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    children = node.get_children()
    node.arrow.extend(children)

evolution_pattern='Linear'
for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    children = node.get_children()
    if len(children) > 1:
       evolution_pattern = 'Branching'
       break
    
    


string = ''
for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    k=0
    for j in range(len(node.get_children())):
      string_k = str(node.temp) + ' -> ' + str(node.arrow[j].temp) + '\n'
      string = string + string_k
      k += 1


all_pos = []
for node in tree_nw.traverse("preorder"):
  all_pos.extend(node.label)

print("Number of labels ",len(all_pos))

with open(data+'dbsnp_nz_nonempty_inferred_tree.gv','w') as f:
  f.write('digraph G {\n')
  for node in tree_nw.traverse("preorder"):
    list_of_label = node.label
    delimiter = ','
    string_of_label = delimiter.join(list_of_label)
    f.write(str(node.temp)+' [label="'+string_of_label+'"];\n')
  f.write(string)
  f.write("}")

with open(data+'evolution_pattern.txt','w') as f:
  f.write(evolution_pattern)

s = Source.from_file(data+'dbsnp_nz_nonempty_inferred_tree.gv')
s.render()

