import numpy as np 
import ete3 
from ete3 import Tree
import argparse
import random
import sys
import copy
import uuid 

def ThereAreDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

def gen_unique_ids(num):
	ids = []
	while True:
		ids = []
		for i in range(num):
			ids.append(str(uuid.uuid4()))
		if not ThereAreDuplicates(ids):
			break
	# ids = random.sample(range(100), num)
	return ids

# def destination_branches(tree, brnch_arr):
# 	branch_arr = []
# 	# print brnch_arr
# 	for node in tree.traverse("postorder"):
# 		if len(node.get_children())!=0:
# 			for child in node.get_children():
# 				if (node.name, child.name) not in brnch_arr:
# 					branch_arr.append((node, child))
# 		else:
# 			continue
# 	if len(branch_arr)==0:
# 		return None
# 	else:
# 		return branch_arr

# def list_branches(tree):
# 	branch_arr = []
# 	for node in tree.traverse("postorder"):
# 		if not node.is_root() and not node.is_leaf():
# 			for child in node.get_children():
# 				branch_arr.append((node, child))
# 		else:
# 			continue
# 	if len(branch_arr)==0:
# 		return None
# 	else:
# 		return branch_arr

def perform_SPR(in_tr, selection, dest):
	returned_trees = []
	selected_node = in_tr&selection
	# if selected_node.up.is_root():
	# 	subtree = selected_node.detach()
	# 	new_internals_names = []
	# 	for node in in_tr.traverse():
	# 		if not node.is_root() and not (node.up).is_root():
	# 			new_internals_names.append(node.name)
	# 	if len(new_internals_names)==0:
	# 		print("No valid SPR rearrangement for this branch")
	# 		returned_trees = None
	# 	else:
	# 		if dest>len(new_internals_names):
	# 			dest = len(new_internals_names)
	# 		samples = random.sample(range(len(new_internals_names)), k=dest)
	# 		for sample in samples:
	# 			tmp_in_tr = copy.deepcopy(in_tr)
	# 			tmp_subtree = copy.deepcopy(subtree)
	# 			s = tmp_in_tr&new_internals_names[sample]
	# 			parent = s.up
	# 			new_node = parent.add_child()
	# 			new_node.add_child(tmp_subtree)
	# 			s.detach()
	# 			new_node.add_child(s)
	# 			to_detach = tmp_in_tr.get_tree_root().get_children()[0]
	# 			returned_trees.append(to_detach.detach())
	# else:
	# flag = selected_node.up.up.is_root()
	to_remove = selected_node.up
	sister = selected_node.get_sisters()[0]
	subtree = selected_node.detach()
	to_remove.delete()
	new_internals_names = []
	# if flag:
	# for node in in_tr.traverse():
			# if not node.up.is_root() and not node.is_root():
				# new_internals_names.append(node.name)
	# else:
	for node in in_tr.traverse():
		if (not node.is_root()) and (node!=sister):
			new_internals_names.append(node.name)
	if len(new_internals_names)==0:
# 		print("No valid SPR rearrangement for this branch")
		returned_trees =  None
	else:
		if dest>len(new_internals_names):
			dest = len(new_internals_names)
		#samples = random.sample(range(len(new_internals_names)), dest)
		samples = range(len(new_internals_names))
		for sample in samples:
			tmp_in_tr = in_tr.copy()
			tmp_subtree = subtree.copy()
			s = tmp_in_tr&new_internals_names[sample]
			parent = s.up
			new_node = parent.add_child()
			new_node.add_child(tmp_subtree)
			s.detach()
			new_node.add_child(s)
			returned_trees.append(tmp_in_tr)

	# trees = []
	# (parent, selected_node) = selection 
	# # print selection
	# selected_node.detach()
	# tuple_ = [(parent.up.name, (parent.get_children()[0]).name)]
	# parent.delete()
	# options_arr = destination_branches(in_tr, tuple_)
	# # print options_arr
	# if options_arr == None:
	# 	return None
	# else:
	# 	''' we pick the other branch randomly '''
	# 	for k in random.sample(range(len(options_arr)), 1):
	# 		in_tr_cpy = copy.deepcopy(in_tr)
	# 		selected_node_cpy = copy.deepcopy(selected_node)
	# 		destination_arr = destination_branches(in_tr_cpy, tuple_)	
	# 		(p,c) = destination_arr[k]
	# 		# print (p,c)
	# 		c.detach()
	# 		in_node = p.add_child()
	# 		in_node.add_child(c)
	# 		in_node.add_child(selected_node_cpy)
	# 		# print in_tr_cpy

	# 		trees.append(in_tr_cpy)
	return returned_trees

def Main(in_tree, N, N_dest):
	''' Since it is not possible to deal with all the SPR rearrangements 
	given a particular topology, we randomly sample N new topologies in each call'''
	tree_list = []
	tr = in_tree
	# print in_tree
	nodes = []
	internal_names = gen_unique_ids(num=len(in_tree)-1)
	selected_names = []
	# lst_branches = list_branches(tr)
	# print len(lst_branches)
	# if len(lst_branches)==0:
	# 	# print("No valid SPR rearrangement for the given tree")
	# 	return tree_list
	if len(in_tree)<=3:
		print("No valid SPR rearrangement for the given tree")
		return tree_list
	else:
		counter = 0
		new_t = Tree(name="artificial_root")
		new_t.add_child(in_tree.copy())
		(new_t.get_children()[0]).name = internal_names[counter]
		counter+=1
		for node in new_t.traverse():
			if not node.is_root():
				if not (node.up).is_root():
					if node.is_leaf():
						nodes.append(node)
						selected_names.append(node.name)
					else:
						node.name = internal_names[counter]
						nodes.append(node)
						selected_names.append(node.name)
						counter+=1
		if N>len(nodes):
			N = len(nodes)
		# samples = random.sample(range(len(nodes)), k=N)
		samples = range(len(nodes))
		for i in samples:
			# top = new_t.get_topology_id()
			tmp_cpy = new_t.copy()
			tt = perform_SPR(in_tr=tmp_cpy, selection=selected_names[i], dest=N_dest)
			if tt!=None:
				# r = tt.get_tree_root().get_children()[0]
				# tree_list.extend(r.detach())
				# tree_list.extend([(x.get_tree_root().get_children()[0]).detach() for x in tt])
				for x in tt:
					if len(x.get_tree_root().get_children())==1:
						(x.get_tree_root().get_children()[0]).delete()
						tree_list.append(x)
					else:
						tree_list.append(x)

		# print(new_t)
		# for i in random.sample(range(len(lst_branches)), k=N):
		# 	tmp_cpy = copy.deepcopy(tr)
		# 	tmp_lst = list_branches(tmp_cpy)
		# 	# print tmp_lst[i]
		# 	# print "-------------------------------------------------------"
		# 	result_trs = perform_SPR(in_tr=tmp_cpy, selection=tmp_lst[i])
		# 	if result_trs == None:
		# 		pass
		# 		# print("This rearrangement resulted in the same tree")
		# 	else:
		# 		tree_list.extend(result_trs)
	# for Tr in tree_list:
		# print Tr
	# print len(tree_list)
	return tree_list

if __name__=="__main__":
	#### NOTE: each of the internal nodes must have a name for this rearrangement
	### Toy example
	### SPR must return 64 different topologies for this toy example
	t = Tree(name="root")
	Z = t.add_child(name="Z")
	Y = t.add_child(name="Y")
	# F = t.add_child(name="F")

	X = Y.add_child(name="X")
	V = Y.add_child(name="V")

	A = V.add_child(name="A")
	B = V.add_child(name="B")

	E = X.add_child(name="E")

	W = X.add_child(name="W")
	D = W.add_child(name="D")
	C = W.add_child(name="C")

	# A = Z.add_child(name="A")
	# B = Z.add_child(name="B")
	print(t)
	lst = Main(in_tree=t, N=100, N_dest=300)
	for tr_ in lst:
		print(tr_)
	print(len(lst))
