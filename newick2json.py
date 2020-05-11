#!/usr/bin/env python3

from sys import argv
import argparse
import UPhO


outg=""
ing=""

class node:
    """Class related to UPhO.split but a single partition, listing the descencedants of that node difines it"""
    def __init__(self):
        self.children=None
        self.branch_length=None
        self.support=None
        self.name=None
        self.parent=None
        self.level=None
        self.size=None
        self.droot=0



def main():
    usage = u"""
 \t\t\u26F1 \u001b[31;1m newick2json.py \u001b[0m \u26F1
\n 
    Convert your newick file in hierafchichal json format.\nThis script is written for python3 and requires the corresponding UPhO (https://github.com/ballesterus/UPhO)for split decoposition.
    
    newick2json.py <NEWICKTREE> <OUTGROUP,TAXA,SEPARATED,WITH,COMMAS>
    """
    print (usage)
    F=open(argv[1], "r")
    Oname =argv[1].split('.')[0] + ".json"
    Fc=F.readline()
    T = UPhO.myPhylo(Fc)
    F.close()
    print("File was read")
    global outg
    global ing
    myoutg=argv[2].split(",")

    #    myoutg=['Epiperipatus_sp', 'Opisthopatus_cinctipes', 'Peripatopsis_capensis']
    outg=set(myoutg)
    ing=set(T.leaves) - outg


    R=clados(T)
    update_parents(R)
    L=ladderize(R)
        
    with open(Oname, "w") as out:
        out.write(build_json(L,R))
    print("Enjoy editing you json. Bye bye")
    

    
    
def belongs_to_og(split_vec):
    s = set(split_vec)
    if outg.issuperset(s) and not ing.issuperset(s):
        return True
    else:
        return False

def belongs_to_in(split_vec):
    s = set(split_vec)
    if ing.issuperset(s) and not outg.issuperset(s):
        return True
    else:
        return False

def clados(phylo):
    """Returns a dictionary of nodes in a UPhO.myPhylo object to reflect polarity based on an outgroup (rooting).
    It also assigns node names and updates the brachlengts leading to the root"""
    counter=0
    result = {}
    rnode=node()
    rnode.name="root"
    rnode.children=phylo.leaves
    rnode.size=len(rnode.children)
    rnode.level=0
    rnode.droot=0
    result[rnode.name]=rnode
    for s in phylo.splits:
        for v in s.vecs:
            if len(v)==1:
                tnode=node()
                tnode.name=v[0]
                tnode.size = 1
                tnode.branch_length=s.branch_length
                tnode.support=None
                tnode.children=v
                result[tnode.name]=tnode
            elif set(v) == outg:
                nbl = float(s.branch_length)/2
                inode=node()                
                inode.children = v
                inode.size=len(v)
                inode.name="outgroup"
                inode.parent="root"
                inode.branch_length=nbl
                inode.support = s.support
                result[inode.name]=inode
                inode=node()                
                inode.children  = list(ing)
                inode.name="ingroup"
                inode.size=len(inode.children)
                inode.parent="root"
                inode.branch_length=nbl
                inode.support = s.support
                result[inode.name]=inode
            elif belongs_to_og(v):
                inode=node()
                inode.children=v
                inode.size=len(v)
                inode.name= "n_" + str(counter)
                counter += 1
                inode.branch_length=s.branch_length
                inode.support=s.support
                result[inode.name]=inode
            elif belongs_to_in(v) and len(v) < len(ing):
                inode=node()
                inode.children=v
                inode.size=len(v)
                inode.name= "n_" + str(counter)
                counter += 1
                inode.branch_length=s.branch_length
                inode.support=s.support
                result[inode.name]=inode
    return result
                
                
def find_mommy(nodeName, nodes_dict):
    """Updates the parent node in a collection of nodes(dictionary)"""
    q=nodes_dict[nodeName]
    qc=set(q.children)
    parent="root"
    min_psize=q.size
    c_psize = nodes_dict["root"].size #This as big as a parent can be a we want the smallest
    for n in nodes_dict.keys():
        cn = set(nodes_dict[n].children)
        if len(cn) > min_psize and len(cn) < c_psize:
            if qc.issubset(cn):
                parent=nodes_dict[n].name
                c_psize=len(cn)
    q.parent=parent

def find_children(nodeName,node_dict):
    """To run after all parents have been identified"""
    result=[]
    asize=node_dict['root'].size
    for k in node_dict.keys():
        if node_dict[k].parent == nodeName:
            if node_dict[k].size < asize:
                result.insert(0, k)
                asize=node_dict[k].size
            else:
                result.append(k)
    return result

def update_parents(nodes_dict):
    for k in nodes_dict:
        if nodes_dict[k].name != "root":
            find_mommy(k,nodes_dict)
            
            
def ladderize(nodes_dict):
    """Updates levels and Return and list of node keys ordered descendingly"""
    ladder=["root"]
    queue=[]
    init=find_children('root', nodes_dict)
    queue= queue + init
    while len(queue) > 0:
        for q in queue:
            queue.remove(q)
            ladder.append(q)
            queue=queue + find_children(q, nodes_dict)
    for n in ladder:
        cp=nodes_dict[n].parent
        if cp != None:
            nodes_dict[n].level = nodes_dict[cp].level + 1
            nodes_dict[n].droot= nodes_dict[cp].droot + float(nodes_dict[n].droot)
    return ladder
            

def json_node(node_name, nodes_dict):
    jstring=""
    node= nodes_dict[node_name]
    desc=find_children(node_name, nodes_dict)
    children_placeholder="#" +'#,#'.join(desc) + "#"
    pad="\t" * node.level
    if node.name == "root":
        jstring="%s{\n%s\"name\" : \"%s\",\n%s\"children\" : [\n%s\n%s]\n%s}\n" %(pad,pad,node.name,pad, children_placeholder, pad,pad)
    else:
        if node.size > 1:
            jstring="\n%s{\n%s\"name\" : \"%s\",\n%s\"parent\" : \"%s\",\n%s\"support\" : %s,\n%s\"branch-length\" : %s,\n%s\"children\" : [\n%s\n%s]\n%s}\n" %(pad,pad,node.name,pad, node.parent, pad, node.support,pad, node.branch_length, pad, children_placeholder, pad,pad)
        else:
            jstring="\n%s{\n%s\"name\" : \"%s\",\n%s\"parent\" : \"%s\",\n%s\"branch-length\" : %s\n%s}\n" %(pad,pad,node.name,pad, node.parent,pad, node.branch_length, pad)
    return jstring


def build_json(ladder, nodes_dict):
    construct=json_node(ladder[0], nodes_dict)  
    for n in ladder:
        search = "#" + n + "#" # makesure
        replac = json_node(n, nodes_dict)
        temp = construct.replace(search,replac)            
        construct = temp
    return construct

if __name__ == "__main__":
    main()
    
