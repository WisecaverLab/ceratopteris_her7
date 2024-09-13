import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import Bio
from Bio import Seq, SeqIO
import ete3
from ete3 import Tree, faces, TreeStyle, PhyloTree, NodeStyle, TextFace, AttrFace, SeqMotifFace
import sys 

treefile = sys.argv[1] #'OG0009093.mpr.tree'
alnfile = sys.argv[2] #'OG0009093.aln'
speciesfile = sys.argv[3] #'species_ids.txt'
colorfile = sys.argv[4] #'lineage_colors.txt'
branch_min = float(sys.argv[5]) #95
outfile = sys.argv[6] #'OG0009093.mpr.pdf'

# Populate lookup dictionaries
seqDict = {}
alnLen = 0
for record in SeqIO.parse(alnfile, "fasta"):
    name = record.id
    sequence = str(record.seq)
    alnLen = len(sequence)
    seqDict[name] = sequence

code2name = {}
code2lineage = {}
fi = open(speciesfile)

for line in fi:
    if line[0] == '#':
        continue
    name, spid, lineage = line.rstrip().split('\t')
    #print(name, spid, lineage)
    code2name[spid] = name
    
    code2lineage[spid] = lineage

fi.close()

lin2color = {}
linorder = {}
fi = open(colorfile)

for line in fi:
    if line[0] == '#':
        continue

    order, color_desc, color, cat, lineage, taxids = line.rstrip().split('\t')  
    order = int(order)
    lin2color[lineage] = color
    linorder[order] = lineage

fi.close()


# add species and lineage information to all leaves
t = Tree(treefile, format=1)
t.ladderize(direction=1)
leafSet = set()

for n in t.get_leaves():
    leafSet.add(n.name)
    #print(n.name)
    tmp = n.name.split("_") 
    spid = tmp.pop(-1)
    tmp.pop(-1)
    speciesname = code2name[spid]
    genename = "_".join(tmp)
    #print(speciesname, genename, code2lineage[spid])
    
    n.add_features(lineage=code2lineage[spid])
    n.add_features(gene=genename)
    n.add_features(species=speciesname)

    #print("spid:", spid, "Species name:", n.species, "Species lineage:", n.lineage, "Color:", lin2color[n.lineage])

    # create a new label with a color attribute
    linF = AttrFace("lineage", fgcolor=lin2color[n.lineage], fsize=1)
    linF.background.color = lin2color[n.lineage]
    linF.margin_top = linF.margin_bottom = linF.margin_left = 10
    
    speciesF = AttrFace("species", fsize=10, fgcolor=lin2color[n.lineage], fstyle="italic")
    speciesF.margin_right = speciesF.margin_left = 10

    if spid == 'cri':
        geneF = AttrFace("gene", fsize=12, fgcolor="#e31a1c", fstyle="bold")
        geneF.margin_right = geneF.margin_left = 5
    
    elif spid == 'ath':
        geneF = AttrFace("gene", fsize=12, fgcolor="#e31a1c", fstyle="bold")
        geneF.margin_right = geneF.margin_left = 5
    
    else:
        geneF = AttrFace("gene", fsize=10, fgcolor="black")
        geneF.margin_right = geneF.margin_left = 5

    # labels aligned to the same level
    n.add_face(speciesF, 0, position='aligned')
    n.add_face(geneF, 0, position='branch-right')
    n.add_face(linF, 1, position='aligned')
    
    my_motifs = [[0, alnLen, "compactseq", 2, 10, None, None, None]]
    seqF = SeqMotifFace(seq=seqDict[n.name], motifs=my_motifs, gap_format="blank")
    seqF.margin_right = seqF.margin_left = 5
    n.add_face(seqF, 2, "aligned")
    

# add lineage information to all internal nodes
style = NodeStyle()

style["size"] = 0
style["hz_line_width"] = 2
style["vt_line_width"] = 2
t.set_style(style)

for n in t.iter_descendants("postorder"):
    #print(n.name)
                
    style["size"] = 0
    style["hz_line_width"] = 2
    style["vt_line_width"] = 2
    n.set_style(style)
    
    lineage_set = set()
    # get descendants, if all descendants are members of same lineage, color lineage color
    #print("NODE CHILDREN:")
    for k in n.iter_descendants("postorder"):
        for l in k.get_leaves():
            lineage_set.add(l.lineage)
            #print("Gene:", l.gene, "Species:", l.species, "Lineage:", l.lineage, "Color:", lin2color[l.lineage])
    
    #print(len(lineage_set), lineage_set)
    if len(lineage_set) == 1:
        node_lin = ''.join(lineage_set)
        #print(len(lineage_set), lineage_set, node_lin, lin2color[node_lin])
    
        newstyle = NodeStyle()
        newstyle["size"] = 0
        newstyle["hz_line_width"] = 2
        newstyle["vt_line_width"] = 2
        newstyle["vt_line_color"] = lin2color[node_lin]
        newstyle["hz_line_color"] = lin2color[node_lin]
        n.img_style = newstyle
        
    #fix branchlengths    
    if n.name not in leafSet and n.name[0] != 'n':
        #print(n.name)
        
        if float(n.name) >= branch_min:
            #branch_support = n.name.split('.')[0]
            branch_support = n.name
            #print(branch_support)
            n.add_features(bootstrap=branch_support)
            
            if len(lineage_set) == 1:
                node_lin = ''.join(lineage_set)
                supF = AttrFace("bootstrap", fgcolor=lin2color[node_lin], fsize=8)
                supF.margin_right = supF.margin_left = 3
                n.add_face(supF, 0, position='branch-bottom')
                
            else:
                supF = AttrFace("bootstrap", fgcolor="#000000", fsize=8)
                supF.margin_right = supF.margin_left = 3
                n.add_face(supF, 0, position='branch-bottom')
        
for n in t.get_leaves():
    
    leafstyle = NodeStyle()
    leafstyle["size"] = 0
    leafstyle["hz_line_width"] = 2
    leafstyle["vt_line_width"] = 2
    leafstyle["vt_line_color"] = lin2color[n.lineage]
    leafstyle["hz_line_color"] = lin2color[n.lineage]
    n.img_style = leafstyle
    
    
ts = TreeStyle()
ts.show_leaf_name = False
#ts.show_branch_support = True
ts.draw_guiding_lines = True

# add legend
ts.title.add_face(TextFace("Taxonomy:", fsize=10), column=0)
for i in range(1, len(lin2color)+1):
    #print(linorder[i], lin2color[linorder[i]])
    ts.title.add_face(TextFace(linorder[i], fsize=10, fgcolor=lin2color[linorder[i]]), column=0)


# render image on notebook or save to file
t.render(outfile, tree_style=ts)


