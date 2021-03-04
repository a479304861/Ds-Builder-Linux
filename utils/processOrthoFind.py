import pandas as pd
import numpy as np
import random
import os.path
from pathlib import Path

# dirInput = 'G:/pythonCode/Ds_extents/program/input/'
# dirOuput = 'G:/pythonCode/Ds_extents/program/output/'
def getSequence(bed,group_dir):
    bed = pd.read_csv(bed,sep='\t',header=None)[[0,1]]
    chrlist = bed[0].unique().tolist()
    new_sequence = []
    new_name_sequence = []
    for i in chrlist:
        split = bed.loc[bed[0] == i][1].tolist()
        chr = []
        chr_name = []
        for j in split:
            # if group_dir.has_key(j):
            if j in group_dir.keys():
               chr.append(group_dir[j])
               chr_name.append(j)
        new_sequence.append(chr)
        new_name_sequence.append(chr_name)
    return new_sequence,new_name_sequence

def outSequence(sequence,outfile,sample_rate):
    outfile = open(outfile,'w')
    for i in sequence:
        for j in i:
            ran = random.random()
            if ran > sample_rate:
                outfile.write(str(j)+' ')
        outfile.write('\n')
    outfile.close()

def outAllSequence(sequence,outfile):
    outfile = open(outfile,'w')
    for i in sequence:
        for j in i:
            outfile.write(str(j)+' ')
        outfile.write('\n')
    outfile.close()

def processOrthoFind(orthogroupsPath,gff_list,dirOuput,size,simpleRage):
    # gff_list = ['PR_genome.gff','PS_genome.gff','PT_genome.gff']
    # sp = ['PR','PS','PT']
    dirOuput+='/processOrthoFind'
    if(not Path(dirOuput).exists()):
        os.makedirs(dirOuput)
    ortho = pd.read_csv(orthogroupsPath,sep='\t')
    ortho = ortho.fillna('')
    columns = ortho.columns.tolist()
    sp=columns[1:]
    ortho = np.asarray(ortho)
    group_dir = {}
    for i in ortho:
        group = i[0]
        species = i[1:]
        group_dir[group] = {}
        for j in range(len(species)):
            genes = species[j].split(', ')
            if genes[0] == '':
                group_dir[group][sp[j]] = []
            else:
                group_dir[group][sp[j]] = genes
    rate_dir = {}
    finalGroup = {}
    for i in group_dir.keys():
        rate_list = []
        for j in group_dir[i].keys():
            rate_list.append(len(group_dir[i][j]))
        ok = 1
        for j in rate_list:
            if j == 0:
                ok = 0
            if j > size:
                ok = 0
        if ok == 0:
            continue
        else:
            rate = ''
            for j in rate_list:
                rate += str(j) + ':'
            rate = rate[:-1]
            finalGroup[i] = group_dir[i]
            if rate not in rate_dir.keys():
                rate_dir[rate] = 1
            else:
                rate_dir[rate] += 1
    # for i in rate_dir.keys():
    #     print(i + '\t' + str(rate_dir[i]))
    outfile = dirOuput + '/group.xls'
    count = 1
    outfile = open(outfile,'w')
    outfile.write('gene\tgroup\n')
    for i in finalGroup.keys():
        for j in finalGroup[i].keys():
            for k in finalGroup[i][j]:
                outfile.write(k+'\t'+str(count)+'\n')
        count += 1
    outfile.close()
    group = pd.read_csv(dirOuput +'/group.xls',sep='\t')
    group = np.asarray(group)
    group_dir = {}
    for i in group:
        group_dir[i[0]] = i[1]
    # sample 50%
    allSequence=[]
    allsequence_name=[]
    speciesAndChrLen={}
    for i in gff_list:
        gff = i
        spliti= i.split('/')
        gff_name=spliti[len(spliti)-1].split('.')[0]
        sequence, sequence_name = getSequence(gff,group_dir)
        speciesAndChrLen[gff_name]=len(sequence)
        allSequence+=sequence
        allsequence_name+=sequence_name
        outallfile = dirOuput +"/"+ gff_name + '.all.sequence'
        outallfilename = dirOuput +"/"+ gff_name + '.all.sequence.genename'

        outAllSequence(sequence, outallfile)
        outAllSequence(sequence_name, outallfilename)
    outSequence(allSequence, dirOuput +"/" + 'sample.sequence', simpleRage) 
    outAllSequence(allSequence, dirOuput +"/" + 'all.sequence') 
    outAllSequence(allsequence_name, dirOuput +"/" + 'all.sequence.genename') 
    return speciesAndChrLen


