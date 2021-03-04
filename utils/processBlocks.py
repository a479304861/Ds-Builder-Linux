from utils import LCS as LCS
from pathlib import Path
import os

def readSequence(file):
    chr = []
    fr = open(file, 'r')
    for line in fr:
        line = line[:-2]
        itemset = line.split(' ')[1:]
        chr.append(itemset)
    chr_Dict = {}
    for i in range(len(chr)):
        chr_Dict[i] = chr[i]
    return chr_Dict

def syntenyDict(file):
    syntenyDict = {}
    with open(file) as sf:
        while True:
            line = sf.readline()[:-2]
            if not line:
                break
            itemset = line.split(' ')
            header = itemset[0].split(':')
            syntenyDict[header[0]] = itemset[1:]
    return syntenyDict

def assambleDrimmSequence(blockSequence,synteny):
    sequences = {}
    sequences_ID = {}
    blockCount = {}
    for i in blockSequence.keys():
        sequence = []
        sequence_ID = []
        for j in blockSequence[i]:
            if j.startswith('-'):
                block = j[1:]
                if block not in synteny.keys():
                    continue
                synteny_sequence = synteny[block][::-1]
                if block not in blockCount.keys():
                    blockCount[block] = 1
                else:
                    blockCount[block] += 1
                for k in range(len(synteny_sequence)):
                    sequence.append(synteny_sequence[k])
                    sequence_ID.append('-'+block+'|'+str(blockCount[block])+'|'+str(k))
            else:
                block = j
                if block not in synteny.keys():
                    continue
                synteny_sequence = synteny[block]
                if block not in blockCount.keys():
                    blockCount[block] = 1
                else:
                    blockCount[block] += 1
                for k in range(len(synteny_sequence)):
                    sequence.append(synteny_sequence[k])
                    sequence_ID.append('+' + block + '|' + str(blockCount[block]) + '|' + str(k))
        sequences[i] = sequence
        sequences_ID[i] = sequence_ID
    return sequences,sequences_ID

def readChrs(file):
    chr_list = []
    fr = open(file, 'r')
    for line in fr:
        line = line[:-1]
        chr_list.append(line)
    return chr_list

def processLCS(species_all_sequences,species_reassamble_sequences,species_all_sequences_name,species_reassamble_sequences_ID):
    block_range = {}
    for i in species_all_sequences.keys():
        species_reassamble_sequence = species_reassamble_sequences[i]
        print(i)
        block_range[i] = {}
        for j in species_all_sequences[i].keys():
            if j not in species_reassamble_sequence.keys():
                continue
            else:
                print(j)
                p = LCS.LCS()
                p.input(species_all_sequences[i][j]
                        ,species_reassamble_sequence[j])
                p.Compute_LCS()
                lcs = p.printOneLCS()
                for k in lcs:
                    genename = species_all_sequences_name[i][j][k[0]]
                    ID = species_reassamble_sequences_ID[i][j][k[1]]
                    ID_split = ID.split('|')
                    block = ID_split[0][1:]
                    block_count = ID_split[1]
                    if block not in block_range[i].keys():
                        block_range[i][block] = {}
                        block_range[i][block][block_count] = [[j,k[0],genename]]
                    else:
                        if block_count not in block_range[i][block].keys():
                            block_range[i][block][block_count] = [[j,k[0],genename]]
                        else:
                            block_range[i][block][block_count].append([j,k[0],genename])
    return block_range

def write(outfile,outfile_name,resultSynteny):
    outfile = open(outfile,'w')
    outfile_name = open(outfile_name,'w')
    for item in resultSynteny:
            outfile.write(str(item.block)+':'+str(item.species)+':'+str(item.block_count)+' ')
            outfile_name.write(str(item.block)+':'+str(item.species)+':'+str(item.block_count)+' ')
            for l in item.genename:
                    outfile_name.write(l+' ')
            for l in item.genesequence:
                    outfile.write(l+' ')     
            outfile.write('\n')
            outfile_name.write('\n')
    outfile.close()
    outfile_name.close()

def writeToFile(block_range,outdir,species_all_sequences_name,species_all_sequences):
    resultAllSynteny=[]
    for i in block_range.keys():
        resultSynteny=[]
        for j in block_range[i].keys():
            block = j
            for k in block_range[i][j].keys():
                block_count = k
                matching_pairs = block_range[i][j][k]
                matching_pairs = sorted(matching_pairs,key=lambda x:x[1])
                chr = matching_pairs[0][0]
                start = matching_pairs[0][1]
                end = matching_pairs[-1][1]
                genename = species_all_sequences_name[i][chr][start:end+1]
                genesequence = species_all_sequences[i][chr][start:end+1]
                resultSynteny.append(OutputBlock(block,block_count,genename,genesequence,i))
        resultSynteny=sorted(resultSynteny)
        resultAllSynteny+=resultSynteny
        outfile = outdir + i + '.synteny'
        outfile_name = outdir + i + '.synteny.genename'
        write(outfile,outfile_name,resultSynteny)
    resultAllSynteny=sorted(resultAllSynteny)
    outfile = outdir +   'all.synteny'
    outfile_name = outdir  + 'all.synteny.genename'
    write(outfile,outfile_name,resultAllSynteny)   

def processBlocks(workdir,speciesAndChrLen):
    # workdir = 'D:/InferAncestorGenome/realData/IAG_V2test/Brassica_ariticle/'
    block_dir = workdir +'/drimmSyntenyOutput/'
    synteny_dir = workdir +'/drimmSyntenyOutput/'
    sequences_dir = workdir + '/processOrthoFind/'
    species=[]    # species = ['Nigra','Oleracea','Rapa']
    for i in speciesAndChrLen:
        species.append(i)
    syntenyFile = synteny_dir + 'sourceSynteny.txt'
    outdir = workdir + '/block_synteny/'
    if(not Path(outdir).exists()):
        os.makedirs(outdir)
    blockSequences = {}
    for i in species:
        blockSequences[i]=readSequence(block_dir+i+'.block')
    synteny = syntenyDict(syntenyFile)
 
    
    species_reassamble_sequences = {}
    species_reassamble_sequences_ID = {}
    for i in blockSequences.keys():
        sequences,sequences_ID = assambleDrimmSequence(blockSequences[i],synteny)
        species_reassamble_sequences[i] = sequences
        species_reassamble_sequences_ID[i] = sequences_ID
   
    # # read allsequence and name
    species_all_sequences = {}
    species_all_sequences_name = {}
    for i in species:
        species_all_sequences[i] = readSequence(sequences_dir + i + '.all.sequence')
        species_all_sequences_name[i] = readSequence(sequences_dir + i + '.all.sequence.genename')
    # match sequence
    block_range= processLCS(species_all_sequences,species_reassamble_sequences,species_all_sequences_name,species_reassamble_sequences_ID)

    writeToFile(block_range,outdir,species_all_sequences_name,species_all_sequences)


class OutputBlock:
    def __init__(self,block,block_count,genename,genesequence,species):
        self.block=int(block)
        self.block_count=block_count
        self.genename=genename
        self.genesequence=genesequence
        self.species=species

    def __eq__(self,other):
        return self.block==other.block and self.block_count==other.block_count and self.species==other.species

    def __le__(self,other):
        if(self.block==other.block):
            if(self.species==other.species):
                return self.block_count<other.block_count
            return self.species<other.species
        return self.block<other.block 
    def __gt__(self,other):
        if(self.block==other.block):
            if(self.species==other.species):
                return self.block_count>other.block_count
            return self.species>other.species
        return self.block>other.block 





