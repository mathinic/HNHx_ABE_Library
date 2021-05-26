# -*- coding: utf-8 -*-
"""
Analysis Script for HNHxABE Library. Sequence reads are initially demultiplexed by using protospacer,
barcode1 (5' of target sequence) and barcode2 (3' of target sequence).
Reads without a match or with ambiguous matches (due to recombinations during PCR or lentiviral packaging)
are discarded and only fully matched reads are further processed.
@author: mathinic
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import gzip
from os import listdir
import numpy as np
import matplotlib.pyplot as plt

def lookup(prototemplate, barcode1template, barcode2template):  # generate lookup dictionaries for all three relevant regions (protospacer, barcode1, barcode2)
    protolookup = {}
    for i, z in enumerate(prototemplate):
            protolookup.setdefault(z, []).append(i)
        
    barcode1lookup = {}
    for i, z in enumerate(barcode1template):
            barcode1lookup.setdefault(z, []).append(i)
            
    barcode2lookup = {}
    for i, z in enumerate(barcode2template):
            barcode2lookup.setdefault(z, []).append(i)
    
    return protolookup, barcode1lookup, barcode2lookup


def importfiles(filename, protolookup, barcode1lookup, barcode2lookup):
    '''
    Parameters
    ----------
    filename : string
        Fastq.gz file with sequencing reads.
    protolookup : dict
        Dictionary which is used as a lookup table for protospacer sequences in sequence reads.
    barcode1lookup : dict
        Dictionary which is used as a lookup table for barcode1 sequences in sequence reads.
    barcode2lookup : dict
        Dictionary which is used as a lookup table for barcode2 sequences in sequence reads.

    Returns
    -------
    diseasedict : dict
        Dictionary containing protospacer-,barcode1- and barcode2-matches for each sequencing read.

    '''
    diseasedict = {}
    shortname = filename[0:-16] #trim long name to short identifier

    filename = shortname+'_Protospacer.fastq.gz'
    with gzip.open(filename, "rt") as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            protospacer = seq_record.seq
            protomatch = protolookup.get(str(protospacer))
            identifier = seq_record.id
            if identifier in diseasedict:
                diseasedict[identifier]["protomatch"] = protomatch
            else:
                diseasedict[identifier] = {'protomatch':protomatch}
    print('Protospacer done')
    
    filename = shortname+'_barcode1.fastq.gz'
    with gzip.open(filename, "rt") as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            barcode1 = str(seq_record.seq)
            barcode1match = barcode1lookup.get(barcode1)
            identifier = seq_record.id
            # print(identifier)
            if identifier in diseasedict:
                diseasedict[identifier]["barcode1match"] = barcode1match
            else:
                diseasedict[identifier] = {'barcode1match':barcode1match}
    print('Barcode1 done')
    
    filename = shortname+'_target.fastq.gz'
    with gzip.open(filename, "rt") as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            barcode2 = str(seq_record.seq)[-6:]
            barcode2match = barcode2lookup.get(barcode2)
            identifier = seq_record.id
            if identifier in diseasedict:
                diseasedict[identifier]["barcode2match"] = barcode2match
            else:
                diseasedict[identifier] = {'barcode2match':barcode2match}
    print('barcode2 done')
    
    return diseasedict


templatedf = pd.read_csv('20210421_HNH_Library_Templatetable.csv') # import csv-file containing library sequence information

def list_files1(directory): # get list of all filenames with "target.fastq.gz" ending
    return [f for f in listdir(directory) if 'target.fastq.gz' in f]

### adapt path to the location of the fastq files!    
path = 'C:\\Users\\xxx\\xxx\\'
###

filelist = list_files1(path)
shortnamelist = []

averageeditingdict = {}
for filename in filelist: # perform analysis for all target.fastq.gz files in the directory (loop)
    shortname = filename[0:-16] # get name w/o "target.fastq.gz" ending
    prototemplate = templatedf[['protospacer']].values.flatten()
    barcode1template =templatedf[['barcode1']].values.flatten()
    barcode2template =templatedf[['barcode2']].values.flatten()
    protolookup, barcode1lookup, barcode2lookup = lookup(prototemplate,barcode1template,barcode2template) # create lookup tables
    diseasedict = importfiles(filename,protolookup, barcode1lookup,barcode2lookup)
    diseasedf = pd.DataFrame.from_dict(diseasedict,orient='index')  # make dataframe from dict
    diseasedf.dropna(inplace=True) # remove all reads which have missing matches
    
    # create column with identical match of protospacer and barcode1 (excluding recombinations)
    diseasedf['match'] = [list(set(a).intersection(set(b))) for a, b in zip(diseasedf.protomatch, diseasedf.barcode1match)]
    final_diseasedf = diseasedf[diseasedf['match'].map(lambda d: len(d)) == 1]  # discard variants with more than 1 match
    
    for x in range(20): # create empty columns which will be filled afterwards
        templatedf[x+1] = np.empty((len(templatedf), 0)).tolist()
        templatedf[str(x+1)+"_percent"] = np.empty((len(templatedf), 0)).tolist()
    
    for x in range(20):
        templatedf[str(x+1)+"_AtoG_editing"] = None
    templatedf["AtoG_editing_list"] = np.empty((len(templatedf), 0)).tolist()
    templatedf["A_positions"] = np.empty((len(templatedf), 0)).tolist()

    
    # Analyse fastqfiles by reading the reads into memory and evaluate editing based on matched library-variants in the importfiles function above.
    readsdict = {}
    targetbasedict = {}
    n = 4  # fastq file has one read for every 4 lines
    with gzip.open(filename, 'rt') as fastqfile:
        lines = []
        for line in fastqfile:
            lines.append(line.rstrip())
            if len(lines) == n:
                identifier = lines[0].split()[0][1:]
                if not identifier in final_diseasedf.index:  # skip reads which have no match (because of recombination or sequencing errors etc)
                    lines = []
                    continue
    
                variantindex = diseasedf.loc[identifier].match[0] # get variantindex based on previous read matching
                raw_identifier = lines[0]
                target = str(Seq(lines[1][-34:-14]).reverse_complement()) # Target sequence within the sequenced reads
                
                # add each base of read to list in dataframe:
                for loc, base in enumerate(target):
                    templatedf.loc[variantindex,loc+1].append(base)
                    
                targetbase = target[templatedf.loc[variantindex,'proto_position']-1]
                # add target to dictionary "readsdict"; all reads from the same target are hereby combined
                if variantindex in readsdict:
                    readsdict[variantindex].append(target)
                else:
                    readsdict[variantindex] = [target]
                    
                if variantindex in targetbasedict:
                    targetbasedict[variantindex].append(targetbase)
                else:
                    targetbasedict[variantindex] = [targetbase]
                lines = []
    
    for entry in readsdict: # analyse all the reads for each variant in the library
        
        count = Counter(readsdict[entry])
        count2 = Counter(targetbasedict[entry])
        nrunchanged = count[templatedf.loc[entry,'protospacer']]
        nrunchangedbase = count2[templatedf.loc[entry,'protospacer'][templatedf.loc[entry,'proto_position']-1]]
        templatedf.loc[entry,'percent_unedited'] = (nrunchanged/len(readsdict[entry]))*100
        templatedf.loc[entry,'percent_uneditedbase'] = (nrunchangedbase/len(targetbasedict[entry]))*100

        templatedf.loc[entry,'nrreads'] = len(readsdict[entry])
        for loc in range(20): # analyse each "A" in targetsequence/targetprotospacer
            templatedf.loc[entry,str(loc+1)+"_percent"] = [Counter(templatedf.loc[entry,loc+1])]
            if templatedf.loc[entry,'protospacer'][loc] == 'A':
                numberofA = templatedf.loc[entry,str(loc+1)+"_percent"]['G']
                totalnumberofreads = len(templatedf.loc[entry,loc+1])
                templatedf.loc[entry,str(loc+1)+"_AtoG_editing"] = numberofA/totalnumberofreads*100
                templatedf.loc[entry,"AtoG_editing_list"].append(numberofA/totalnumberofreads*100)
                templatedf.loc[entry,"A_positions"].append(loc+1)
            else:
                templatedf.loc[entry,str(loc+1)+"_AtoG_editing"] = None
        
        averageAtoG = []
        for baseposition in range(1,21):
            meanA = templatedf[str(baseposition)+"_AtoG_editing"].mean()
            averageAtoG.append(meanA)
    
    #create small plot for quick visualization of editing efficiencies at different positions
    plt.plot(list(range(1,21)), averageAtoG)
    plt.xticks(list(range(1,21)))
    plt.title("Editing window overview\nHNH library")
    plt.xlabel("Position in Protospacer")
    plt.ylabel("AtoG editing (%)")
    shortnamelist.append(shortname)
    plt.ylim([0, 40])
    averageeditingdict[filename] = {'overall_modified':100-templatedf['percent_unedited'].mean()}
    averageeditingdict[filename]['target_base_modified'] = 100-templatedf['percent_uneditedbase'].mean()
    
    #store final dataframe containing sequence analysis
    templatedf.to_csv('20210502_HNH_dataframe_'+shortname+'_Rep1.csv')

aveditdf = pd.DataFrame.from_dict(averageeditingdict,orient='index')
plt.legend(shortnamelist)
plt.savefig("Average_Editing_Variants.png")

    

