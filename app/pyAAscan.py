def PrimerTm(primer):
    '''
    Calculate the Tm for a given primer based on the GC content. 
    this function was missing in the source code, but based on the paper, this is set to be
    64.9°C + 41°C x (number of G’s and C’s in the primer – 16.4)/length of the primer
    The function was verified to produce the same output as the original AAscan application.  
    '''
    GC_content = len([i for i in primer if i in ['G', 'C', 'g', 'c']])
    skip_count = len([i for i in primer if i not in ['G', 'C', 'g', 'c', 'A', 'a', 'T', 't']])
    return round(64.9+ 41*(GC_content - 16.4)/(len(primer)-skip_count), 1)

def GC_content(primer):
    '''
    gives GC content as a ratio
    '''
    GC_content = len([i for i in primer if i in ['G', 'C', 'g', 'c']])
    return round(GC_content/len(primer), 2)

def TranslateCodon(codon):
    '''
    Translates a codon string into its matching amino acid (given in single letter code)
    '''
    c = codon.upper()
    result = '@'  #undefined codon - with n -
    if (c=='GCT') or (c=='GCC') or (c=='GCA') or (c=='GCG'): result='A'
    if (c=='CGT') or (c=='CGC') or (c=='CGA') or (c=='CGG') or (c=='AGA') or (c=='AGG'): result='R'
    if (c=='AAT') or (c=='AAC'): result='N'
    if (c=='GAT') or (c=='GAC'): result='D'
    if (c=='TGT') or (c=='TGC'): result='C'
    if (c=='CAA') or (c=='CAG'): result='Q'
    if (c=='GAA') or (c=='GAG'): result='E'
    if (c=='GGT') or (c=='GGC') or (c=='GGA') or (c=='GGG'): result='G'
    if (c=='CAT') or (c=='CAC') :result='H'
    if (c=='ATT') or (c=='ATC') or (c=='ATA'): result='I'
    if (c=='TTA') or (c=='TTG') or (c=='CTT') or (c=='CTC') or (c=='CTA') or (c=='CTG'): result='L'
    if (c=='AAA') or (c=='AAG'): result='K'
    if (c=='ATG'): result='M'
    if (c=='TTT') or (c=='TTC'): result='F'
    if (c=='CCT') or (c=='CCC') or (c=='CCA') or (c=='CCG'): result='P'
    if (c=='TCT') or (c=='TCC') or (c=='TCA') or (c=='TCG') or (c=='AGT') or (c=='AGC'): result='S'
    if (c=='ACT') or (c=='ACC') or (c=='ACA') or (c=='ACG'): result='T'
    if (c=='TGG'): result='W'
    if (c=='TAT') or (c=='TAC'): result='Y'
    if (c=='GTT') or (c=='GTC') or (c=='GTA') or (c=='GTG'): result='V'
    if (c=='TAA') or (c=='TGA'): result='*'
    if (c=='TAG'): result='$' #AMBER stop codon
    return result

def TranslateOligo(sequence, startpos=0, ignore_incomplete=True):
    '''
    Translates a given sequence in one go using the TranslateCodon function.  
    sequence : the input sequence
    startpos : the starting position from which to read the codons, by default set to 0
    ignore_ incomplete : ignore the last incomplete codon if present (default:True)    
    '''
    out = ''
    for i in range(startpos, len(sequence),3):
        codon = sequence[i:i+3]
        if len(codon) != 3 and ignore_incomplete:
            continue
        out+=str(TranslateCodon(codon))
    return out

def BestCodons(res):
    '''
    most likely two codons of each AA for e.coli, taken from the genscript website:
    https://www.genscript.com/tools/codon-frequency-table
    '''
    topcodons= {'R':['CGT','CGC'], 'H':['CAT','CAC'], 'K':['AAA','AAG'], 'D':['GAT','GAC'], 'E':['GAA','GAG'], 
                'S':['AGC','TCT'], 'T':['ACC','ACG'], 'N':['AAC','AAT'], 'Q':['CAG', 'CAA'], 
                'C':['TGC','TGT'], 'G':['GGC','GGT'], 'P':['CCG','CCA'], 
                'A':['GCG','GCC'], 'I':['ATT','ATC'], 'L':['CTG','TTA'], 'M':['ATG', 'ATG'], 
                'F':['TTT','TTC'], 'V':['GTG','GTT'], 'W':['TGG','TGG'], 'Y':['TAT','TAC']}
    if res in topcodons.keys():
        return topcodons[res]
    else:
        print('specified mutation is not in the canonical 20 amino acids: ', mutation)
        return

def CodonMatchScore(c1,c2):
    '''
    Given two codon strings, scores their similarity, 
    giving 1 point for matching A/T, and 2 points for matching G/C
    '''
    result = 0
    for i, j in zip(c1, c2):
        if i.upper() == j.upper():
            if i in ['g','G','c','C']:
                result+=2
            else:
                result+=1
    return result

def GCClampCount(primer, reverse=False):
    '''
    count how many residues from the end the first non-GC is found
    combination of the orignial GCClampScore and GCClampScoreRC2 functions,
    use reverse=False for GCClampScore and reverse=True for GCClampScoreRC2
    '''
    result=0
    # go through in reverse order if reverse is True
    rev = -1
    if reverse:
        rev = 1
    for i in primer[::rev]:
        if i in ['g','G','c','C']:
            result+=1
        else:
            return result
    return result

def GCClampScore(primer, reverse=False):
    '''
    scores CG claps, with 0 corresponding to the worst GC clamp 
    and 3 corresponding to the best, such that:
    [GC][GC][GC] = 0; 
    [ATGC][ATGC][AT] = 1 
    [ATGC][AT][GC] = 2 
    [AT][GC][GC] = 3, 
    '''
    result = -1
    clamp_count = GCClampCount(primer, reverse)
    if clamp_count >=3:
        result = 0
    if clamp_count == 0:
        result = 1
    if clamp_count == 1:
        result = 2
    if clamp_count == 2:
        result = 3
    return result

def ComplementNucleotide(ch):
    '''
    Match each character in the nucleotide sequence with its complement
    '''
    complement_dict = {'A':'T', 'a':'t', 'T':'A', 't':'a',
                       'G':'C', 'g':'c', 'C':'G', 'c':'g', 
                       'X':'Y', 'x':'y', 'Y':'X', 'y':'x', 
                       'N':'N', 'n':'n'}
    if ch in complement_dict.keys():
        return complement_dict[ch]
    else:
        return ch
    
def RCOligo(primer):
    '''
    creates reverse complement (RC) for a given primer
    '''
    result = ''
    for i in primer[::-1]:
        result+=ComplementNucleotide(i)
    return result

def PrimerLongFormat1(mutpos, strand, GCclamp, AnnLen, primer, fullprimer):
    '''
    Prints out all information about primers in format 'long1'
    '''
    result=str(mutpos)+'_'+strand+' '+primer+' len='+str(len(primer))+' Tm='+str(PrimerTm(primer))
    result+=' Tmfull='+str(PrimerTm(fullprimer))+' GCclamp='+str(GCclamp)+' AnnLenF='+str(AnnLen)
    result+=' GCcontent='+str(GC_content(fullprimer))+' '+fullprimer
    return result

def PrimerLongFormat2(mutpos, strand, GCclamp, AnnLen, primer, fullprimer):
    '''
    Prints out all information about primers in format 'long2'
    '''
    result=str(mutpos)+'_'+strand+' '+str(mutpos)+' '+strand+' '+primer+' '+str(len(primer))
    result+=str(PrimerTm(primer))+' '+str(PrimerTm(fullprimer))+' '+str(GCclamp)+str(AnnLen)
    result+=str(GC_content(fullprimer))+' '+fullprimer
    return result

def PrimerShortFormat(mutcodon, mutpos, newcodon, strand, fullprimer):
    '''
    Prints out all information about primers in format 'short'
    '''
    result=str(TranslateCodon(mutcodon))+str(mutpos)+str(TranslateCodon(newcodon))+'_'+strand+' '+fullprimer
    return result

def PrimermFastaFormat(mutpos, strand, fullprimer):
    '''
    Prints out all information about primers in format 'mFASTA'
    '''
    result='>'+str(mutpos)+'_'+strand+' '+fullprimer
    return result

def Mutate(seq_in, cod1pos, mutpos, codon1='GCG', codon2='GCA', codon3='GGT', codon4='GGC', 
           altset=['GCA','GCT','GCG','GCC'], minlen=18, maxlen=60, minGCcl=2, 
           mintm=60, maxtm=70, maxdtm=5, maxsugg=1, minanlen=15, minover=13, maxover=15,
           optGCcl=True, outputmode='long1', verbose=False):
    '''
    In the original pascal version, this function would take all input from the GUI, 
    and calculate a desired set of primers. Here the GUI and primer design are seperated
    
    seq_in  : DNA seq incl. 40-50 nt flanking regions
    cod1pos : Nt position of the 1st triplet
    mutpos  : position of mutation (as AA number)
    codon1  : codon to mutate to at desired position
    codon2  : alternative codon at desired position
    codon3  : in case the orignal codon is an in the altset, mutate to this
    codon4  : alternative codon at desired position for altset
    altset  : for these codons, mutate to codon3 or codon4 (default: if Ala, mutate to Gly)
    minlen  : minimum length of the final primer
    maxlen  : maximum length of the final primer
    minGCcl : minimum GCclamp score
    mintm   : minimum primer Tm
    maxtm   : maximum primer Tm
    maxsugg : maximum number of suggested primers
    minanlen: minimum annealing length
    minover : minimum overlap
    maxover : maximum overlap
    optGCcl : optimized GCclamp
    sepFR   : give forward and reverse primers as seperate outputs
    verbose : return additional information
    '''
    lines_out = []
    # prune to keep just GCAT bases
    seq = ''.join([i.lower() for i in seq_in if i in ['g','G','c','C','a','A','t','T']])
    # get position of mutation, remove -1 due to 0 vs 1 indexing between pascal and python
    ntmutpos = (cod1pos+3*(mutpos-1))-1
    mutcodon = seq[ntmutpos:ntmutpos+3]
    # if codon is in altset, switch to second set of codons
    if mutcodon[:2].upper() in altset:
        if CodonMatchScore(mutcodon, codon3) < CodonMatchScore(mutcodon, codon4):
            newcodon = codon4
        else:
            newcodon = codon3
    # if mutated colon is not in altset, use codon1 or codon2
    else:
        if CodonMatchScore(mutcodon, codon1) < CodonMatchScore(mutcodon, codon2):
            newcodon = codon2
        else:
            newcodon = codon1
    # create sequence with and without new mutated codon
    seqtemp1, seqtemp2 = seq, seq
    for i in range(0,3):
        if mutcodon[i].upper() == newcodon[i].upper():
            seqtemp1= seqtemp1[:ntmutpos+i]+newcodon[i]+seqtemp1[ntmutpos+i+1:]
        else:
            seqtemp1 = seqtemp1[:ntmutpos+i]+'X'+seqtemp1[ntmutpos+i+1:]
        seqtemp2= seqtemp2[:ntmutpos+i]+newcodon[i]+seqtemp2[ntmutpos+i+1:]
    # cycle through all possible primers and find the first pair that fits the criteria
    FW_found, RC_found = False,False
    for plength in range(minlen, maxlen+1):        
        for startpos in range(ntmutpos+3+minanlen-plength, ntmutpos+1):
            FW_found, RC_found = False,False
            # designing forward primer
            primer = seqtemp1[startpos:startpos+plength]
            fullprimer = seqtemp2[startpos:startpos+plength]
            annlen_f = startpos+plength-ntmutpos-3
            primertm = PrimerTm(primer)
            if optGCcl:
                GC = GCClampScore(primer)
            else:
                GC = GCClampCount(primer)
            # check if all conditions are met for the FW primer    
            if (primertm >= mintm) and (primertm <= maxtm) and (GC>=minGCcl):
                FW_found = True
            else:
                continue
            # designing reverse primer
            for plength_r in range(minlen, maxlen+1):
                for overlap in range(minover, maxover+1):
                    startpos_r = startpos-plength_r+overlap
                    primer_r = seqtemp1[startpos_r:startpos_r+plength_r]
                    fullprimer_r = seqtemp2[startpos_r:startpos_r+plength_r]
                    annlen_r = ntmutpos-startpos_r
                    primertm_r = PrimerTm(primer_r)
                    if optGCcl:
                        GC_r = GCClampScore(RCOligo(primer_r))
                    else:
                        GC_r = GCClampCount(RCOligo(primer_r))
                    # check if all conditions are met for the RC primer
                    if (primertm_r >= mintm) and (primertm_r <= maxtm) and (GC_r>=minGCcl):
                        RC_found = True
                        break
                if RC_found:
                    break
            break   
        if FW_found and RC_found:
            break
    # if no suitable primer is found, return *******
    if FW_found == False or RC_found == False:
        annlen_, annlen_r, GC, GC_r = 0,0,0,0
        primer, fullprimer = '*******', '*******'
        primer_r, fullprimer_r = '*******', '*******'
    
    # format results into proper output mode        
    if outputmode == 'long1':
        lines_out+=[PrimerLongFormat1(mutpos, 'F', GC, annlen_f, primer, fullprimer)]
        lines_out+=[PrimerLongFormat1(mutpos, 'R', GC_r, annlen_r, primer_r, RCOligo(fullprimer_r))]
    elif outputmode == 'long2':
        lines_out+=[PrimerLongFormat2(mutpos, 'F', GC, annlen_f, primer, fullprimer)]
        lines_out+=[PrimerLongFormat2(mutpos, 'R', GC_r, annlen_r, primer_r, RCOligo(fullprimer_r))]
    elif outputmode == 'short':
        lines_out+=[PrimerShortFormat(mutcodon, mutpos, newcodon, 'F', fullprimer)]
        lines_out+=[PrimerShortFormat(mutcodon, mutpos, newcodon, 'R', RCOligo(fullprimer_r))]
    elif outputmode == 'mFASTA':
        lines_out+=[PrimermFastaFormat(mutpos, 'F', fullprimer)]
        lines_out+=[PrimermFastaFormat(mutpos, 'R', RCOligo(fullprimer_r))]  
    else:
        print('ERROR: nonexistent output mode specified: {}'.format(outputmode))
        print('please use any of the following: long1, long2, short, mFASTA')
        return
    # if verbose is set as True, specify additional information
    if verbose:
        lines_verbose=['sequence length {} nt'.format(len(seq))]
        lines_verbose+=['mutating position {}, codon {} ({})'.format(mutpos, mutcodon, TranslateCodon(mutcodon))]
        if mutcodon[:2].upper() == 'GC':
            lines_verbose+=['match score c1: {} ({})'.format(CodonMatchScore(mutcodon, codon3), codon3)]
            lines_verbose+=['match score c2: {} ({})'.format(CodonMatchScore(mutcodon, codon4), codon4)]
        else:
            lines_verbose+=['match score c1: {} ({})'.format(CodonMatchScore(mutcodon, codon1), codon1)]
            lines_verbose+=['match score c2: {} ({})'.format(CodonMatchScore(mutcodon, codon2), codon2)]
        lines_verbose+=['using {} ({}) to replace {} ({})'.format(newcodon, TranslateCodon(newcodon), 
                                                                  mutcodon, TranslateCodon(mutcodon))]
        lines_verbose+=['mutation: {}{}{}'.format(TranslateCodon(mutcodon), mutpos, TranslateCodon(newcodon))]
        lines_out = lines_verbose+lines_out
    return lines_out

def Parse_Seq(seqfile):
    '''
    Parses a file containing a sequence and returns just the sequence. 
    it saves just the ATCG content (without spaces) from each line if:
    - the line contains only ATCG (covers plain sequence, FASTA)
    - the line contains only ATCG or numbers/spaces (covers EMBL, GCG, GB formats)
    '''
    with open(seqfile, 'r') as seqf:
        seq = seqf.readlines()
    seq = [i[:-1] for i in seq]
    seq_out = ''
    for line in seq:
        line_s = line.replace(" ", "")
        if set(line_s).issubset(set('AaTtCcGg')):
            seq_out+=line_s
        elif set(line_s).issubset(set('AaTtCcGg0123456789')):
            line_s = ''.join([i for i in line_s if i in set('AaTtCcGg')])
            seq_out+=line_s
        else:
            continue
    return seq_out

if __name__ == "__main__":
    import argparse
    import sys
    parser = argparse.ArgumentParser(prog='pyAAscan',
                    description='a tool for making primers based on the original AAscan implemented in pascal')
    parser.add_argument('-cod', '--firstcodonpos', required=True, type=int)
    parser.add_argument('-seq', '--sequencefile', required=True)      
    parser.add_argument('-mutf', '--mutfile', default='')
    parser.add_argument('-mut', '--mutation', default='')
    parser.add_argument('-minl', '--minlen', default=18)
    parser.add_argument('-maxl', '--maxlen', default=60)
    parser.add_argument('-mincg', '--minGCclamp', default=2)
    parser.add_argument('-mintm', '--mintm', default=60)
    parser.add_argument('-maxtm', '--maxtm', default=70)
    parser.add_argument('-maxdtm', '--maxdtm', default=5)
    parser.add_argument('-maxsug', '--maxsuggested', default=1)
    parser.add_argument('-minanl', '--minanneallen', default=15)
    parser.add_argument('-mino', '--minoverlap', default=13)
    parser.add_argument('-maxo', '--maxoverlap', default=15)
    parser.add_argument('-aascan', '--aascan', action='store_true', default=False)
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('-o', '--out', default='AAscan.out')
    parser.add_argument('-f', '--format', default='short')
    args = parser.parse_args()       
    # if a single mutation is specified, use that
    if args.mutation != '':
        muts = [args.mutation]
    # if a file with mutations is specified, use that
    elif args.mutfile != '':
        with open(args.mutfile, 'r') as mutf:
            muts = mutf.readlines()
        muts = [i[:-1] for i in muts]
    # if no mutations are given, exit
    else:
        print('specifiy either a single mutation with -mut or a file of mutations with -mutf')
        sys.exit()
    # get sequence from file and check for start codon
    sequence = Parse_Seq(args.sequencefile)
    if TranslateCodon(sequence[args.firstcodonpos-1:args.firstcodonpos+2]) != 'M':
        print('WARNING: the specified first codon does not translate to M')
    lines_out = []
    # if AAscan mode, go over each mutation and mutate to A or G
    if args.aascan:
        for mut in muts:
            mutpos = int(''.join([i for i in mut if i in '0123456789']))
            lines_out+= Mutate(cod1pos=args.firstcodonpos, seq_in=sequence,  mutpos=mutpos, 
                               minlen=args.minlen, maxlen=args.maxlen, mintm=args.mintm, maxtm=args.maxtm, 
                               maxdtm=args.maxdtm, minGCcl=args.minGCclamp,minover=args.minoverlap, 
                               maxover=args.maxoverlap, outputmode=args.format, verbose=args.verbose)
    # if not AAscan mode, go over each mutation and mutate to specified mutation
    else:    
        for mut in muts:
            mutpos = int(''.join([i for i in mut if i in '0123456789']))
            codon1, codon2 = BestCodons(mut[-1])
            lines_out+= Mutate(cod1pos=args.firstcodonpos, seq_in=sequence,  mutpos=mutpos, 
                               codon1=codon1, codon2=codon2, altset=[], minlen=args.minlen, maxlen=args.maxlen, 
                               mintm=args.mintm, maxtm=args.maxtm, maxdtm=args.maxdtm, minGCcl=args.minGCclamp,
                               minover=args.minoverlap, maxover=args.maxoverlap, outputmode=args.format, verbose=args.verbose)
    # print and save result
    if args.verbose:
        for line in lines_out:
            print(line) 
    with open(args.out, 'w') as out:
        out.writelines([i+'\n' for i in lines_out])
    

    
    
    
