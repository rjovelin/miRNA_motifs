def make_mirna_set_with_without_motif(mast_output, mirna_upstream):
    '''
    (file, file) -> dict
    Return a diictionnary of mirnas with and without motif by comparing the mirnas with motifs
    in mast_output and the entire set of mirnas in mirna_upstream    
    '''

    # make a set of mirna that have a motif
    from get_positions_from_mast import motif_position_from_mast
    mast_positions = motif_position_from_mast(mast_output)
    mast_motif = mast_positions[1]
    mirna_with_motif = set()
    for mirna in mast_motif:
        mirna_with_motif.add(mirna)



    # make a set of mirna that have no motif
    from mirna_data_for_meme import convert_fasta
    upstreams = convert_fasta(mirna_upstream)
    mirna_no_motif = set()
    for mirna in upstreams:
        if not mirna in mirna_with_motif:
            mirna_no_motif.add(mirna)

    # make a dictionnary to store information about whether mirna has motif or not
    mirna_motif_status = {}
    for mirna in mirna_with_motif:
        mirna_motif_status[mirna] = 'has_motif'
    for mirna in mirna_no_motif:
        mirna_motif_status[mirna] = 'no_motif'

    return mirna_motif_status



def expression_mirna_withmotif_nomotif(mast_output, mirna_expression, mirna_upstream, outputfile):
    '''
    (file, file, file, file) -> file
    Parition mirna expression into 2 categories for mirnas in mirna_upstream that have or not a motif in mast_output 
    '''

    # make a dictionnary to record the motif status
    mirna_motif_status = make_mirna_set_with_without_motif(mast_output, mirna_upstream)
    
    # make a dictionnary storing mirna expression {mirna:expression}
    expression_file = open(mirna_expression, 'r')
    expression = {}
    line = expression_file.readline()
    for line in expression_file:
        line = line.rstrip()
        if line != '':
            line = line.split()
            mirna = line[0]
            expression[mirna] = line[1]
            
    expression_file.close()

    # write the mirna expression values for mirna with and without motif into the output file
    myfile = open(outputfile, 'w')
    myfile.write('mirna' + '\t' + 'has_motif' + '\t' + 'expression' + '\n')
     
    for mirna in expression:
        for mirna_name in mirna_motif_status:
            if mirna == mirna_name[0: mirna_name.index('_')]:
                myfile.write(mirna_name + '\t' + mirna_motif_status[mirna_name] + '\t' + expression[mirna] + '\n')
               
    myfile.close()



def expression_mirna_no_motif(mast_output1, mast_output2, mast_output3, mast_output4, mast_output5, mirna_expression, mirna_upstream, outputfile):
    '''
    (file, file, file, file, file, file, file) -> file
    Save in outputfile the expression level of miRNAs lacking all motifs m1 to m5
    '''

    # make dictionnaries to record the motif status of each motif
    mirna_motif1 = make_mirna_set_with_without_motif(mast_output1, mirna_upstream)
    mirna_motif2 = make_mirna_set_with_without_motif(mast_output2, mirna_upstream)
    mirna_motif3 = make_mirna_set_with_without_motif(mast_output3, mirna_upstream)
    mirna_motif4 = make_mirna_set_with_without_motif(mast_output4, mirna_upstream)
    mirna_motif5 = make_mirna_set_with_without_motif(mast_output5, mirna_upstream)

    # make a set of mirnas with no motif
    no_motif = set()

    # make a set of mirnas with any motif
    any_motif = set()

    # get the mirnas with no motif and mirna with motifs
    for mirna in mirna_motif1:
        if mirna_motif1[mirna] == 'has_motif':
            any_motif.add(mirna)
        elif mirna_motif1[mirna] == 'no_motif':
            no_motif.add(mirna)

    for mirna in mirna_motif2:
        if mirna_motif2[mirna] == 'has_motif':
            any_motif.add(mirna)
        elif mirna_motif2[mirna] == 'no_motif':
            no_motif.add(mirna)

    for mirna in mirna_motif3:
        if mirna_motif1[mirna] == 'has_motif':
            any_motif.add(mirna)
        elif mirna_motif3[mirna] == 'no_motif':
            no_motif.add(mirna)

    for mirna in mirna_motif4:
        if mirna_motif1[mirna] == 'has_motif':
            any_motif.add(mirna)
        elif mirna_motif4[mirna] == 'no_motif':
            no_motif.add(mirna)

    for mirna in mirna_motif5:
        if mirna_motif1[mirna] == 'has_motif':
            any_motif.add(mirna)
        elif mirna_motif5[mirna] == 'no_motif':
            no_motif.add(mirna)

    # remove any mirnas from no_motif set that are also in the any_motif set
    for mirna in any_motif:
        if mirna in no_motif:
            no_motif.discard(mirna)
    
    # make a dictionnary storing mirna expression {mirna:expression}
    expression_file = open(mirna_expression, 'r')
    expression = {}
    line = expression_file.readline()
    for line in expression_file:
        line = line.rstrip()
        if line != '':
            line = line.split()
            mirna = line[0]
            expression[mirna] = line[1]

    # write the mirna expression values for mirna in the no_motif set
    if len(no_motif) != 0:
        myfile = open(outputfile, 'w')
        myfile.write('mirna' + '\t' + 'has_motif' + '\t' + 'expression' + '\n')

        for mirna in expression:
            for mirna_name in no_motif:
                if mirna == mirna_name[0: mirna_name.index('_')]:
                    myfile.write(mirna_name + '\t' + 'no_motif' + '\t' + expression[mirna] + '\n')
        myfile.close()
    else:
        print('no_motif is empty')

def expression_mirna_two_motifs(mast_output1, mast_output2, mirna_expression, mirna_upstream, outputfile):
    '''
    (file, file, file, file, file) -> file
    Partition mirna expression into 4 categories for mirnas in mirna_upstream that have:
    motif1 + motif2, motif1_alone, motif2_alone, no_motif
    '''

    # make a dictionnary to record the motif status
    mirna_motif1 = make_mirna_set_with_without_motif(mast_output1, mirna_upstream)
    mirna_motif2 = make_mirna_set_with_without_motif(mast_output2, mirna_upstream)

    # make sets of mirnas that have or not the motifs
    motif1_motif2 = set()
    motif1 = set()
    motif2 = set()
    no_motif = set()

    for mirna in mirna_motif1:
        if mirna_motif1[mirna] == 'has_motif' and mirna_motif2[mirna] == 'has_motif':
            motif1_motif2.add(mirna)
        elif mirna_motif1[mirna] == 'has_motif' and mirna_motif2[mirna] == 'no_motif':
            motif1.add(mirna)
        elif mirna_motif1[mirna] == 'no_motif' and mirna_motif2[mirna] == 'has_motif':
            motif2.add(mirna)
        elif mirna_motif1[mirna] == 'no_motif' and mirna_motif2[mirna] == 'no_motif':
            no_motif.add(mirna)
      
    
    # make a dictionnary storing mirna expression {mirna:expression}
    expression_file = open(mirna_expression, 'r')
    expression = {}
    line = expression_file.readline()
    for line in expression_file:
        line = line.rstrip()
        if line != '':
            line = line.split()
            mirna = line[0]
            expression[mirna] = line[1]
            
    expression_file.close()

    # write the mirna expression values for mirna with and without motif into the output file
    myfile = open(outputfile, 'w')
    myfile.write('mirna' + '\t' + 'motif_status' + '\t' + 'expression' + '\n')

    # the mirna name in expression is a substring of the mirna name in the motif sets
    for mirna in expression:
        for mirna_name in motif1_motif2:
            if mirna == mirna_name[0: mirna_name.index('_')]:
                myfile.write(mirna_name + '\t' + 'both_motif' + '\t' + expression[mirna] + '\n')

    for mirna in expression:
        for mirna_name in motif1:
            if mirna == mirna_name[0: mirna_name.index('_')]:
                myfile.write(mirna_name + '\t' + 'alone_1' + '\t' + expression[mirna] + '\n')

    for mirna in expression:
        for mirna_name in motif2:
            if mirna == mirna_name[0: mirna_name.index('_')]:
                myfile.write(mirna_name + '\t' + 'alone_2' + '\t' + expression[mirna] + '\n')

    for mirna in expression:
        for mirna_name in no_motif:
            if mirna == mirna_name[0: mirna_name.index('_')]:
                myfile.write(mirna_name + '\t' + 'no_motif' + '\t' + expression[mirna] + '\n')

              
    myfile.close()


def count_motifs_expression(mast_output, mirna_expression, mirna_upstream, outputfile):
    '''
    (file, file, file, file) -> file
    Save in the outputfile the number of occurence of a motif found upstream of a given mirna and the expression level of that mirna stored in the expression file.
    '''

    # get the number of motif occurence
    from get_positions_from_mast import motif_position_from_mast
    mast_positions = motif_position_from_mast(mast_output)
    mast_motif = mast_positions[1]

    # get the all the mirnas into a dictionary
    from mirna_data_for_meme import convert_fasta
    upstreams = convert_fasta(mirna_upstream)

    
    # make a dictionnary storing mirna expression {mirna:expression}
    expression_file = open(mirna_expression, 'r')
    expression = {}
    line = expression_file.readline()
    for line in expression_file:
        line = line.rstrip()
        if line != '':
            line = line.split()
            mirna = line[0]
            expression[mirna] = line[1]
            
    expression_file.close()

    # write the mirna expression values for mirna with and without motif into the output file
    myfile = open(outputfile, 'w')
    myfile.write('mirna' + '\t' + 'N_motif' + '\t' + 'expression' + '\n')

    for mirna in upstreams:
        if mirna[0:mirna.index('_')] in expression:
            if not mirna in mast_motif:
                myfile.write(mirna + '\t' + '0' + '\t' + expression[mirna[0:mirna.index('_')]] + '\n')
            else:
                myfile.write(mirna + '\t' + str(len(mast_motif[mirna])) + '\t' + expression[mirna[0:mirna.index('_')]] + '\n')
                  
    myfile.close()



def get_mirna_targets(mirna_motif_status):
    '''
    (dict) -> dict
    Take the dictionnary of mirnas with motif status and return a dictionnary including the set of target genes for mirnas with and without motif
    '''


    all_targets = ['mir50_targets.txt', 'mir244_targets.txt', 'mir235_targets.txt', 'mir1_targets.txt', 'mir4816_targets.txt', 'mir245_targets.txt', 'mir79_targets.txt',
                   'mir1019_targets.txt', 'mir2_targets.txt', 'mir71_targets.txt', 'mir2211_targets.txt', 'mir1828_targets.txt', 'mir1817_targets.txt', 'mir795_targets.txt',
                   'mir794_targets.txt', 'mir1823_targets.txt', 'mir1824_targets.txt', 'mir2218a_targets.txt', 'mir2218b_targets.txt', 'mir1818_targets.txt', 'mir4805_targets.txt',
                   'mir72_targets.txt', 'mir1830_targets.txt', 'lin4_targets.txt', 'mir1822_targets.txt', 'mir60_targets.txt', 'mir236_targets.txt', 'mir57_targets.txt', 'mir85_targets.txt',
                   'mir2216_targets.txt', 'mir252_targets.txt', 'mir35_targets.txt', 'mir36_targets.txt', 'mir37_targets.txt', 'mir38_targets.txt', 'mir39_targets.txt', 'mir40_targets.txt',
                   'mir41_targets.txt', 'mir355_targets.txt', 'mir45_targets.txt', 'mir42_targets.txt', 'mir43_targets.txt', 'mir44_targets.txt', 'mir4806_targets.txt', 'mir77_targets.txt',
                   'mir2215_targets.txt', 'mir234_targets.txt', 'mir229_targets.txt', 'mir64_targets.txt', 'mir65_targets.txt', 'mir66_targets.txt', 'mir76_targets.txt', 'mir4813_targets.txt',
                   'mir67_targets.txt', 'mir2213_targets.txt', 'mir231_targets.txt', 'mir356a_targets.txt', 'mir356b_targets.txt', 'mir80_targets.txt', 'mir238_targets.txt', 'mir90_targets.txt',
                   'mir1020_targets.txt', 'mir1832a_targets.txt', 'mir86_targets.txt', 'mir46_targets.txt', 'mir4814_targets.txt', 'mir2209b_targets.txt', 'mir2208a_targets.txt', 'mir2208b_targets.txt',
                   'mir2209a_targets.txt', 'mir2209c_targets.txt', 'mir58_targets.txt', 'mir242_targets.txt', 'mir243_targets.txt', 'mir2217_targets.txt', 'mir228_targets.txt', 'mir790_targets.txt',
                   'mir4815_targets.txt', 'mir58b_targets.txt', 'mir1820_targets.txt', 'mir83_targets.txt', 'mir2210_targets.txt', 'mir232_targets.txt', 'mir51_targets.txt', 'mir53_targets.txt',
                   'mir59_targets.txt', 'mir124_targets.txt', 'mir52_targets.txt', 'mir1833_targets.txt', 'mir789_1_targets.txt', 'mir1832b_targets.txt', 'mir789_2_targets.txt', 'mir792_targets.txt',
                   'mir1821_targets.txt', 'mir255_targets.txt', 'mir253_targets.txt', 'mir70_targets.txt', 'mir357_targets.txt', 'mir358_targets.txt', 'mir259_targets.txt', 'lsy6_targets.txt',
                   'mir250_targets.txt', 'mir61_targets.txt', 'mir87_targets.txt', 'mir48_targets.txt', 'mir241_targets.txt', 'mir800_targets.txt', 'mir2221_targets.txt', 'mir796_targets.txt',
                   'mir4807_targets.txt', 'mir4808_targets.txt', 'mir4809_targets.txt', 'mir2220_targets.txt', 'mir4810_targets.txt', 'mir1018_targets.txt', 'mir248_targets.txt', 'mir73_targets.txt',
                   'mir74_targets.txt', 'mir75_targets.txt', 'mir81_targets.txt', 'mir82_targets.txt', 'mir34_targets.txt', 'mir359_targets.txt', 'mir785_targets.txt', 'mir249_targets.txt',
                   'mir247_targets.txt', 'mir797_targets.txt', 'mir230_targets.txt', 'mir791_targets.txt', 'mir360_targets.txt', 'mir1022_targets.txt', 'mir4811_targets.txt', 'mir240_targets.txt',
                   'mir786_targets.txt', 'mir784_targets.txt', 'mir237_targets.txt', 'mir788_targets.txt', 'mir799_targets.txt', 'mir49_targets.txt', 'mir251_targets.txt', 'mir787_targets.txt',
                   'mir254_targets.txt', 'mir239b_targets.txt', 'mir239a_targets.txt', 'mir233_targets.txt', 'mir62_targets.txt', 'mir56_targets.txt', 'mir55_targets.txt', 'mir54_targets.txt',
                   'mir392_targets.txt', 'mir793_targets.txt', 'mir47_targets.txt', 'let7_targets.txt', 'mir1829b_targets.txt', 'mir1829c_targets.txt', 'mir4812_targets.txt', 'mir1829a_targets.txt',
                   'mir84_targets.txt', 'mir2212_targets.txt', 'mir1819_targets.txt', 'mir63_targets.txt']




    # make a dictionnary with keys being the motif status and values being the set of filenames containing the mirna targets 
    mirna_motif_filenames = {}

    for mirna_name in mirna_motif_status:
        mirna = mirna_name[4:mirna_name.index('_')]
        mirna = mirna[0:3] + mirna[4:]
        mirna = mirna.replace('-', '_')
        mirna = mirna + '_targets.txt'
        if mirna in all_targets:
            if mirna_motif_status[mirna_name] in mirna_motif_filenames:
                mirna_motif_filenames[mirna_motif_status[mirna_name]].add(mirna)
            else:
                mirna_motif_filenames[mirna_motif_status[mirna_name]] = {mirna}

    # combine targets of mirnas with and without motif in 2 separate sets

    targets_mirna_with_motif = set()
    for mirna_filename in mirna_motif_filenames['has_motif']:
        myfile = open(mirna_filename, 'r')
        for line in myfile:
            line = line.strip()
            if line != '':
                targets_mirna_with_motif.add(line)
        myfile.close()
    targets_mirna_no_motif = set()
    for mirna_filename in mirna_motif_filenames['no_motif']:
        myfile = open(mirna_filename, 'r')
        for line in myfile:
            line = line.strip()
            if line != '':
                targets_mirna_no_motif.add(line)
        myfile.close()

    # make a dictionnary with motif status has key and set of targets as value
    mirna_targets = {}
    mirna_targets['has_motif'] = targets_mirna_with_motif
    mirna_targets['no_motif'] = targets_mirna_no_motif

    return mirna_targets

    

def mirna_targets_printer(mirna_targets, motif_mirna_targets, no_motif_mirna_targets):
    '''
    (dict, str) -> file
    Save the mirna targets in dictionnary mirna_targets into file motif_mirna_targets if the mirnas
    have the motif and into file no_motif_mirna_targets if the mirna don't have a motif
    '''

    motif_file = open(motif_mirna_targets, 'w')
    for target in mirna_targets['has_motif']:
        motif_file.write(target + '\n')
    motif_file.close()

    no_motif_file = open(no_motif_mirna_targets, 'w')
    for target in mirna_targets['no_motif']:
        no_motif_file.write(target + '\n')
    no_motif_file.close()

    
        
        

    

    

    
    
    
