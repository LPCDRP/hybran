#!/usr/bin/env python
from Bio import SeqIO                           # import needed modules from BioPython
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import Reference
import operator                                 # imports operator for sorting feature objects
import csv                                   	# imports csv module, to read the csv file into an iterator
import os        	                            # imports os to remove temporary files in python scrypt
import argparse                                 # imports arparse so ursers or makefiles can specify in/out filenames


def tssaddotate(reference_genome_embl, csvfile, tsscol, strandcol,  out_put,
                csvauthors='Teresa Cortes, Olga T. Schubert, Graham Rose, Kristine B. Arnvig, '
                           + 'Inaki aki Comas, Ruedi Aebersold, and Douglas B. Young',
                ptitle='Genome-wideMapping of Transcriptional Start Sites Defines an Extensive Leaderless '
                       + 'Transcriptome in Mycobacterium tuberculosis',
                pjournal='Cell Reports', pubmed='24268774',
                exp='COORDINATES: TSS from 5\' RNAseq',
                sameref=0, flipstrand=False, note_columns='', usrnote='', dontsort=0, moreref=0):
    # defines function that takes in as arguments:
    # reference_genome_embl defines the emble file to be the starting reference genome
    # csvfile defines the source csv file to be added to reference genome
    # tsscol defines the colum with TSS genome position. each column is its own TSS, and will become it's own feature.
        # call this function seperately for each columm
    # strandcol defines column with the TSS strand data.
    # out_put defines name of output file
    # csvauthors defines the authors of the source
    # ptitle defines the title of the paper the TSS data came from
    # pjournal defines the name of the journal that published the paper the TSS data came from
    # pubmed defines pubmedid for csv source to be added (the journal article it came from)
    # exp defines string for experiment qualifier. Defaults to
        # 'COORDINATES: TSS from 5\' RNAseq'    \' should be read into string as '
    # sameref set this to True when calling with data from the same author as last call. Default is 0
    # flipstrand tells the function to use the opposite strand from that indicated in the Strand Column for each TSS
    # note_columns defines which column from the csv will have their values added as note qualifiers
        # Default no column is used
    # usrnote provides a string to be added to the feature under the note qualifier if present.
        # Default nothing is added
    # dontsort tells the function not to sort the features (just sort in the last call).
        # Default will sort.
    # moreref adds some additional citations. default doesn't add them.

    # reads in the emble file. makes one SeqRecord object with many SeqFeatures for each gene/cds/etc.
    ref_genome = SeqIO.read(reference_genome_embl, "embl")

    # add citation
    if moreref != 0:
        moresource = Reference()
        moresource2 = Reference()
        moresource.authors = 'Lew J., Ecole Polytechnique Federale de Lausanne, CH-1015, Lausanne,'
        moresource2.authors = 'Lew J.M.'
        moresource2.pubmed_id = '1-4411532'
        ref_genome.annotations['references'] = [moresource, moresource2]

    source_num = 0       # sets a counter for the citation qualifier
    if 'references' in ref_genome.annotations:
        # moves the citation counter past sources already present in the reference genome
        source_num = source_num + len(ref_genome.annotations['references'])
    if sameref == 0:									   # adds citation if this argument is default
        my_source = Reference()                            # creates new Reference object
        my_source.authors = csvauthors
        my_source.pubmed_id = pubmed
        my_source.title = ptitle
        my_source.journal = pjournal
        if 'references' in ref_genome.annotations:
            # adds the new reference to the reference class. in embl file its the RN's etc. at the top
            ref_genome.annotations['references'].append(my_source)
        else:
            ref_genome.annotations['references'] = my_source
        # gives the index of the reference, to be used in the citation qualifier for the new features.
        source_num = source_num + 1

    # add TSS data as new features
    with open(csvfile) as the_file:
        reader = csv.DictReader(the_file, delimiter='\t')
        for row in reader:       		                            # each row in the csv file should be a new TSS
                if row[tsscol] != '':						    	# avoids blank entries in csv

                    # determine TSS's genome location
                    tssstart = int(row[tsscol]) - 1                 # -1 for python counting
                    if str(row[strandcol]) == '+':
                        tssstrand = 1
                    else:
                        tssstrand = -1

                    # creates the new feature to store the TSS
                    tp1 = ExactPosition(tssstart)
                    tp2 = ExactPosition(tssstart + 1)
                    my_feature_location = FeatureLocation(tp1, tp2, tssstrand)
                    my_feature = SeqFeature(my_feature_location, type='mRNA')

                    # adds some qualifiers
                    my_feature.qualifiers['locus_tag'] = row['RvNumber']

                    # provides the integer of the corrosponding source in the Reference class
                    my_feature.qualifiers['citation'] = source_num

                    # provides experiment qualifier from a function argument
                    my_feature.qualifiers['experiment'] = exp

                    # add misc data to note qualifier
                    note_holder = 'Transcription Start Site.'
                    if usrnote != '':
                        note_holder = note_holder + ' ' + usrnote
                    if note_columns != '':  # adds any misc csv column data specified in function argument
                        note_holder = note_holder + ' ' + note_columns + ': ' + str(row[note_columns]) + '.'
                    my_feature.qualifiers['note'] = note_holder  # creates a note qualifier for the new feature

                    # adds the feature to the recod object
                    ref_genome.features.append(my_feature)

    # sorts features by start location in genome, unless told not to, also fixes translation attribute of CDSs
    if dontsort == 0:
        for cds in ref_genome.features:
            if cds.type == 'CDS' and 'pseudogene' not in cds.qualifiers:
                genesequence = cds.extract(ref_genome.seq)
                if len(genesequence)%3 != 0:
                    print(cds.qualifiers['locus_tag'])
                    print(len(genesequence))
                cds.qualifiers['translation'] = genesequence.translate(table="Bacterial")
        ref_genome.features.sort(key=operator.attrgetter("location.start"))

    # write out new file containing added annotations
    SeqIO.write(ref_genome, out_put, 'embl')


# calls TSSaddotate with the CSV data
def addotatcaller(inemble='ittas-man-h37.embl', csvpath='csvdata/', outemble='itass-fincatotat.embl'):
    #bodgy debug of thing i don't understand:
    csvpath = csvpath.replace('\t', '')
    csvpath = csvpath.replace(' ', '')

    # call 1
    print(csvpath + 'mod_cortez13_grow_internal_tss.csv')
    tssaddotate(reference_genome_embl=inemble,
                csvfile=(csvpath + 'mod_cortez13_grow_internal_tss.csv'),
                out_put='catotat1.embl',
                sameref=0, pubmed='24268774',
                csvauthors='Teresa Cortes, Olga T. Schubert, Graham Rose, Kristine B. Arnvig, '
                           + 'Inaki aki Comas, Ruedi Aebersold, and Douglas B. Young',
                ptitle='Genome-wideMapping of Transcriptional Start Sites Defines an Extensive Leaderless '
                       + 'Transcriptome in Mycobacterium tuberculosis',
                pjournal='Cell Reports',
                tsscol='Internal TSS 1 Genome position',
                strandcol='Internal TSS Strand under Exponential Growth',
                note_columns='Internal TSS 1 Expression During Growth', dontsort=1, moreref=1)

    # calls 2-13
    for filedex in range(2, 14):
        in_file = 'catotat' + str(filedex - 1) + '.embl'
        out_file = 'catotat' + str(filedex) + '.embl'
        t_col = 'Internal TSS ' + str(filedex) + ' Genome position'
        n_col = 'Internal TSS ' + str(filedex) + ' Expression During Growth'
        tssaddotate(reference_genome_embl=in_file, csvfile=(csvpath + 'mod_cortez13_grow_internal_tss.csv'),
                    out_put=out_file, sameref=1, tsscol=t_col,
                    strandcol='Internal TSS Strand under Exponential Growth',
                    note_columns=n_col, dontsort=1)

    # deleate temporary files:
    for filedex in range(1, 13):
        dead_file = 'catotat' + str(filedex) + '.embl'
        os.remove(dead_file)

    # calls 14-24
    for filedex in range(1, 12):
        in_file = 'catotat' + str(filedex + 12) + '.embl'
        out_file = 'catotat' + str(filedex + 13) + '.embl'
        t_col = 'Internal TSS ' + str(filedex) + ' Genome position'
        n_col = 'Internal TSS ' + str(filedex) + ' Expression During Starvation'
        tssaddotate(reference_genome_embl=in_file, csvfile=(csvpath + 'mod_cortez13_starv_internal_tss.csv'),
                    out_put=out_file, sameref=1, tsscol=t_col,
                    strandcol='Internal TSS Strand under Starvation',
                    note_columns=n_col, dontsort=1)

    # deleate more temporary files:
    for filedex in range(13, 24):
        dead_file = 'catotat' + str(filedex) + '.embl'
        os.remove(dead_file)

    # call 25
    tssaddotate(reference_genome_embl='catotat24.embl', csvfile=(csvpath + 'mod_cortez13_starv_RPKM_tss.csv'),
                out_put='catotat25.embl', sameref=1,
                tsscol='Starvation Expressed Primary TSS Genome Position',
                strandcol='Starvation Expressed TSS Strand', usrnote='Expressed During Starvation',
                dontsort=1)
    # call 26
    tssaddotate(reference_genome_embl='catotat25.embl', csvfile=(csvpath + 'mod_cortez13_starv_RPKM_tss.csv'),
                out_put='catotat26.embl', sameref=1,
                tsscol='Starvation Expressed Alternate TSS Genome Position',
                strandcol='Starvation Expressed TSS Strand', usrnote='Expressed During Starvation',
                dontsort=1)

    # call 27
    tssaddotate(reference_genome_embl='catotat26.embl', csvfile=(csvpath + 'mod_cortez13_starv_RPKM_tss.csv'),
                out_put='catotat27.embl', sameref=1,
                tsscol='Starvation Expressed 2nd Alternate TSS Genome Position',
                strandcol='Starvation Expressed TSS Strand', usrnote='Expressed During Starvation',
                dontsort=1)

    # call 28
    tssaddotate(reference_genome_embl='catotat27.embl', csvfile=(csvpath + 'mod_cortez13_starv_RPKM_tss.csv'),
                out_put='catotat28.embl', sameref=1,
                tsscol='Starvation Expressed 3rd Alternate TSS Genome Position',
                strandcol='Starvation Expressed TSS Strand', usrnote='Expressed During Starvation',
                dontsort=1)

    # call 29
    tssaddotate(reference_genome_embl='catotat28.embl', csvfile=(csvpath + 'mod_cortez13_starv_RPKM_tss.csv'),
                out_put='catotat29.embl', sameref=1,
                tsscol='Starvation Expressed 4th Alternate TSS Genome Position',
                strandcol='Starvation Expressed TSS Strand', usrnote='Expressed During Starvation',
                dontsort=1)

    # deleate more temporary files:
    for filedex in range(24, 29):
        dead_file = 'catotat' + str(filedex) + '.embl'
        os.remove(dead_file)

    # call 30
    tssaddotate(reference_genome_embl='catotat29.embl',
                csvfile=(csvpath + 'mod_cortez13_starv_RPKM_tss_antisense.csv'),
                out_put='catotat30.embl', sameref=1,
                tsscol='Starvation Expressed Strongest Peak Height Antisense TSS Genome Position',
                strandcol='Starvation Expressed Antisense TSS Strand', flipstrand=False,
                usrnote='Expressed During Starvation', dontsort=1)

    # call 31
    tssaddotate(reference_genome_embl='catotat30.embl',
                csvfile=(csvpath + 'mod_cortez13_starv_RPKM_tss_antisense.csv'),
                out_put='catotat31.embl', sameref=1,
                tsscol='Starvation Expressed 2nd Antisense TSS Genome Position',
                strandcol='Starvation Expressed Antisense TSS Strand', flipstrand=False,
                usrnote='Expressed During Starvation', dontsort=1)

    # call 32
    tssaddotate(reference_genome_embl='catotat31.embl',
                csvfile=(csvpath + 'mod_cortez13_starv_RPKM_tss_antisense.csv'),
                out_put='catotat32.embl', sameref=1,
                tsscol='Starvation Expressed 3rd Antisense TSS Genome Position',
                strandcol='Starvation Expressed Antisense TSS Strand', flipstrand=False,
                usrnote='Expressed During Starvation', dontsort=1)

    # calls 33-35
    for filedex in range(1, 4):
        in_file = 'catotat' + str(filedex + 31) + '.embl'
        out_file = 'catotat' + str(filedex + 32) + '.embl'
        t_col = 'Starvation Expressed ' + str(filedex + 3) + 'th Antisense TSS Genome Position'
        tssaddotate(reference_genome_embl=in_file,
                    csvfile=(csvpath + 'mod_cortez13_starv_RPKM_tss_antisense.csv'),
                    out_put=out_file, sameref=1, tsscol=t_col,
                    strandcol='Starvation Expressed Antisense TSS Strand',
                    flipstrand=False, usrnote='Expressed During Starvation', dontsort=1)

    # this last call uses a different source paper:
    # Pubmed: 26536359
    # Authors: Scarlet S. Shell1, Jing Wang, Pascal Lapierre, Mushtaq Mir, Michael R. Chase, Margaret M. Pyle,
    # Richa Gawande, Rushdy Ahmad, David A. Sarracino, Thomas R. Ioerger, Sarah M. Fortune, Keith M. Derbyshire,
    # Joseph T. Wade, Todd A. Gray

    # call 36
    tssaddotate(reference_genome_embl='catotat35.embl', csvfile=csvpath + 'mod_shell15_TSS_all.csv',
                out_put=outemble, pubmed='26536359',
                csvauthors='Scarlet S. Shell1, Jing Wang, Pascal Lapierre, Mushtaq Mir, Michael R. Chase, '
                + 'Margaret M. Pyle, Richa Gawande, Rushdy Ahmad, David A. Sarracino, Thomas R. Ioerger, '
                + 'Sarah M. Fortune, Keith M. Derbyshire, Joseph T. Wade, Todd A. Gray',
                ptitle='Leaderless Transcripts and Small Proteins Are Common Features of the Mycobacterial'
                + 'Translational Landscape',
                pjournal='PLoS Genetics',
                tsscol='TSS Coordinate (genome version NC_000962)', strandcol='TSS Strand')

    # deleate the last temporary files:
    for filedex in range(29, 36):
        dead_file = 'catotat' + str(filedex) + '.embl'
        os.remove(dead_file)


# function to merge duplicate TSS mRNA features
# defined as TSS's on the same strand, within the duprange parameter (default 2 bp) of bp from each other
def featmerge(refembl='fincatotat.embl', outembl='v3TSSmannotat.embl', duprange=2):
    # reads in the emble file.
    ref_genome = SeqIO.read(refembl, "embl")

    listmergedfeatures = []
    for feat in ref_genome.features:
        if feat.type == 'mRNA' and str(feat.qualifiers.get('note', 0)) != 'featuremergepurge' \
                and 'Transcription Start Site' in str(feat.qualifiers.get('note', 0)):
            mergedfeature = SeqFeature(feat.location, type='mRNA')
            mfloctag = feat.qualifiers['locus_tag']
            mfcit = []
            mfexp = ''
            mfnote = ''
            for ofeat in ref_genome.features:   # this will count itself, but that's ok.
                                                # every unique TSS will get "merged" with itself
                if ofeat.type == 'mRNA' and 'Transcription Start Site' in str(ofeat.qualifiers.get('note', 0)) \
                        and (feat.location.start + duprange) > ofeat.location.start > (feat.location.start - duprange) \
                        and (feat.location.end + duprange) > ofeat.location.end > (feat.location.end - duprange) \
                        and ofeat.location.strand == feat.location.strand:
                    # store data for merging
                    ofcitlist = ofeat.qualifiers.get('citation', 0)
                    if ofcitlist != 0:
                        mfcit.append(str(ofcitlist[0]))

                    ofexplist = ofeat.qualifiers.get('experiment', 0)
                    if ofexplist != 0:
                        mfexp = mfexp + str(ofexplist[0]) + ', '

                    ofnotlist = ofeat.qualifiers.get('note', 0)
                    # also removing 'Transcription Start Site.' so it doesn't get added multiple times to
                    # same merged feature. code will add this back to each merged feature at the end.
                    if ofnotlist != 0:
                        ofnote = str(ofnotlist[0]).replace('Transcription Start Site.', '', 1)
                    else:
                        ofnote = ''
                    if ofnote != '':
                        mfnote = mfnote + ofnote + ', '

                    # mark for purge (and tells meta loop not to check this one)
                    ofeat.qualifiers['note'] = 'featuremergepurge'

            mergedfeature.qualifiers['locus_tag'] = mfloctag
            mergedfeature.qualifiers['citation'] = mfcit
            mergedfeature.qualifiers['experiment'] = 'COORDINATES: TSS from 5\' RNAseq'
            mergedfeature.qualifiers['note'] = 'Transcription Start Site. ' + mfnote

            # store merged feature object for later inclusion in SeqRecord
            listmergedfeatures.append(mergedfeature)

    # remove features marked for purging (aka all TSS mRNA features) from SeqRecord
    ref_genome.features = [item for item in ref_genome.features
                           if str(item.qualifiers.get('note', 0)) != 'featuremergepurge']

    # add merged features to SeqRecord
    for mfeat in listmergedfeatures:
        ref_genome.features.append(mfeat)

    # sort features and write out
    ref_genome.features.sort(key=operator.attrgetter("location.start"))
    SeqIO.write(ref_genome, outembl, 'embl')


def main():
    # take in user/makefile's arguments as filenames for function calls
    parser = argparse.ArgumentParser(description='Adds TSS data to embl file',
                                     epilog='Transcription Start Site data comes from specific csv files')
    parser.add_argument('-r', '--reference', help='Input Reference embl File', default='ittas-man-h37.embl')
    parser.add_argument('-o', '--output', help='Output embl file', default='ittas-man-h37-tss-merge.embl')
    parser.add_argument('-c', '--csvdir', help='Directory with csv files with TSS data', default='csvdata/')
    args = parser.parse_args()

    print(args.reference)
    print(args.output)
    print(args.csvdir)

    tempfile = 'ittas-fincatotat.embl'  # created by addotatcaller and used by featmerge

    # call addotatcaller (which calls tssaddotate multiple times for each column of TSSs in csv files)
    addotatcaller(inemble=args.reference, csvpath=args.csvdir, outemble=tempfile)

    # calls featmerge to combine duplicate TSS data into same feature
    featmerge(refembl=tempfile, outembl=args.output, duprange=2)

    # deleate tempfile
    os.remove(tempfile)


if __name__ == '__main__':
    main()
