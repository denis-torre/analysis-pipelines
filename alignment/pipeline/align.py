import os, sys
import pandas as pd

def align(readFiles, outfile, method='STAR', genome='Homo_sapiens.GRCh38', threads=5, paired=False, kallisto_length=200, kallisto_sd=30):

    # Get base dir
    basedir = os.path.realpath(__file__).replace('/pipeline/align.py', '')

    # STAR
    if method == 'STAR':

        # Get genome dir
        genomeDir = '{basedir}/s2-star.dir/{genome}/genomeDir'.format(**locals())
        

        # Command
        command = ''' STAR \
            --genomeDir {genomeDir}  \
            --outFileNamePrefix {outfile}  \
            --readFilesIn {readFiles}  \
            --readFilesCommand gzcat \
            --quantMode GeneCounts \
            --limitBAMsortRAM 10000000000  \
            --limitIObufferSize 50000000 \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical  \
            --outSAMtype BAM SortedByCoordinate  \
            --outReadsUnmapped Fastx \
            --runThreadN {threads}
        '''.format(**locals())

    # kallisto
    elif method == 'kallisto':

        # Get genome index
        genomeIndex = '{basedir}/s1-kallisto.dir/{genome}.cdna.all.idx'.format(**locals())

        # Single or paired
        single = '--single' if paired == False else ''

        # Statement
        command = ''' kallisto quant \
            {single} \
            -t {threads} \
            -l {kallisto_length} \
            -s {kallisto_sd} \
            -i {genomeIndex} \
            -o {outfile} \
            {readFiles}
        '''.format(**locals())

    # other
    else:
        raise ValueError('Please indicate STAR or kallisto as the alignment method.')

    # Run
    os.system(command)

def merge(countFiles, method='STAR', collapse=True, annotation='GRCh38.p12', filter_biotype=False):

    # Get base dir
    basedir = os.path.realpath(__file__).replace('/pipeline/align.py', '')
    annotation_dataframe = pd.read_table('{basedir}/ensembl.dir/biomart/{annotation}-genes.txt'.format(**locals()))

    # Filter biotype
    if filter_biotype:
        annotation_dataframe = annotation_dataframe[annotation_dataframe['Gene type'].isin(filter_biotype)]

    # Initialize results
    results = []

    # STAR
    if method == 'STAR':

        # Read infiles
        for countFile in countFiles:

            # Read data
            dataframe = pd.read_table(countFile, header=None, names=['gene_symbol', 'counts_unstranded', 'counts_1', 'counts_2'])

            # Get sample name
            dataframe['sample'] = countFile.split('/')[-2]

            # Append
            results.append(dataframe)

        # Pivot
        count_dataframe = pd.concat(results).pivot(index='gene_symbol', columns='sample', values='counts_unstranded').drop(['N_ambiguous', 'N_multimapping', 'N_noFeature', 'N_unmapped'])

        # Collapse
        if collapse:
            count_dataframe = count_dataframe.merge(annotation_dataframe[['Gene stable ID', 'Gene name']].drop_duplicates('Gene stable ID'), left_index=True, right_on='Gene stable ID').drop('Gene stable ID', axis=1).groupby('Gene name').sum().rename_axis('gene_symbol')


    # kallisto
    elif method == 'kallisto':

        # Read infiles
        for countFile in countFiles:

            # Read data
            dataframe = pd.read_table(countFile)

            # Get sample name
            dataframe['sample'] = countFile.split('/')[-2]

            # Append
            results.append(dataframe)

        # Pivot
        count_dataframe = pd.concat(results).pivot(index='target_id', columns='sample', values='est_counts').astype(int).rename_axis('gene_symbol')
        count_dataframe.index = [x.split('.')[0] for x in count_dataframe.index]

        # Collapse
        if collapse:
            count_dataframe = count_dataframe.merge(annotation_dataframe[['Transcript stable ID', 'Gene name']].drop_duplicates('Transcript stable ID'), left_index=True, right_on='Transcript stable ID').drop('Transcript stable ID', axis=1).groupby('Gene name').sum().rename_axis('gene_symbol')


    # other
    else:
        raise ValueError('Please indicate STAR or kallisto as the alignment method.')

    # Return
    return count_dataframe
