import os

current_directory = os.getcwd()
files = os.listdir(current_directory)

for file_name in files:
    if file_name.startswith('Averaged_ECA') and file_name.endswith('_map.txt'):
        infile = open(file_name, 'rt')
        line = infile.readline()
        if len(file_name) > 21:
            chrom = file_name[12:14]
            outfile = open('BEAGLE_Averaged_ECA' + str(chrom) + '_map.txt', 'wt')
            for line in infile:
                line = line.rstrip()
                split = line.split('\t')
                print('chr' + str(chrom), '.', split[2], split[0], sep='\t', file=outfile)
        elif len(file_name) == 21:
            chrom = file_name[12]
            outfile = open('BEAGLE_Averaged_ECA' + str(chrom) + '_map.txt', 'wt')
            for line in infile:
                line = line.rstrip()
                split = line.split('\t')
                print('chr' + str(chrom), '.', split[2], split[0], sep='\t', file=outfile)
