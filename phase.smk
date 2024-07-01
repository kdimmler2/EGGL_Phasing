localrules: vcf_list

chroms = [str(i) for i in range(1,32)] + ['X']

rule all:
    input:
        expand('results/unphased/split_vcf/chr{chrom}/output.vcf.gz', chrom=chroms),
        expand('results/unphased/split_vcf/chr{chrom}/output.vcf.gz.tbi', chrom=chroms),
        expand('results/unphased/split_multis/chr{chrom}/output.vcf.gz', chrom=chroms),
        expand('results/unphased/split_multis/chr{chrom}/output.vcf.gz.tbi', chrom=chroms),
        expand('results/phased/chr{chrom}/output.vcf.gz', chrom=chroms),
        expand('results/phased/chr{chrom}/output.vcf.gz.tbi', chrom=chroms),
        'results/phased/vcfs.list',
        'results/final/goldenPath.Phased.vcf.gz',
        'results/final/goldenPath.Phased.vcf.gz.tbi'

rule split_vcf:
    input:
        vcf_to_phase = config['vcf_to_phase'],
        vcf_index = config['vcf_index']
    output:
        split_vcf = 'results/unphased/split_vcf/chr{chrom}/output.vcf.gz',
        split_tbi = 'results/unphased/split_vcf/chr{chrom}/output.vcf.gz.tbi'
    threads: 4
    resources:
        time    = 360,
        mem_mb  = 60000,
        cpus    = 4
    shell:
        '''
            bcftools view \
            -r chr{wildcards.chrom} \
            -Oz -o {output.split_vcf} \
            {input.vcf_to_phase}

            gatk IndexFeatureFile -I {output.split_vcf}
        '''

rule split_multis:
    input:
        split_vcf = 'results/unphased/split_vcf/chr{chrom}/output.vcf.gz', 
        split_tbi = 'results/unphased/split_vcf/chr{chrom}/output.vcf.gz.tbi',
    output:
        split_multi = 'results/unphased/split_multis/chr{chrom}/output.vcf.gz',
        split_multi_tbi = 'results/unphased/split_multis/chr{chrom}/output.vcf.gz.tbi'
    threads: 4
    resources:
        time    = 360,
        mem_mb  = 60000,
        cpus    = 4
    shell:
        '''
            bcftools norm -m - -o {output.split_multi} {input.split_vcf}

            gatk IndexFeatureFile -I {output.split_multi}
        '''
rule phase:
    input:
        split_multi = 'results/unphased/split_multis/chr{chrom}/output.vcf.gz',
    output:
        phased_vcf = 'results/phased/chr{chrom}/output.vcf.gz',
        indexed_vcf = 'results/phased/chr{chrom}/output.vcf.gz.tbi'
    threads: 24
    resources:
        time    = 5760,
        mem_mb  = 200000,
        cpus    = 4
    params:
        map_directory = config['map_directory']
    shell:
        '''
            java -Xmx200g -jar beagle.22Jul22.46e.jar \
                gt={input.split_multi} \
                chrom=chr{wildcards.chrom} \
                map={params.map_directory}/BEAGLE_Averaged_ECA{wildcards.chrom}_map.txt \
                out='results/phased/chr{wildcards.chrom}/output' \
                nthreads=24

            gatk IndexFeatureFile -I {output.phased_vcf}
        '''

rule vcf_list:
    input:
        vcfs = expand('results/phased/{chrom}/output.vcf.gz', chrom=['chr' + str(i) for i in range(1, 32)] + ['chrX'])
    output:
        vcf_list = 'results/phased/vcfs.list'
    threads: 4
    resources:
        time    = 60,
        mem_mb  = 24000,
        cpus    = 4
    run:
        outfile = open(output.vcf_list, 'wt')
        chrom=['chr' + str(i) for i in range(1,32)] + ['chrX']
        for chrm in chrom:
            print('results/phased/' + chrm + '/output.vcf.gz', file=outfile)
        outfile.close()

rule gather_vcfs:
    input:
        vcf_list = {rules.vcf_list.output.vcf_list}
    output:
        vcf = 'results/final/goldenPath.Phased.vcf.gz',
        tbi = 'results/final/goldenPath.Phased.vcf.gz.tbi'
    threads: 4
    resources:
        time    = 360,
        mem_mb  = 24000,
        cpus    = 4
    shell:
        '''
            bcftools concat \
            -Oz -o {output.vcf} \
            -f {input.vcf_list}

            gatk IndexFeatureFile -I {output.vcf}
        '''




