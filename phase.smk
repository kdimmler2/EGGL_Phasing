localrules: vcf_list

rule all:
	input:
		expand('phasing/unphased/chr{chrom}/output.vcf.gz', chrom=[str(i) for i in range(1,32)] + ['X']),
		expand('phasing/unphased/chr{chrom}/output.vcf.gz.tbi', chrom=[str(i) for i in range(1,32)] + ['X']),
		expand('phasing/phased/chr{chrom}/output.vcf.gz', chrom=[str(i) for i in range(1,32)] + ['X']),
		expand('phasing/phased/chr{chrom}/output.vcf.gz.tbi', chrom=[str(i) for i in range(1,32)] + ['X']),
		'phasing/phased/vcfs.list',
		'goldenPath.Phased.vcf.gz',
		'goldenPath.Phased.vcf.gz.tbi'

rule split:
	input:
		vcf_to_phase = config['vcf_to_phase'],
		vcf_index = config['vcf_index']
	output:
		split_vcf = 'phasing/unphased/chr{chrom}/output.vcf.gz',
		split_tbi = 'phasing/unphased/chr{chrom}/output.vcf.gz.tbi'
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

rule phase:
	input:
		split_vcf = {rules.split.output.split_vcf}
	output:
		phased_vcf = 'phasing/phased/chr{chrom}/output.vcf.gz',
		indexed_vcf = 'phasing/phased/chr{chrom}/output.vcf.gz.tbi'
	threads: 24
	resources:
		time    = 720,
		mem_mb  = 200000,
		cpus    = 4
	params:
		map_directory = config['map_directory']
	shell:
		'''
			java -Xmx200g -jar beagle.22Jul22.46e.jar \
				gt={input.split_vcf} \
				chrom=chr{wildcards.chrom} \
				map={params.map_directory}/BEAGLE_Averaged_ECA{wildcards.chrom}_map.txt \
				out='phasing/phased/chr{wildcards.chrom}/output' \
				nthreads=24

			gatk IndexFeatureFile -I {output.phased_vcf}
		'''

rule vcf_list:
	input:
		vcfs = expand('phasing/phased/{chrom}/output.vcf.gz', chrom=['chr' + str(i) for i in range(1, 32)] + ['chrX'])
	output:
		vcf_list = 'phasing/phased/vcfs.list'
	threads: 4
	resources:
		time    = 60,
		mem_mb  = 24000,
		cpus    = 4
	run:
		outfile = open(output.vcf_list, 'wt')
		chrom=['chr' + str(i) for i in range(1,32)] + ['chrX']
		for chrm in chrom:
			print('phasing/phased/' + chrm + '/output.vcf.gz', file=outfile)
		outfile.close()

rule gather_vcfs:
	input:
		vcf_list = {rules.vcf_list.output.vcf_list}
	output:
		vcf = 'goldenPath.Phased.vcf.gz',
		tbi = 'goldenPath.Phased.vcf.gz.tbi'
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




