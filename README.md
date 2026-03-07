# EGGL_Phasing

This is a Snakemake workflow to phase whole genome sequence data using Beagle. This was designed for equine, so the publically available equine recombination maps are included. However, these can be replaced with any other maps.

This workflow splits a VCF into indiviudal chromosomes, phases them, and merges them back into a single vcf.

# Installation

```bash
git clone https://github.com/kdimmler2/EGGL_Phasing.git
cd EGGL_Phasing

conda env create -f phase.yaml
conda activate phase
```

The path to the VCF and recombination maps can be indicated in the config.yaml
