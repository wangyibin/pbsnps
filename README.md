# pbsnp
## Instroduce
### The pipeline of pacbio reads snps calling by gatk4
 
## Dependencies
Following is a list of thirty-party progams that will be used in pbsnp pipeline.
- [minimap2](https://github.com/lh3/minimap2)
- [gatk4 v4.0+](https://software.broadinstitute.org/gatk/)
- [picard](https://github.com/broadinstitute/picard)
- [samtools](https://github.com/samtools/samtools)
[](- [snakemake](https://snakemake.readthedocs.io/)(optional))

## Install
1. Download
```bash
git clone https://github.com/wangyibin/pbsnps.git
```
2. Configure
```bash
export PATH=/path-to-pbsnps:$PATH
export picard="/path-to-picard/picard.jar"
```

## Usage
```bash
pacbio-snps-gatk4.sh test.fasta test.fq.gz test
```
or run with multithreads
```bash
pacbio-snps-gatk4.sh test.fasta test.fq.gz test 8
```

## Reference
[Comprehensive variant detection in a human genome with PacBio high-fidelity reads. William J. Rowell, Paul Peluso, John Harting, Yufeng Qian, Aaron Wenger, Richard Hall, David R. Rank. PacBio, 1305 O'Brien Drive, Menlo Park, CA 94025.](https://www.pacb.com/wp-content/uploads/Rowell-CSHLBioData-2018-Comprehensive-Variant-Detection-in-a-Human-Genome-with-PacBio-High-Fidelity-Reads.pdf)

