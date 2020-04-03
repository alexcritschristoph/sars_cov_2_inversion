### Reproducible analyses for rejecting rare genomic inversions in SARS-CoV-2

**Conclusion from below: I believe the recombination events seen here are sequencing artifacts and do not represent real intrapatient viral variation.**


A recent preprint proposed genomic evidence for intra-patient recombination of the SARS-CoV-2 virus, available at the link below:
https://www.biorxiv.org/content/10.1101/2020.03.27.009480v1

In particular, the authors proposed 5 genomic recombination inversions based on denovo assemblies of very small (<200 bp) contigs in SPADES. The hypothesis was that while these recombinations are not the dominant genotype in the samples (they are at very low coverage), they represent some small fraction of the viral population that underwent recombination during infection of the patient. The authors make the case that these inversions are indicative of a high frequency of recombination for the virus that may have epidemiological implications.

Two of these inversions were well supported by > 20 reads and occurred at similar locations in the S protein gene in two different patients (the location for both of these inversions is about 23956-24088 in the reference genome available in this repository). This was an intriguing possibility, especially considering they would likely have an effect on the S protein, and worth investigating further.

I copied the inversion sequences from their paper, and these two proposed inversion sequences are in the file `inversion_sequences.fna` in this repository.

The WGS reads for the two patient samples are available at the links below:
https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10903401
https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10903402

As far as I can tell, these samples were deeply sequenced without any amplification or targeted capture - human and bacterial reads were (mostly) removed. 

The first thing I wanted to do was check to see if I could replicate the preprint's main result, which was assembling a denovo contig with the proposed inversion. This was easy enough with two assemblers, megahit and metaSPADES:
```
megahit -1 SRR10903401.1.fq -2 SRR10903401.2.fq -t 16 -o megahit_3401
megahit -1 SRR10903402.1.fq -2 SRR10903402.2.fq -t 16 -o megahit_3402

metaspades.py -1 SRR10903401.1.fq -2 SRR10903401.2.fq -o spades_3401 -t 16
metaspades.py -1 SRR10903402.1.fq -2 SRR10903402.2.fq -o spades_3402 -t 16
```

We can then BLAST results of these assemblies against the reference:

```
alexcc@biotite:$ blastn -query megahit_3401/final.contigs.fa -subject reference.fna -outfmt 6

k141_718	MN908947	100.000	29846	0	0	111	29956	2	29847	0.0	55116
k141_718	MN908947	100.000	116	0	0	1	116	190	75	1.35e-56	215
```

```
alexcc@biotite:$ blastn -query megahit_3402/final.contigs.fa -subject reference.fna -outfmt 6
k141_16	MN908947	100.000	123	0	0	106	228	19769	19891	3.53e-62	228
k141_16	MN908947	96.471	85	1	2	36	118	19886	19802	1.70e-35	139
k141_994	MN908947	99.997	29837	1	0	88	29924	31	29867	0.0	55094
k141_994	MN908947	97.895	95	2	0	1	95	202	108	1.45e-41	165
```

For both samples, we don't see any contig matching the proposed inversion at the ~24 Kb location. The same was true of the metaSPADES analysis - so, I was unable to replicate obtaining the inversion by denovo assembly.

However, I was curious as to whether the authors had used different assembly settings that allowed them to assemble the inversion (they don't specify). And, I was curious about the reads that they show supporting the inversions - they do show many reads mapping to the inversion in their figure. 

**However, this was paired read sequencing data, and the locations of the read pairs for the reads mapping to the inversion are not shown**. These read pairs could provide crucial information about the validity of the inversion - if it was real, these reads should often be paired with nearby reads on the genome of the virus, as the read pair would have come from whole viral genomes that contained the inversion.

So, I modified two versions of the reference, *adding the proposed inversion sequences in to the reference genomes*. These are available in `inverted_reference.fna` and `inverted_reference2.fna` for each of the two (slightly different) inversions. To map reads to these references, first, I trim reads for quality and adapters:
```
sickle pe -f SRR10903401.1.fq -r SRR10903401.2.fq -t sanger -o SRR10903401.trim.1.fq -p SRR10903401.trim.2.fq -s SRR10903401.trim.unk.fq
sickle pe -f SRR10903402.1.fq -r SRR10903402.2.fq -t sanger -o SRR10903402.trim.1.fq -p SRR10903402.trim.2.fq -s SRR10903402.trim.unk.fq
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o SRR10903401.cutadapt.1.fq -p SRR10903401.cutadapt.2.fq SRR10903401.trim.1.fq SRR10903401.trim.2.fq
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o SRR10903402.cutadapt.1.fq -p SRR10903402.cutadapt.2.fq SRR10903402.trim.1.fq SRR10903402.trim.2.fq
```

And then I map reads with bowtie2:
```
bowtie2 -1 SRR10903401.cutadapt.1.fq -2 SRR10903401.cutadapt.2.fq -x inverted_reference.fna -p 16 > 3401_inversion.sam
bowtie2 -1 SRR10903402.cutadapt.1.fq -2 SRR10903402.cutadapt.2.fq -x inverted_reference2.fna -p 16 > 3402_inversion.sam
```

[Figure1]: 3401_inversion.png
[Figure2]: 3401_ref.png
[Figure3]: 3402_inversion.png


When we look at the proposed 3401 sample inversion site, this is what we see:

![Figure 1][Figure1]


Crucially, I've colored reads here by their *distances to their pairs*. Dark red reads have no pair; bright red reads have a pair, but are in the backwards orientation (inconsistent with an inversion); and light red to brown reads have perfectly overlapping pairs. Green reads have pairs that are occurring at expected distances - a few dozen up to maybe a couple hundred bp away. 

**As we can see, none of the read pairs overlapping the inversion have pairs nearby on the viral genome, indicating that there is no evidence that this inversion was seen in viral genomes, and is likely artifactual.**

What does this site look like in the correct reference? Completely different:

![Figure 2][Figure2]


The same is true for the 3402 sample inversion:

![Figure 3][Figure3]

**Conclusion: I believe the recombination events seen here are sequencing artifacts and do not represent real intrapatient viral variation.**
