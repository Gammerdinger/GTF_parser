# GTF_parser.pl

## Aim

Take "complete" transcripts from Transdecoder and parse the GTF file using a FASTA file that has a subset of the transcripts in the GTF file. In the process, GTF_parser.pl will also remove entries that have formatting issues, so that the output GTF should be a clean output for [SnpEff](https://pcingola.github.io/SnpEff/)'s `build` database command. I am not sure why Transdecoder has these errors, but sometimes transcripts will have:

1. Transcripts ending in the GTF file after the end of the chromosome/scaffold/contig
2. Start codons that are not three bases long
3. Stop codons that are not three bases long

Also, SnpEff doesn't like having the "Gene" line in the GTF file for multiple transcripts. SnpEff also doesn't seem to like to have the annotations for the 5'UTR, 3'UTR, start codon and stop codon. I was running SnpEffv4.3t for this work, so subsequent versions may handle these annotations differently. However, GTF_parser.pl cleans all of these issues for at least SnpEffv4.3t.

## Usage

An example of a command to use GTF_parser.pl is:

```
perl GTF_parser.pl  --gtf_file=gtf_file.gtf --transcriptome_input=transcriptome_input.fasta --genome_input=genome_input.fasta --output_file=output_file.gtf
```

The  `--gtf_input` option accepts the GTF output of Transdecoder after it has been filtered by transcripts annotated as "complete".

The  `--transcriptome_input` option accepts your subsetted transcriptome FASTA file. I think (while I haven't explicitly tested this) if you have a transcriptome FASTA that contains all of the transcripts entries in the GTF input, it will still clean up the file and remove the problematic GTF entries and superfluous annotations that SnpEff doesn't like.

The `--genome_input` option accepts the genome FASTA file. We just need this so we know if transcript annotations run off of the chromosome/scaffold/contig.

The `--output_file` option accepts the cleaned GTF output filename.

An example of a GTF file input from Transdecoder might have entries that look like:

```
C1174700	transdecoder	5UTR	1077	1197	0	-	0	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	start_codon	1074	1076	0	-	0	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	CDS	855	1076	0	-	0	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	exon	855	1197	0	-	.	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	CDS	636	779	0	-	0	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	stop_codon	633	635	0	-	0	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	3UTR	628	632	0	-	0	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	exon	628	779	0	-	.	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	3UTR	2	563	0	-	0	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	exon	2	563	0	-	.	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	transcript	2	1197	0	-	.	gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700	transdecoder	gene	2	1197	0	-	.	gene_id "CUFF.9826"; Name "ORF type:complete len:123 (+),score=29.48";
```

Notice how this entry has `type:complete`. In order to parse out the "complete" transcripts prior to running GTF_parser.pl, I just ran:

```
grep "type:complete" transcoder_output.gtf > gtf_file.gtf
```

After GTF_parser.pl, this entry will look like:

```
C1174700        transdecoder    gene    2       1197    0       -       .       gene_id "CUFF.9826"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700        transdecoder    transcript      2       1197    0       -       .       gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700        transdecoder    CDS     855     1076    0       -       0       gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700        transdecoder    exon    855     1197    0       -       .       gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700        transdecoder    CDS     636     779     0       -       0       gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700        transdecoder    exon    628     779     0       -       .       gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
C1174700        transdecoder    exon    2       563     0       -       .       gene_id "CUFF.9826"; transcript_id "CUFF.9826.2.p2"; Name "ORF type:complete len:123 (+),score=29.48";
```
