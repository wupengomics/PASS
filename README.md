# PASS

- **Description**: a **P**roteomics **A**lternative **S**plicing **S**creening pipeline.


<img src="./workflow.png" width = "80%" height = "80%" />


- **Version**: 1.0.0


- **System requirements**

  - Install [perl](https://www.perl.org)
  - Install [R](https://www.r-project.org)
  - Install [proABMr R package](https://bioconductor.org/packages/release/bioc/html/proBAMr.html)
  - Install [python](https://www.python.org)
  - Install [EMBOSS](https://emboss.sourceforge.net/download/)
  - Install [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/), unless other alignment tool is used
  - Install [Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml), unless other alignment tool is used
  - Install [Cufflinks](https://cole-trapnell-lab.github.io/cufflinks/), unless other assemble tool is used
  - Install [MaxQuant](https://www.coxdocs.org/doku.php?id=maxquant:common:download_and_installation), unless other MS searching tool is used


- **processRNASEQ**:

    - **Description**: Align RNA-Seq reads to the reference genome and reconstruct transcripts.

    - **Usage**: `perl processRNASEQ.pl [options] [-g <genome>] [-f <.gtf>] [-r <.fastq>]`

        ```
	    -g    Genome bowtie2 index name.
	    -f    Gene annotation file, .gtf format.
	    -r    File names for sequencing reads, .fastq format.
	          - Compressed files (.fastq.gz) are also supported.
	          - Paired-end files separated by commas.
	    -t    Path to tophat, eg. /home/user/bin/tophat
	          - By default, we try to search tophat in system PATH.
	    -c    Path to cufflinks, eg. /home/user/bin/cufflinks
	          - By default, we try to search cufflinks in system PATH. 
	    -p    Number of used threads. [Default: 12]
	    -o    Output folder. [Default: ./RNASEQ_output]
	    -h    Help message.
        ```

    - **Example**: `perl processRNASEQ.pl -p 12 -g human.genome -f gencode.v27.gtf -r Sample_R1.fq.gz,Sample_R2.fq.gz`


- **sortGTF**:   

    - **Description**: Sort annotation file.

    - **Usage**: `perl sortGTF.pl [-f input<.fastq>]`

        ```
        -f    File name of gene annotation, .gtf format.
              - Recommend cufflinks to generate this file.  
        -h    Help message.
        ```

    - **Example**: `perl sortGTF.pl -f genes.cuff.gtf`


- **getORF**:

    - **Description**: Perform local alignment on clean reads.

    - **Usage**: `perl getORF.pl [-f <.gtf>] [-g <genome.fa>]`

        ```
        -f    File name of gene annotation, .gtf format.
              - Recommend cufflinks to generate this file. 
        -g    Reference genome file name, fasta format.
        -h    Help message.
        ```

    - **Example**: `perl getORF.pl -f gene.cuff.sorted.gtf -g human.genome.fa`

    - **Output**:
    
    	```
        genes.gtf             Must add the CDS information
        transcript.fa         The sequence header includes protein ID, the corresponding gene ID and the CDS position. 
                              - Example: >CUFF.11.2_1|ENSMUSG00000033845.9|CDS:63-218
        protein.fa            The sequence header includes protein ID and the corresponding gene ID. 
                              - Example: >CUFF.11.2_1|ENSMUSG00000033845.9
	```

- **searchMS**:

    - **Description**: Search MS file against protein sequence database.

    - **Usage**: `perl searchMS.pl [-b <MaxQuantCmd.exe>] [-p <mqpar.xml>]`

        ```
        -p    Parameter configuration file.
              - Recommend to preconfigure the mqpar.xml file in MaxQuant GUI.
        -b    Path to MaxQuantCmd.exe eg. /home/user/MaxQuant/bin/MaxQuantCmd.exe.
              - By default, we try to search MaxQuantCmd.exe in system PATH.
        -h    Help message.
        ```

    - **Example**: `perl searchMS.pl -b /usr/local/bin/MaxQuantCmd.exe -p mqpar.xml`


- **preparePSM**:

    - **Description**: Prepare PSM from MaxQuant msms.txt results

    - **Usage**: `perl preparePSM.pl [-m <msms.txt>]`

        ```
        -m    File name of peptide spectral matches.
              - Recommend MaxQuant to generate this file. 
        -h    Help message.
        ```

    - **Example**: `perl preparePSM.pl -m sample.msms.txt`

    - **Output**
        - Default output file: sample.psm.tab
        - Output format:
        The output requires 9 columns:
	
          No.|Column
          -|-
          1|spectrum
          2|spectrumNativeID
          3|assumed_charge
          4|hit_rank
          5|peptide
          6|num_missed_cleavages
          7|mvh
          8|modification
          9|NTT


- **generateSAM**:

    - **Description**: Convert peptide spectal matches to alignment file.

    - **Usage**: `perl generateSAM.pl [-m <psm.tab>] [-f <genes.gtf>] [-t <transcript.fa>] [-p <protein.fa>] [-o <.sam>]`

        ```
        -m    Peptide spectral matches.
              - Original msms files require conversion format via preparePSM.pl.
        -f    File name of gene annotation, .gtf format.
        -t    File name of transcript sequences, .fa format.
        -p    File name of protein sequences, .fa format.
        -o    Output file, .sam format.
        -h    Help message.
        ```

    - **Example**: `perl generateSAM.pl -m sample.psm.tab -f genes.gtf -t transcript.fa -p protein.fa -o sample.sam`

    - **Output**:
        - Output .sam format could reference to [proBAM](http://www.psidev.info/probam) format.


- **detectAS**:

    - **Description**: Detect AS events from annotation and alignment files.

    - **Note**: This function code is sourced from [MATS](https://rnaseq-mats.sourceforge.net).

    - **Usage**: `python detectAS.py [<genes.gtf>] [<outputPrefix>] [<SAMfiles>]`

        ```
        genes.gtf           File name of gene annotation, .gtf format.
        outputPrefix        Prefix of output AS event files
        SAMfiles            File name of alignment files. 
                            - Multiple samples are separated by commas.
        ```

    - **Example**: `python detectAS.py genes.gtf fromGTF sample1.sam,sample2.sam ./temp`

    - **Output**
      - Output files list:
      ```
      fromGTF.SE.txt
      fromGTF.MXE.txt
      fromGTF.A5SS.txt
      fromGTF.A3SS.txt
      fromGTF.RI.txt
      fromGTF.AFE.txt
      fromGTF.ALE.txt
      ```
      - Format:
      
      Column|Description
      -|-
      ID|AS event id
      GeneID|Gene id
      geneSymbol|Gene name
      chr|Chromosome
      strand|Strand of the gene
      event coordinates|Coordinates of the exons in the events, with multiple columns
      

- **quantifyAS**:

    - **Description**: Quantify AS by the alignment file.
    - **Note**: This function code is sourced from [MATS](https://rnaseq-mats.sourceforge.net).

    - **Usage**: `python quantifyAS.py [<ASPrefix>] [<SAMfiles>] [<outdir>]`
        ```
        ASPrefix         Prefix of AS event files
        SAMfiles         File name of alignment files. 
                         - Multiple samples are separated by commas.
        outdir           Folder of output files.
        ```

    - **Example**: `python quantifyAS.py fromGTF sample1.sam,sample2.sam ./`

    - **Output**
      - Output files list:
      ```
      JCEC.proteome.SE.txt
      JCEC.proteome.MXE.txt
      JCEC.proteome.A5SS.txt
      JCEC.proteome.A3SS.txt
      JCEC.proteome.RI.txt
      JCEC.proteome.AFE.txt
      JCEC.proteome.ALE.txt
      ```
      
      - Format:
      
      Column|Description
      -|-
      ID|AS event id
      IC|inclusion counts
      SC|skipping counts
      

- **Contact**:

    Peng Wu; wupeng1@ihcams.ac.cn
