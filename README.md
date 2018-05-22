# de.STAIR RNA-Seq scripts

Place for some of our scripts.

For an introduction to git you might want to explore a [Github tutorial](https://www.atlassian.com/git/tutorials/comparing-workflows)

Feel free to contribute, fork or add! For editing mark down, please see [cheat sheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)


### Source code

R code written and executed in RStudio (R version 3.4.3)

Python code written in Spyder3 (Python version 3.6)

Unix stuff executed with (Ubuntu 16)


### de.STAIR

See project details at [de.STAIR project page](http://destair.bioinf.uni-leipzig.de/)


### Tools

* Demultiplexing [flexbar](https://github.com/seqan/flexbar)
* Quality check with good old [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [qualimap](http://qualimap.bioinfo.cipf.es/)
* Pre-processing [fastx](http://hannonlab.cshl.edu/fastx_toolkit/) - alternatively [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or [trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* Mapping [kallisto](https://pachterlab.github.io/kallisto/download.html)
  * See for reference genomes: [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) or [Ensembl](http://www.ensembl.org/info/data/ftp/index.html)
* DE Analysis and statistics with kallisto sub-sequent tool [sleuth](https://pachterlab.github.io/kallisto/download.html)
  * Preferred setup for larger data sets, because kallisto is incredibly fast
  * Sleuth needs R/Rstudio
  * Alternatively using traditional mappers like [Star](https://github.com/alexdobin/STAR) or [Hisat2](http://www.ccb.jhu.edu/software/hisat/index.shtml)
  * [TRAPLINE](https://usegalaxy.org/u/mwolfien/p/trapline---manual) only for showcases =P
