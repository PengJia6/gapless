# Gapless

Gapless is developed for genome assemblies merging with a hybrid assembly strategy 


## Authors
  * Pengjia (pengjia@stu.xjtu.edu.cn)


## License

Gapless is free for non-commercial use
by academic, government, and non-profit/not-for-profit institutions. A
commercial version of the software is available and licensed through
Xi’an Jiaotong University. For more information, please contact with
Peng Jia (pengjia@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn).


## Dependence

* System: Unix like system (test on ubuntu and centOS)
* Language: Python and snakemake
* Packages: snakemake, pysam, numpy and pandas 
* Software: [bwa](https://github.com/lh3/bwa), [minimap2](https://github.com/lh3/minimap2), [samtools](http://www.htslib.org/), [minidot](https://github.com/thackl/minidot), [jellyfish](https://github.com/gmarcais/Jellyfish), [ragtag](https://github.com/malonge/RagTag), and [seqtk](https://github.com/lh3/seqtk)








## How to run Gapless?


### 1. Prepare you environment and install dependence   

   You can install the dependent packages and software with conda.     

### 2. Download the pipeline and config the environment
   Download the pipeline. 

```shell script
cd /path/to/gapless
git clone https://github.com/PengJia6/gapless.git
```
Changed the software or packages path in conf/software.smk.  
Changed the parameters of the pipeline in conf/config.yaml.
Changed the input path of you scaffolds and contigs in conf/scaffold.yaml.

### 3. Run the pipeline. 

```shell
cd /path/to/gapless
snakemake  -s Snakefile -j 10 -k --ri -k 
```

If you want to run the pipeline in other directory, please copy the conf file in your directory and run with following command:

```shell
cd /m/work/directory
cp -r /path/to/gapless/gapless/conf .
snakemake -s /path/to/gapless/gapless/Snakefile -j 10 -k --ri 
```



## Contact

If you have any questions, please contact with Peng Jia (pengjia@stu.xjtu.edu.cn) at Xi'an Jiaotong University.
