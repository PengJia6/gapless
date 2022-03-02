# ======================================================================================================================
# Project: Project_Human_iharbor
# Script : scaffolding_PJ.smk TODO check 
# Author : Peng Jia
# Date   : 2021.07.21
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
### TODO List
## * documentation
## * get gaps from seqtk cutN

import random

configfile: "scaffold.yaml"
minimap2 = "/data/home/pengjia/mybin/minimap2"
samtools = "/data/home/pengjia/miniconda3/envs/default/bin/samtools"
bwa = "/data/home/pengjia/miniconda3/envs/assm/bin/bwa"
ref_path = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/collapse_assm/02_contig_eval/fastas/chm13.draft_v1.1_GRCh38_chrY.fa"
eps2pdf = "/data/home/pengjia/mybin/eps2pdf"
minidot = "/data/home/pengjia/mysoftware/assm/miniasm/minidot"
jellyfish = "/data/home/pengjia/miniconda3/envs/assm/bin/jellyfish"
ragtag = "/data/home/pengjia/miniconda3/envs/assm/bin/ragtag.py"
seqtk = "/data/home/pengjia/miniconda3/envs/assm/bin/seqtk"

dir_data = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/hap_assm_post/scaffolds_PJ_new_200K_prefix/"
dir_data = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/04_gap_close_v1.1/"
dir_primary = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/04_gap_close_v1.1/primary/"
dir_contigs = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/04_gap_close_v1.1/contigs/"
dir_assembly_ragtag = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/03_contigs_eval/ragtag/"
segment_length = 100
prefix_len = 200000
suffix_len = 200000


def get_pri(wildcards):
    # tool = config["gap_conf"][wildcards.assm][wildcards.hap]["primary"]
    return dir_assembly_ragtag + "{assm}.{hap}.{tool}.ragtag.scaffold.pj.fa".format(assm=wildcards.assm,
        hap=wildcards.hap,tool=wildcards.tool_pri)


def get_contig(wildcards):
    return dir_assembly_ragtag + "{assm}.{hap}.{tool}.ragtag.corContigs.pj.fa".format(assm=wildcards.assm,
        hap=wildcards.hap,tool=wildcards.tool_contig)


# def get_fq(wildcards):
#     return config["contigs"][wildcards.assm]["fq"]

tool_contigs = []
tool_primary = []
for assm, assm_info in config["gap_conf"].items():
    for hap, hap_info in assm_info.items():
        for tool in hap_info["contigs"]:
            tool_contigs.append(tool)
        for tool in hap_info["primary"]:
            tool_primary.append(tool)

assms = list(config["gap_conf"].keys())
wildcard_constraints:
    assm="|".join(assms),
    hap="hap1|hap2|collapsed|h1|h2",
    tool_pri="|".join(tool_primary),
    tool_contig="|".join(tool_contigs)


# tool="|".join(config[""])


# assms = ["CQ_hap1"]
#
# def get_all_contigs(wildcards):
#     res = []
#     for assm, assm_info in config["contigs"].items():
#         for hap, hap_info in assm_info.items():
#             for tool in hap_info["contig"]:
#                 res.append(dir_data + "contigs/{assm}.{hap}.{tool}/{assm}.{hap}.{tool}.ok".format(
#                     assm=assm,hap=hap,tool=tool))
#                 res.append(dir_data + "bwa_bams/{assm}.{hap}.{tool}.sort.bam".format(
#                     assm=assm,hap=hap,tool=tool))
#                 res.append(dir_data + "gaps/{assm}.{hap}.{tool}.gaps.info".format(
#                     assm=assm,hap=hap,tool=tool))
#                 res.append(dir_data + "close/{assm}.{hap}.{tool}.map".format(
#                     assm=assm,hap=hap,tool=tool))
# res.append(dir_data + "new_asm/{assm}.{hap}.{tool}/{assm}.{hap}.{tool}.fa".format(
#     assm=assm,hap=hap,tool=tool))
# return res


# def get_paired_fq(wildcards):
#     res = []
#     for assm, assm_info in config["contigs"].items():
#         for hap, hap_info in assm_info.items():
#             res.append(dir_data + "paired_fq/{assm}.{hap}.L.fq".format(assm=assm,hap=hap))
#     return res


def get_all(wildcards):
    res = []
    for assm, assm_info in config["gap_conf"].items():
        for hap, hap_info in assm_info.items():
            for tool in hap_info["primary"]:
                res.append(dir_data + "final/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.closeGap.fa".format(
                    assm=assm,hap=hap,tool_pri=tool))
    return res


rule all:
    input:
        get_all


# get_bwa()
# expand(dir_data + "kmer/{assm}.{hap}.kmer.fasta",assm="ChineseQuartet",hap=["hap1", "hap2"]),
# expand(dir_data + "kmer/{assm}.{hap}.kmer.bam",assm="ChineseQuartet",hap=["hap1", "hap2"]),
# expand(dir_data + "paired_sunk/{assm}.{hap}.L.fq",assm="ChineseQuartet",hap=["hap1", "hap2"]),
# expand(dir_data + "closed/{assm}.{hap}/{assm}.{hap}.ragtag.closeGap.fa",assm="ChineseQuartet",
#     hap=["hap1", "hap2", "collapsed"]),


# get_all_contigs_sunk_align


# get_all_contigs_sunk_align,
# get_all_contigs


# get_paired_fq,
# get_all_contigs,


# get_all_contigs


def get_contig_fa(wildcards):
    return config["contigs"][wildcards.assm][wildcards.hap]["contig"][wildcards.tool]


rule jellyfish:
    input:
        fa=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa",
        fai=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa.fai"
    output:
        kmer=dir_data + "kmer/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.kmer.fasta",
    threads: 20
    run:
        js = str(output.kmer)[:-5] + "js"
        tsv = str(output.kmer)[:-5] + "tsv"
        shell("{jellyfish} count -m {segment_length} -o {js} -t {threads} -s 3G {input.fa}")
        shell("{jellyfish} dump -c -t {js} > {tsv}")
        fasta = open(str(output.kmer),"w")
        linenum = 0
        for line in open(str(tsv)):
            lineinfo = line[:-1].split("\t")
            if lineinfo[1] != "1":
                continue
            linenum += 1
            fasta.write(">sunk_{cont}\n{substr}\n".format(cont=linenum,substr=lineinfo[0]))
        fasta.close()

#
# rule cutN:
#     input:
#         fa=dir_data + "contigs/{assm}.{hap}.{tool}/{assm}.{hap}.{tool}.scaffold.fa",
#     output:
#         fa=dir_data + "contigs/{assm}.{hap}.{tool}/{assm}.{hap}.{tool}.fa",
#     run:
#         shell("{seqtk} cutN -n 10 {input.fa} >{output.fa}")

rule get_gaps:
    input:
        dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa",
    output:
        dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.gap.bed",
    run:
        shell("{seqtk} cutN -n1 -g {input} >{output}")

rule sunk_location:
    input:
        kmer=dir_data + "kmer/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.kmer.fasta",
        ref=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa",
        ref_fai=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa.fai"
    output:
        kmer=dir_data + "kmer/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.kmer.bam",
    threads: 48
    run:
        ref = str(input.ref)
        shell("{bwa} mem -M -t {threads} {ref} {input.kmer} | "
              "{samtools} view -Shb -@ {threads} | "
              "{samtools} sort -@ {threads} -m 2G -T {output}_tmp -o {output} -O BAM ")


rule get_sunk_between_gaps:
    input:
        ref=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa",
        ref_fai=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa.fai",
        gap=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.gap.bed",
        kmer=dir_data + "kmer/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.kmer.bam",
        kmer_index=dir_data + "kmer/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.kmer.bam.bai",
    output:
        fq1=dir_data + "kmer/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.kmer.L.fq",
        fq2=dir_data + "kmer/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.kmer.R.fq",
    priority: 50
    run:
        import pysam

        sam_file = pysam.Samfile(str(input.kmer))
        fasta_file = pysam.Fastafile(str(input.ref))
        fasta_length = {j: i for i, j in zip(fasta_file.lengths,fasta_file.references)}
        fasta_file.close()
        fq1_file = open(str(output.fq1),"w")
        fq2_file = open(str(output.fq2),"w")
        # apg = str(input.fa)[:-3] + '.agp'
        gaps = {}
        for line in open(str(input.gap)):
            if line.startswith("#"): continue
            lineinfo = line[:-1].split("\t")
            # if lineinfo[4] in ["U", "N"]:
            #     gap_name = ".".join(lineinfo[:3])
            gap_chr, gap_start, gap_end = lineinfo[:3]
            gap_name = ".".join(lineinfo[:3])
            gaps[gap_name] = [gap_chr, int(gap_start), int(gap_end)]
        for gap, gap_info in gaps.items():
            if gap_info[1] - prefix_len > 0:
                L_start = gap_info[1] - prefix_len
                L_end = gap_info[1] - segment_length
            else:
                continue
            if gap_info[2] + prefix_len < fasta_length[gap_info[0]]:
                R_start = gap_info[2] + segment_length
                R_end = gap_info[2] + suffix_len
            else:
                continue
            idx = 0
            for read in sam_file.fetch(gap_info[0],L_start,L_end):
                if read.mapping_quality < 30: continue
                if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
                cigar_tuples = read.cigartuples
                if cigar_tuples[0][0] in [4, 5, 6, 3] or cigar_tuples[-1][0] in [4, 5, 6, 3]: continue
                idx += 1
                # left_sunks.append(read)
                read_name = gap + "_" + str(read.reference_start) + "_" + str(idx)
                left = read.query_sequence
                left_record = "@{read_name}\n{read_str}\n+\n{quals}\n".format(read_name=read_name,read_str=left,quals="!" * 100)
                fq1_file.write(left_record)
            idx = 0
            for read in sam_file.fetch(gap_info[0],R_start,R_end):
                if read.mapping_quality < 30: continue
                if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
                cigar_tuples = read.cigartuples
                if cigar_tuples[0][0] in [4, 5, 6, 3] or cigar_tuples[-1][0] in [4, 5, 6, 3]: continue
                idx += 1
                read_name = gap + "_" + str(read.reference_start) + "_" + str(idx)
                right = read.query_sequence
                right_record = "@{read_name}\n{read_str}\n+\n{quals}\n".format(read_name=read_name,read_str=right,quals="!" * 100)
                fq2_file.write(right_record)
        sam_file.close()
        fq1_file.close()
        fq2_file.close()


rule sunk_align_to_contig:
    input:
        ref=dir_contigs + "{assm}.{hap}.{tool_contig}/{assm}.{hap}.{tool_contig}.fa",
        ref_index=dir_contigs + "{assm}.{hap}.{tool_contig}/{assm}.{hap}.{tool_contig}.fa.fai",
        fq=dir_data + "kmer/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.kmer.{side}.fq",
    output:
        dir_data + "maps/{assm}.{hap}.{tool_pri}/{tool_contig}/{assm}.{hap}.{tool_pri}.{side}.{tool_contig}.sort.bam"
    threads: 48
    run:
        ref = str(input.ref)
        shell("{bwa} mem -M -t {threads} {ref} {input.fq} | "
              "{samtools} view -Shb -@ {threads} | "
              "{samtools} sort -@ {threads} -T {output}_tmp -o {output} -O BAM "
        )


rule get_best_sunk_two:
    input:
        bam_L=dir_data + "maps/{assm}.{hap}.{tool_pri}/{tool_contig}/{assm}.{hap}.{tool_pri}.L.{tool_contig}.sort.bam",
        bai_L=dir_data + "maps/{assm}.{hap}.{tool_pri}/{tool_contig}/{assm}.{hap}.{tool_pri}.L.{tool_contig}.sort.bam.bai",
        bam_R=dir_data + "maps/{assm}.{hap}.{tool_pri}/{tool_contig}/{assm}.{hap}.{tool_pri}.R.{tool_contig}.sort.bam",
        bai_R=dir_data + "maps/{assm}.{hap}.{tool_pri}/{tool_contig}/{assm}.{hap}.{tool_pri}.R.{tool_contig}.sort.bam.bai",
    output:
        dir_data + "maps/{assm}.{hap}.{tool_pri}/{tool_contig}/{assm}.{hap}.{tool_pri}.{tool_contig}.best.map",
    # dir_data + "sunk/{assm}.{hap}.{tool}.best.map"
    threads: 1
    run:
        import pysam

        out_info = open(str(output),"w")
        # from collections import Counter

        samfiles = {"L": pysam.AlignmentFile(str(input.bam_L)), "R": pysam.AlignmentFile(str(input.bam_R))}

        # out_sam = pysam.AlignmentFile(str(output),"wb",template=samfile)
        contig_length_dict = {i: j for i, j in zip(samfiles["L"].references,samfiles["L"].lengths)}
        reads_dict = {}
        for side, side_sam in samfiles.items():
            for read in side_sam:
                if read.is_unmapped: continue
                if read.get_tag("NM") > 15: continue
                # if read.mapping_quality < 20: continue
                # if read.has_tag("XA"): continue
                # if read.has_tag("SA"): continue
                cigar_tuples = read.cigartuples
                if cigar_tuples[0][0] in [4, 5, 6, 3] or cigar_tuples[-1][0] in [4, 5, 6, 3]: continue
                name = read.query_name
                gap_name = "_".join(name.split("_")[:2])
                contig_name = read.reference_name
                gap_chrom, gap_start, gap_end = gap_name.split(".")
                if side == "L":
                    dis2gap = int(gap_start) - int(read.query_name.split("_")[2])
                    if not read.is_reverse:
                        dis2contig_end = contig_length_dict[read.reference_name] - read.reference_start
                    else:
                        dis2contig_end = read.reference_end
                else:
                    dis2gap = int(read.query_name.split("_")[2]) - int(gap_end)
                    if not read.is_reverse:
                        dis2contig_end = read.reference_start
                    else:
                        dis2contig_end = contig_length_dict[read.reference_name] - read.reference_end
                if dis2contig_end - dis2gap < 10: continue
                if gap_name not in reads_dict:
                    reads_dict[gap_name] = {}
                if contig_name not in reads_dict[gap_name]:
                    reads_dict[gap_name][contig_name] = {"Forward": {"L": [], "R": []}, "Reversed": {"L": [], "R": []}}
                strand = "Reversed" if read.is_reverse else "Forward"
                reads_dict[gap_name][contig_name][strand][side].append(read)
        # if "chrX_RagTag.113059349.113073530" in reads_dict:
        #     print(reads_dict["chrX_RagTag.113059349.113073530"])
        for gap, gap_info in reads_dict.items():
            # total_ctg_reads_num = {"Forward": {"L": [], "R": []}, "Reversed": {"L": [], "R": []}}
            # per_ctg_reads_num = {}
            gap_score_dict = {}
            total_reads_num = {"L": 0, "R": 0}

            for ctg_name, ctg_reads in gap_info.items():
                for strand, strand_info in ctg_reads.items():
                    for side, side_info in strand_info.items():
                        total_reads_num[side] += len(side_info)

            for ctg_name, ctg_reads in gap_info.items():
                for strand, strand_info in ctg_reads.items():
                    ctg_read_num = {side: len(side_reads) for side, side_reads in strand_info.items()}
                    if gap == "chr11_RagTag.50914747.50914846":
                        print(ctg_name,ctg_read_num["L"],ctg_read_num["R"])
                    #     print(strand,total_reads_num)
                    #     print(ctg_name,ctg_read_num)
                    #     print("=================================")
                    L_min = max([20, total_reads_num["L"] * 0.1])
                    R_min = max([20, total_reads_num["R"] * 0.1])
                    if ctg_read_num["L"] < L_min or ctg_read_num["R"] < R_min: continue
                    score1 = (ctg_read_num["L"] / prefix_len + ctg_read_num["R"] / suffix_len) / 2
                    score2 = (ctg_read_num["L"] / total_reads_num["L"] + ctg_read_num["R"] / total_reads_num["R"]) / 2
                    gap_score_dict["_".join([ctg_name, strand])] = [score1, score2, (score1 + score2) / 2]
            # print(score1,score2)
            if len(gap_score_dict) < 1: continue
            max_score_ctg = max(gap_score_dict,key=lambda i: gap_score_dict[i][2])
            # print(max_score_ctg)
            best_ctg_name = "_".join(max_score_ctg.split("_")[:-1])
            best_ctg_strand = max_score_ctg.split("_")[-1]
            best_reads = reads_dict[gap][best_ctg_name][best_ctg_strand]
            best_read_paired = {}
            for side in best_reads:
                positions_reads = {}
                for read in best_reads[side]:
                    if read.reference_start not in positions_reads:
                        positions_reads[read.reference_start] = read
                    else:
                        continue
                positions_reads_order = list(sorted(positions_reads))
                best_reads_position = positions_reads_order[len(positions_reads_order) // 2]
                best_read_paired[side] = positions_reads[best_reads_position]
            scafford_name = best_read_paired['L'].query_name.split(".")[0]
            scafford_start = int(best_read_paired["L"].query_name.split("_")[2])
            scafford_end = int(best_read_paired["R"].query_name.split("_")[2])

            if best_ctg_strand == "Forward":
                contig_start = int(best_read_paired["L"].reference_start)
                contig_end = int(best_read_paired["R"].reference_start)
            else:
                contig_start = int(best_read_paired["R"].reference_end)
                contig_end = int(best_read_paired["L"].reference_end)
            if contig_end - contig_start < 10 or contig_end - contig_start > 10000000: continue

            out_info.write("{chr}\t{chr_start}\t{chr_end}\t{contig}\t{contig_start}\t{contig_end}\t{direction}\t{gap_name}\t{score}\n".format(
                chr=scafford_name,chr_start=scafford_start,chr_end=scafford_end,
                contig=best_ctg_name,contig_start=contig_start,contig_end=contig_end,
                gap_name=gap,direction=best_ctg_strand,score=gap_score_dict["_".join([best_ctg_name, best_ctg_strand])][
                    2]

            ))
        out_info.close()
        for side in samfiles:
            samfiles[side].close()


# rule bam_sort:
#     input:
#         dir_data + "sunk/{assm}.{hap}.{tool}.{side}.best.sunk.bam",
#     output:
#         dir_data + "sunk/{assm}.{hap}.{tool}.{side}.best.sunk.sort.bam",
#     run:
#         shell("{samtools} sort -m 2G -T {output}_tmp -o {output} -O BAM  {input}")


def get_reverse_complementary(my_str):
    alphblte = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
    new_str = ""
    for i in my_str[::-1]:
        new_str = new_str + alphblte[i]
    return new_str


def get_contigs_map(wildcards):
    tools = list(config["gap_conf"][wildcards.assm][wildcards.hap]["contigs"])
    return expand(dir_data + "maps/{assm}.{hap}.{tool_pri}/{tool_contig}/{assm}.{hap}.{tool_pri}.{tool_contig}.best.map",assm=wildcards.assm,hap=wildcards.hap,
        tool_contig=tools,tool_pri=wildcards.tool_pri)


def get_contigs_fa(wildcards):
    tools = list(config["gap_conf"][wildcards.assm][wildcards.hap]["contigs"])
    return expand(dir_contigs + "{assm}.{hap}.{tool_contig}/{assm}.{hap}.{tool_contig}.fa",assm=wildcards.assm,hap=wildcards.hap,tool_contig=tools)


rule merged_and_get_fa:
    input:
        closer=get_contigs_map,
        scafold=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa",
        contigs=get_contigs_fa
    output:
        dir_data + "final/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.closeGap.fa"
    # dir_data+ "sunk/{assm}.{hap}.{tool}.best.map"
    run:
        import pysam

        scaffold_fa = pysam.Fastafile(str(input.scafold))
        tools = list(config["gap_conf"][wildcards.assm][wildcards.hap]["contigs"])
        tools_samfile = {
            tool: pysam.Fastafile(str(dir_contigs + "{assm}.{hap}.{tool}/{assm}.{hap}.{tool}.fa").format(
                assm=wildcards.assm,hap=wildcards.hap,tool=tool
            )) for tool in tools}
        gaps_info_all = {}
        gap_filling_score = {}
        for tool in tools:
            for line in open(dir_data + "maps/{assm}.{hap}.{tool_pri}/{tool}/{assm}.{hap}.{tool_pri}.{tool}.best.map".format(
                    assm=wildcards.assm,hap=wildcards.hap,tool=tool,tool_pri=wildcards.tool_pri)):
                scaffold_id, scaffold_start, scafford_end, contig_id, contig_start, contig_end, direction, gap_name, score = \
                    line[:-1].split("\t")
                if gap_name not in gaps_info_all:
                    gaps_info_all[gap_name] = {}
                scaffold_start = int(scaffold_start)
                scafford_end = int(scafford_end)
                contig_start = int(contig_start)
                contig_end = int(contig_end)
                score = float(score)
                gaps_info_all[gap_name][tool] = [scaffold_id, scaffold_start, scafford_end, contig_id,
                                                 contig_start,
                                                 contig_end, direction, tool, score]
        gaps_info = {}
        for gap_name in gaps_info_all:
            # print(gaps_info_all[gap_name])
            this_gap_info = gaps_info_all[gap_name]
            gaps_info[gap_name] = gaps_info_all[gap_name][max(this_gap_info,key=lambda i: this_gap_info[i][-1])]

        new_gaps = {}
        for gap in gaps_info:
            gap_chr = gaps_info[gap][0]
            if gap_chr not in new_gaps:
                new_gaps[gap_chr] = []
            new_gaps[gap_chr].append(gaps_info[gap] + [gap])
        gap_order = {}
        for gap_chr, g_info in new_gaps.items():
            gap_start = {g[1]: g for g in g_info}
            g_info_sort_by_gap_start = sorted(list(gap_start.keys()))
            st = 0
            close_pos_this_time = []
            for s in g_info_sort_by_gap_start:
                if st < gap_start[s][1] < gap_start[s][2]:
                    close_pos_this_time.append([s, gap_start[s]])
                    st = gap_start[s][2]
                else:
                    pass
            gap_order[gap_chr] = close_pos_this_time
        scaffold_chroms_length_dict = {i: j for i, j in zip(scaffold_fa.references,scaffold_fa.lengths)}
        output_scafflod = open(str(output),"w")
        for chrom, chrom_len in scaffold_chroms_length_dict.items():
            header = ">{}\n".format(chrom)
            if chrom not in gap_order:
                # continue
                chrom_str = str(scaffold_fa.fetch(chrom)) + "\n"
            # print(chrom,0,chrom_len)
            else:
                chrom_str = ""
                close_start = [i[1][1] for i in gap_order[chrom]]
                close_end = [i[1][2] for i in gap_order[chrom]]
                no_close_start = [0] + [i[1][2] for i in gap_order[chrom]][:-1]
                no_close_end = [i[1][1] for i in gap_order[chrom]]
                contigs_info = [i[1] for i in gap_order[chrom]]
                # print("close_start",close_start)
                # print("close_end",close_end)
                for a, b, c, d, contig_info in zip(no_close_start,no_close_end,close_start,close_end,contigs_info):
                    chrom_str += scaffold_fa.fetch(chrom,a,b)
                    # print(contig_info[6])
                    if contig_info[6] == "Forward":
                        chrom_str += tools_samfile[contig_info[7]].fetch(contig_info[3],contig_info[4],contig_info[5])
                    # print("F")
                    # print(tools_samfile[contig_info[7]].fetch(contig_info[3],contig_info[4],contig_info[5]))
                    # print(scaffold_fa.fetch(chrom,c,d))
                    # print("================================================")
                    else:
                        chrom_str += get_reverse_complementary(
                            tools_samfile[contig_info[7]].fetch(contig_info[3],contig_info[4],contig_info[5]))
                        chrom_str += tools_samfile[contig_info[7]].fetch(contig_info[3],contig_info[4],contig_info[5])
                # print("R")
                # print(get_reverse_complementary(
                #     tools_samfile[contig_info[7]].fetch(contig_info[3],contig_info[4],contig_info[5])))
                # print(scaffold_fa.fetch(chrom,c,d))
                # print("================================================")
                # print(chrom,[a, b, c, d])
                if close_end[-1] < chrom_len:
                    chrom_str += scaffold_fa.fetch(chrom,close_end[-1],chrom_len)
                chrom_str += "\n"
            output_scafflod.write("{}".format(header))
            output_scafflod.write("{}".format(chrom_str))
        for tool in tools_samfile:
            tools_samfile[tool].close()
        scaffold_fa.close()
        output_scafflod.close()

# rule res_vis:
#     input:
#         scafold= dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa",
#         scafold_index= dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa.fai",
#         contig=dir_contigs + "{assm}.{hap}.{tool_contig}/{assm}.{hap}.{tool_contig}.fa",
#         contig_index=dir_contigs+ "{assm}.{hap}.{tool_contig}/{assm}.{hap}.{tool_contig}.fa.fai"
#     output:
#         paf=dir_data + "pafs/{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.{tool_contig}.paf",
#     threads: 48
#     run:
#         shell("{minimap2} -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 --secondary=no -t {threads} --eqx -Y "
#               "-O 5,56 -E 4,1 -B 5 -o{output.paf} {input.scafold} {input.contig}")
#
#
rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    threads: 48
    shell:
        "{samtools} index -@ {threads} {input}"

#
# rule bwa_index:
#     input:
#         dir_data + "contigs/{assm}.{hap}.{tool}/{assm}.{hap}.{tool}.fa",
#     output:
#         bwa=dir_data + "contigs/{assm}.{hap}.{tool}/{assm}.{hap}.{tool}.ok"
#     run:
#         # output_pre = str(output.bwa)[:-3]
#         # shell("ln -sf {input} {output_pre}.fa")
#         shell("{bwa} index {input}")
#         shell("touch {output}")

# rule samtools_index_fai:
#     input:
#         "{prefix}.fa"
#     output:
#         "{prefix}.fa.fai"
#     run:
#         shell("{samtools} faidx {input}")

rule bwa_index_primary:
    input:
        get_pri
    output:
        fa=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa",
        fai=dir_primary + "{assm}.{hap}.{tool_pri}/{assm}.{hap}.{tool_pri}.fa.fai"
    run:
        shell("ln {input} {output.fa}")
        shell("{bwa} index {output.fa}")
        shell("{samtools} faidx {output.fa}")


rule bwa_index_contig:
    input:
        get_contig
    output:
        fa=dir_contigs + "{assm}.{hap}.{tool_contig}/{assm}.{hap}.{tool_contig}.fa",
        fai=dir_contigs + "{assm}.{hap}.{tool_contig}/{assm}.{hap}.{tool_contig}.fa.fai"
    run:
        # fa_out = str(output)[:-4]
        shell("ln {input} {output.fa}")
        shell("{bwa} index {output.fa}")
        shell("{samtools} faidx {output.fa}")
