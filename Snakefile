configfile: "conf/config.yaml"
configfile: "conf/scaffold.yaml"
include: "conf/software.smk"
include: "rules/utils.smk"
include: "rules/gapless.smk"

wildcard_constraints:
    assm="|".join(config["assm_names"]),
    tool_pri="|".join(config["primary_names"]),
    tool_contig="|".join(config["contig_names"])


def get_all(wildcards):
    res = []
    for assm, assm_info in config["gap_conf"].items():
        for tool in assm_info["primary"]:
            res.append(config["work_dir"] + "final/{assm}.{tool_pri}/{assm}.{tool_pri}.closeGap.fa".format(
                assm=assm,tool_pri=tool))
    return res


rule all:
    input:
        get_all
