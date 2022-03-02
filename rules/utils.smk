import os


def init():
    if ("work_dir" not in config) or config["work_dir"] == "":
        config["work_dir"] = os.path.abspath(os.curdir) + "/"
        print(config["work_dir"])

    if not config["work_dir"].endswith("/"):
        config["work_dir"] = config["work_dir"] + "/"
    if ("contig_dir" not in config) or config["contig_dir"] == "":
        config["contig_dir"] = config["work_dir"] + "contigs/"
    if ("primary_dir" not in config) or config["primary_dir"] == "":
        config["primary_dir"] = config["work_dir"] + "primaries/"
    if not config["contig_dir"].endswith("/"):
        config["contig_dir"] = config["contig_dir"] + "/"
    if not config["primary_dir"].endswith("/"):
        config["primary_dir"] = config["primary_dir"] + "/"

    contig_names = []
    primary_names = []
    assms = []
    for assm, assm_info in config["gap_conf"].items():
        assms.append(assm)
        for tool in assm_info["contigs"]:
            contig_names.append(tool)
        for tool in assm_info["primary"]:
            primary_names.append(tool)
    config["assm_names"] = assms
    config["contig_names"] = contig_names
    config["primary_names"] = primary_names
    print(config["contig_names"])
    print(config["primary_names"])


init()


def get_final():
    return config[""]
