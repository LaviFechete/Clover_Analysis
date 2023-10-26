from gwf import Workflow

gwf = Workflow()

## Reference is a fasta file of the genes.

def bowtie2_index(reference, reference_name):
    inputs = [reference]
    outputs = [reference_name+".1.bt2"]
    options = {
        'cores': 1,
        'memory': '16g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = '''
    source activate HyLiTE

    bowtie2-build {ref} {ref_n}
    '''.format(ref=reference, ref_n=reference_name)

    return inputs, outputs, options, spec

def star_index(reference, output):
    inputs = [reference]
    outputs = [output]
    options = {"cores":8, "memory":"64g", "account":"NChain", "walltime": "12:00:00"}

    directory = "/".join(output.split("/")[:-1])

    spec = """
    source /com/extra/STAR/2.5.2b/load.sh
    STAR --runMode genomeGenerate --runThreadN 8 --genomeDir {dir} --genomeFastaFiles {ref} --limitGenomeGenerateRAM=64000000000 --genomeSAindexNbases 3
    """.format(ref=reference, dir=directory)

    return inputs, outputs, options, spec

reference_file = "./references/To/To.v5.gDNA.fasta"
index_ref_file = "./references/To/To.ref"

raw_gDNA = {"To": ["/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_occidentale_180bp_1.fastq.gz",
                   "/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_occidentale_180bp_2.fastq.gz"],
            "Tp": ["/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_pallescens_180_1.fastq.gz",
                   "/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_pallescens_180_2.fastq.gz"],
            "TrR": ["/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_repens_180bp_1.fastq.gz",
                    "/home/marnit/NChain/faststorage/20181120_clover_180bp_gDNA/Trifolium_repens_180bp_2.fastq.gz"]}

raw_RNA = {"To": {"F_TO_1_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_1_1/F_TO_1_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_1_1/F_TO_1_1_2.fq.gz"],
                  "F_TO_1_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_1_2/F_TO_1_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_1_2/F_TO_1_2_2.fq.gz"],
                  "F_TO_1_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_1_3/F_TO_1_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_1_3/F_TO_1_3_2.fq.gz"],
                  "F_TO_2_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_2_1/F_TO_2_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_2_1/F_TO_2_1_2.fq.gz"],
                  "F_TO_2_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_2_2/F_TO_2_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_2_2/F_TO_2_2_2.fq.gz"],
                  "F_TO_2_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_2_3/F_TO_2_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_2_3/F_TO_2_3_2.fq.gz"],
                  "F_TO_3_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_3_1/F_TO_3_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_3_1/F_TO_3_1_2.fq.gz"],
                  "F_TO_3_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_3_2/F_TO_3_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TO_3_2/F_TO_3_2_2.fq.gz"],
                  "F_TO_3_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TO_3_3/F_TO_3_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TO_3_3/F_TO_3_3_2.fq.gz"],
                  "F_TO_4_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TO_4_1/F_TO_4_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TO_4_1/F_TO_4_1_2.fq.gz"],
                  "F_TO_4_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TO_4_2/F_TO_4_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TO_4_2/F_TO_4_2_2.fq.gz"],
                  "F_TO_4_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TO_4_3/F_TO_4_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TO_4_3/F_TO_4_3_2.fq.gz"]},

           "Tp": {"F_TP_1_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_1_1/F_TP_1_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_1_1/F_TP_1_1_2.fq.gz"],
                  "F_TP_1_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_1_2/F_TP_1_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_1_2/F_TP_1_2_2.fq.gz"],
                  "F_TP_1_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_1_3/F_TP_1_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_1_3/F_TP_1_3_2.fq.gz"],
                  "F_TP_2_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_2_1/F_TP_2_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_2_1/F_TP_2_1_2.fq.gz"],
                  "F_TP_2_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_2_2/F_TP_2_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_2_2/F_TP_2_2_2.fq.gz"],
                  "F_TP_2_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_2_3/F_TP_2_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_2_3/F_TP_2_3_2.fq.gz"],
                  "F_TP_3_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_3_1/F_TP_3_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_3_1/F_TP_3_1_2.fq.gz"],
                  "F_TP_3_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_3_2/F_TP_3_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_3_2/F_TP_3_2_2.fq.gz"],
                  "F_TP_3_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_3_3/F_TP_3_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_3_3/F_TP_3_3_2.fq.gz"],
                  "F_TP_4_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_4_1/F_TP_4_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_4_1/F_TP_4_1_2.fq.gz"],
                  "F_TP_4_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_4_2/F_TP_4_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_4_2/F_TP_4_2_2.fq.gz"],
                  "F_TP_4_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_4_3/F_TP_4_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_05/raw_data/F_TP_4_3/F_TP_4_3_2.fq.gz"]},

           "88": {"F_88_1_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/F_88_1_1/F_88_1_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/F_88_1_1/F_88_1_1_2.fq.gz"],
                  "F_88_1_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/F_88_1_2/F_88_1_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/F_88_1_2/F_88_1_2_2.fq.gz"],
                  "F_88_1_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_1_3/F_88_1_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_1_3/F_88_1_3_2.fq.gz"],
                  "F_88_2_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_2_1/F_88_2_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_2_1/F_88_2_1_2.fq.gz"],
                  "F_88_2_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_2_2/F_88_2_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_2_2/F_88_2_2_2.fq.gz"],
                  "F_88_2_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_2_3/F_88_2_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_2_3/F_88_2_3_2.fq.gz"],
                  "F_88_3_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_3_1/F_88_3_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_3_1/F_88_3_1_2.fq.gz"],
                  "F_88_3_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_3_2/F_88_3_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_3_2/F_88_3_2_2.fq.gz"],
                  "F_88_3_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_3_3/F_88_3_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_3_3/F_88_3_3_2.fq.gz"],
                  "F_88_4_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_4_1/F_88_4_1_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_4_1/F_88_4_1_2.fq.gz"],
                  "F_88_4_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_4_2/F_88_4_2_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_4_2/F_88_4_2_2.fq.gz"],
                  "F_88_4_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_4_3/F_88_4_3_1.fq.gz",
                               "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_88_4_3/F_88_4_3_2.fq.gz"]},

           "S10": {"A_S1_1_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_1_1/A_S1_1_1_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_1_1/A_S1_1_1_2.fq.gz"],
                   "A_S1_1_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_1_2/A_S1_1_2_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_1_2/A_S1_1_2_2.fq.gz"],
                   "A_S1_1_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_1_3/A_S1_1_3_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_1_3/A_S1_1_3_2.fq.gz"],
                   "A_S1_2_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_2_1/A_S1_2_1_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_2_1/A_S1_2_1_2.fq.gz"],
                   "A_S1_2_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_2_2/A_S1_2_2_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_2_2/A_S1_2_2_2.fq.gz"],
                   "A_S1_2_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_2_3/A_S1_2_3_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_2_3/A_S1_2_3_2.fq.gz"],
                   "A_S1_3_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_3_1/A_S1_3_1_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_3_1/A_S1_3_1_2.fq.gz"],
                   "A_S1_3_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_3_2/A_S1_3_2_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_3_2/A_S1_3_2_2.fq.gz"],
                   "A_S1_3_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_3_3/A_S1_3_3_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_3_3/A_S1_3_3_2.fq.gz"],
                   "A_S1_4_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_4_1/A_S1_4_1_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_4_1/A_S1_4_1_2.fq.gz"],
                   "A_S1_4_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_4_2/A_S1_4_2_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_4_2/A_S1_4_2_2.fq.gz"],
                   "A_S1_4_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_4_3/A_S1_4_3_1.fq.gz",
                                        "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_01/raw_data/A_S1_4_3/A_S1_4_3_2.fq.gz"],
                   "F_S1_1_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_S1_1_1/F_S1_1_1_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_S1_1_1/F_S1_1_1_2.fq.gz"],
                   "F_S1_1_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_S1_1_2/F_S1_1_2_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_S1_1_2/F_S1_1_2_2.fq.gz"],
                   "F_S1_1_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_S1_1_3/F_S1_1_3_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_S1_1_3/F_S1_1_3_2.fq.gz"],
                   "F_S1_2_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_S1_2_1/F_S1_2_1_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_02/raw_data/F_S1_2_1/F_S1_2_1_2.fq.gz"],
                   "F_S1_2_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_2_2/F_S1_2_2_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_2_2/F_S1_2_2_2.fq.gz"],
                   "F_S1_2_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_2_3/F_S1_2_3_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_2_3/F_S1_2_3_2.fq.gz"],
                   "F_S1_3_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_3_1/F_S1_3_1_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_3_1/F_S1_3_1_2.fq.gz"],
                   "F_S1_3_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_3_2/F_S1_3_2_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_3_2/F_S1_3_2_2.fq.gz"],
                   "F_S1_3_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_3_3/F_S1_3_3_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_3_3/F_S1_3_3_2.fq.gz"],
                   "F_S1_4_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_4_1/F_S1_4_1_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_4_1/F_S1_4_1_2.fq.gz"],
                   "F_S1_4_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_4_2/F_S1_4_2_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_4_2/F_S1_4_2_2.fq.gz"],
                   "F_S1_4_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_4_3/F_S1_4_3_1.fq.gz",
                                "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_S1_4_3/F_S1_4_3_2.fq.gz"]},

           "Tienshan": {"F_TI_1_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_1_1/F_TI_1_1_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_1_1/F_TI_1_1_2.fq.gz"],
                        "F_TI_1_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_1_2/F_TI_1_2_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_1_2/F_TI_1_2_2.fq.gz"],
                        "F_TI_1_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_1_3/F_TI_1_3_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_1_3/F_TI_1_3_2.fq.gz"],
                        "F_TI_2_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_2_1/F_TI_2_1_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_2_1/F_TI_2_1_2.fq.gz"],
                        "F_TI_2_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_2_2/F_TI_2_2_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_2_2/F_TI_2_2_2.fq.gz"],
                        "F_TI_2_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_2_3/F_TI_2_3_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_2_3/F_TI_2_3_2.fq.gz"],
                        "F_TI_3_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_3_1/F_TI_3_1_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_3_1/F_TI_3_1_2.fq.gz"],
                        "F_TI_3_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_3_2/F_TI_3_2_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_3_2/F_TI_3_2_2.fq.gz"],
                        "F_TI_3_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_3_3/F_TI_3_3_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_3_3/F_TI_3_3_2.fq.gz"],
                        "F_TI_4_1": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_4_1/F_TI_4_1_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_4_1/F_TI_4_1_2.fq.gz"],
                        "F_TI_4_2": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_4_2/F_TI_4_2_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_4_2/F_TI_4_2_2.fq.gz"],
                        "F_TI_4_3": ["/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_4_3/F_TI_4_3_1.fq.gz",
                                     "/home/marnit/NChain/faststorage/WhiteCloverColdResponse/X201SC19112319-Z01-F003_03/raw_data/F_TI_4_3/F_TI_4_3_2.fq.gz"]}}

all_species_gDNA = ["To", "Tp", "TrR"]
all_species_RNA = ["To", "Tp", "88", "S10", "Tienshan"]

tissues = {"To": ["F_TO_1_1", "F_TO_1_2", "F_TO_1_3",
                  "F_TO_2_1", "F_TO_2_2", "F_TO_2_3",
                  "F_TO_3_1", "F_TO_3_2", "F_TO_3_3",
                  "F_TO_4_1", "F_TO_4_2", "F_TO_4_3"],
           "Tp": ["F_TP_1_1", "F_TP_1_2", "F_TP_1_3",
                  "F_TP_2_1", "F_TP_2_2", "F_TP_2_3",
                  "F_TP_3_1", "F_TP_3_2", "F_TP_3_3",
                  "F_TP_4_1", "F_TP_4_2", "F_TP_4_3"],
           "88": ["F_88_1_1", "F_88_1_2", "F_88_1_3",
                  "F_88_2_1", "F_88_2_2", "F_88_2_3",
                  "F_88_3_1", "F_88_3_2", "F_88_3_3",
                  "F_88_4_1", "F_88_4_2", "F_88_4_3"],
           "S10": ["F_S1_1_1", "F_S1_1_2", "F_S1_1_3",
                   "F_S1_2_1", "F_S1_2_2", "F_S1_2_3",
                   "F_S1_3_1", "F_S1_3_2", "F_S1_3_3",
                   "F_S1_4_1", "F_S1_4_2", "F_S1_4_3",
                   "A_S1_1_1", "A_S1_1_2", "A_S1_1_3",
                   "A_S1_2_1", "A_S1_2_2", "A_S1_2_3",
                   "A_S1_3_1", "A_S1_3_2", "A_S1_3_3",
                   "A_S1_4_1", "A_S1_4_2", "A_S1_4_3"],
           "Tienshan": ["F_TI_1_1", "F_TI_1_2", "F_TI_1_3",
                        "F_TI_2_1", "F_TI_2_2", "F_TI_2_3",
                        "F_TI_3_1", "F_TI_3_2", "F_TI_3_3",
                        "F_TI_4_1", "F_TI_4_2", "F_TI_4_3"]}


gwf.target_from_template("ToIndex",
                         bowtie2_index(reference_file, index_ref_file))

def star_mapping(read1, output, output_prefix, genomeDir):
    inputs = [read1, genomeDir+"/genomeParameters.txt"]
    outputs = [output]
    options = {
	"cores": 8,
	"memory": "64g",
	"account": "NChain",
	"walltime": "12:00:00"}
    OFilePrefix = output_prefix

    spec = """
    source /com/extra/STAR/2.5.2b/load.sh
    STAR --runThreadN 8 --genomeDir {dir} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {of} --readFilesIn {r1} --limitGenomeGenerateRAM=64000000000""".format(r1 = read1, of=output_prefix, dir=genomeDir)

    if read1.split(".")[-1]=="gz":
        spec += " --readFilesCommand zcat"

    return inputs, outputs, options, spec

def bowtie2_mapping(read1, output, genomeDir):
    inputs = [read1, genomeDir+".1.bt2"]
    outputs = [output]
    options = {
	"cores": 8,
	"memory": "24g",
	"account": "NChain",
	"walltime": "12:00:00"}

    spec = """
/home/marnit/miniconda3/envs/HyLiTE/bin/bowtie2 -N 1 -x {dir} -U {reads} -S {sam} --local -p {cores}
    """.format(reads = read1, sam=output, dir=genomeDir, cores=options["cores"])

    return inputs, outputs, options, spec


for species in all_species_gDNA:
    for i, readfile in enumerate(raw_gDNA[species]):
        gwf.target_from_template("Smpl"+species+"gDNAmap"+str(i+1),
                                 bowtie2_mapping(readfile, "./gDNA/{s}/{s}.{i}.gDNA.sam".format(s=species, i=i+1),
                                                 index_ref_file))

for species in all_species_RNA:
    for tissue in tissues[species]:
        for i, readfile in enumerate(raw_RNA[species][tissue]):
            gwf.target_from_template("Smpl"+species+"_"+tissue+"_"+"RNAmap"+str(i+1),
                                     bowtie2_mapping(readfile, "./RNA/{s}/{s}.{t}.{i}.gDNA.sam".format(s=species, t=tissue, i=i+1),
                                                     index_ref_file))

def samtools_merge(infile1, infile2, outfile):
    inputs = [infile1, infile2]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '4g',
        'account': 'NChain',
        'walltime': '24:00:00'
    }

    spec = """
    source /com/extra/samtools/1.3/load.sh
    samtools merge {} {} {}
    """.format(outfile, infile1, infile2)

    return inputs, outputs, options, spec


for species in all_species_gDNA:
    gwf.target_from_template("Smpl"+species+"gDNAmerge",
                             samtools_merge("./gDNA/{s}/{s}.{i}.gDNA.sam".format(s=species, i=1),
                                            "./gDNA/{s}/{s}.{i}.gDNA.sam".format(s=species, i=2),
                                            "./gDNA/{s}/{s}.gDNA.sam".format(s=species)))

for species in all_species_RNA:
    for tissue in tissues[species]:
        gwf.target_from_template("Smpl"+species+"_"+tissue+"_"+"RNAmerge",
                                 samtools_merge("./RNA/{s}/{s}.{t}.{i}.gDNA.sam".format(s=species, t=tissue, i=1),
                                                "./RNA/{s}/{s}.{t}.{i}.gDNA.sam".format(s=species, t=tissue, i=2),
                                                "./RNA/{s}/{s}.{t}.gDNA.sam".format(s=species, t=tissue)))


def samtools_process(infile, outfile):
    inputs = [infile]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = """
    source /com/extra/samtools/1.3/load.sh

    samtools sort {infile} -o {outfile}
    samtools index {outfile}
    """.format(infile=infile, outfile=outfile)

    return inputs, outputs, options, spec

for species in all_species_gDNA:
    gwf.target_from_template("Smpl"+species+"gDNAsort",
                             samtools_process("./gDNA/{s}/{s}.gDNA.sam".format(s=species),
                                            "./gDNA/{s}.gDNA.s.bam".format(s=species)))

for species in all_species_RNA:
    for tissue in tissues[species]:
        gwf.target_from_template("Smpl"+species+"_"+tissue+"_"+"RNAsort",
                                 samtools_process("./RNA/{s}/{s}.{t}.gDNA.sam".format(s=species, t=tissue),
                                                "./RNA/{s}.{t}.RNA.s.bam".format(s=species, t=tissue)))

def make_protocol_file(in_files, individual, samples, ploidy, filetype, outfile):
    inputs = in_files
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = ''
    for inf, ind, sample, p, f in zip(in_files, individual, samples, ploidy, filetype):
        spec += 'echo "'+"\t".join([ind, p, sample, f, inf])+'" >> '+outfile+"\n"

    return inputs, outputs, options, spec


Name = {"To": "P1", "Tp": "P2", "TrR": "Ch"}
Ploidy = {"To": "1", "Tp": "1", "TrR": "2"}

bamfiles = []
individual = []
ploidy_levels = []
seq_type = []

#### MAKE SAMPLES LIST
samples = []


translate = {"To": "To", "Tp": "Tp",
             "88": "TrR", "S10": "TrR",
             "Tienshan": "TrR"}

samples_map = {"To": [["F_TO_1_1", "F_TO_1_2", "F_TO_1_3",
                      "F_TO_2_1", "F_TO_2_2", "F_TO_2_3",
                       "F_TO_3_1", "F_TO_3_2", "F_TO_3_3",
                       "F_TO_4_1", "F_TO_4_2", "F_TO_4_3"]],
               "Tp": [["F_TP_1_1", "F_TP_1_2", "F_TP_1_3",
                       "F_TP_2_1", "F_TP_2_2", "F_TP_2_3",
                       "F_TP_3_1", "F_TP_3_2", "F_TP_3_3",
                       "F_TP_4_1", "F_TP_4_2", "F_TP_4_3"]],
               "S10": [["F_S1_1_1", "F_S1_1_2", "F_S1_1_3",
                       "F_S1_2_1", "F_S1_2_2", "F_S1_2_3",
                       "F_S1_3_1", "F_S1_3_2", "F_S1_3_3",
                       "F_S1_4_1", "F_S1_4_2", "F_S1_4_3",
                       "A_S1_1_1", "A_S1_1_2", "A_S1_1_3",
                       "A_S1_2_1", "A_S1_2_2", "A_S1_2_3",
                       "A_S1_3_1", "A_S1_3_2", "A_S1_3_3",
                       "A_S1_4_1", "A_S1_4_2", "A_S1_4_3"]],
               "88": [["F_88_1_1", "F_88_1_2", "F_88_1_3",
                      "F_88_2_1", "F_88_2_2", "F_88_2_3",
                      "F_88_3_1", "F_88_3_2", "F_88_3_3",
                      "F_88_4_1", "F_88_4_2", "F_88_4_3"]],
               "Tienshan": [["F_TI_1_1", "F_TI_1_2", "F_TI_1_3",
                            "F_TI_2_1", "F_TI_2_2", "F_TI_2_3",
                            "F_TI_3_1", "F_TI_3_2", "F_TI_3_3",
                            "F_TI_4_1", "F_TI_4_2", "F_TI_4_3"]]}

#species_included = {"To": False, "Tp": False, "TrR": False}

for species in [["88", "S10", "Tienshan"], "To", "Tp"]:
    if type(species)==list:
        sample_n = 1
        for specie in species:
            #print(specie, samples_map[specie])
            for current_samples in samples_map[specie]:
                for tissue in current_samples:
                    #print(specie, tissue)
                    bamfiles.append("./RNA/{s}.{t}.RNA.s.bam".format(s=specie, t=tissue))
                    seq_type.append("RNAseq")
                    individual.append(Name[translate[specie]])
                    ploidy_levels.append(Ploidy[translate[specie]])
                    samples.append("sample"+str(sample_n))
                    sample_n += 1
        for specie in [species[0]]:
            bamfiles.append("./gDNA/{s}.gDNA.s.bam".format(s=translate[specie]))
            seq_type.append("gDNA")
            individual.append(Name[translate[specie]])
            ploidy_levels.append(Ploidy[translate[specie]])
            samples.append("sample"+str(sample_n))
            sample_n += 1
    else:
        sample_n = 1
        for current_samples in samples_map[species]:
            for tissue in current_samples:
                bamfiles.append("./RNA/{s}.{t}.RNA.s.bam".format(s=species, t=tissue))
                seq_type.append("RNAseq")
                individual.append(Name[translate[species]])
                ploidy_levels.append(Ploidy[translate[species]])
                samples.append("sample"+str(sample_n))
                sample_n += 1
        bamfiles.append("./gDNA/{s}.gDNA.s.bam".format(s=translate[species]))
        seq_type.append("gDNA")
        individual.append(Name[translate[species]])
        ploidy_levels.append(Ploidy[translate[species]])
        samples.append("sample"+str(sample_n))


#print(bamfiles)
#print(individual)
#print(ploidy_levels)
#print(seq_type)

gwf.target_from_template("HyLiTEprotocolfile",
                         make_protocol_file(bamfiles,
                                            individual,
                                            samples,
                                            ploidy_levels,
                                            seq_type,
                                            "repens_protocol.txt"))

def make_bam_file(in_files,outfile):
    inputs = in_files
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '1g',
        'account': 'NChain',
        'walltime': '01:00:00'
    }

    spec = ""
    for ind in in_files:
        spec += 'echo "'+ind+'" >> '+outfile+"\n"

    return inputs, outputs, options, spec

gwf.target_from_template("mpileup_bamfile",
                         make_bam_file(bamfiles, "bamfiles.txt"))

def mpileup(reference, bamfiles, outfile):
    inputs = [reference, bamfiles]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '24g',
        'account': 'NChain',
        'walltime': '48:00:00'
    }

    spec = '''
    source /com/extra/samtools/1.3/load.sh

    samtools mpileup -BQ 0 -d 1000000 -f {ref} -b {bamfiles} > {outfile}
    '''.format(ref=reference, bamfiles=bamfiles, outfile=outfile)

    return inputs, outputs, options, spec


gwf.target_from_template("mpileup",
                         mpileup(reference_file, "bamfiles.txt",
                                 "repens_combo.pileup"))


def fix_mpileup(pileup_file, outfile):
    inputs = [pileup_file]
    outputs = [outfile]
    options = {
        'cores': 1,
        'memory': '2g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    spec = '''
    python mpileupfix.py {pileup} > {out}
    '''.format(pileup=pileup_file, out=outfile)

    return inputs, outputs, options, spec

gwf.target_from_template("mpileupfix",
                         fix_mpileup("repens_combo.pileup", "repens_fixed.pileup"))


def run_HyLiTE(reference, protocol, pileup_file):
    inputs = [reference, protocol, pileup_file]
    outputs = ["HyLiTE_output/HyLiTE_output.expression.txt"]
    options = {
        'cores': 1,
        'memory': '124g',
        'account': 'NChain',
        'walltime': '64:00:00'
    }

    spec = '''
/home/marnit/miniconda3/envs/HyLiTE/bin/HyLiTE -v -f {proto} -p {pileup} -n HyLiTE_output -r {ref}
    '''.format(ref=reference, proto=protocol, pileup=pileup_file)

    return inputs, outputs, options, spec


def HacknHyLiTE(reference, protocol, pileup_file, segments):
    inputs = [reference, protocol, pileup_file]
    outputs = []
    options = {
        'cores': 1,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '12:00:00'
    }

    for i in range(segments):
        outputs.append("./slicing_directories/slicing.subset{}.sh".format(i))

    spec = '''
 /home/marnit/miniconda3/envs/HyLiTE/bin/HacknHyLiTE -n {seg} -o slicing_directories --name slicing -f {proto} -p {pileup} --options "-r {ref}"
    '''.format(seg=segments, ref=reference, proto=protocol, pileup=pileup_file)

    return inputs, outputs, options, spec


#gwf.target_from_template("HyLiTE",
#                         run_HyLiTE(reference_file, "repens_protocol.txt", "repens_fixed.pileup"))


n_segments = 50

gwf.target_from_template("HacknHyLiTE",
                         HacknHyLiTE(reference_file, "repens_protocol.txt", "repens_fixed.pileup",
                                     n_segments))

subsets = []
subsets_out = []
for i in range(n_segments):
    subsets.append("./slicing_directories/slicing.subset{}.sh".format(i))
    subsets_out.append("./slicing_directories/subset{}/slicing.subset{}.expression.txt".format(i, i))


def run_subset(subset, output):
    inputs = [subset]
    outputs = [output]
    options = {
        'cores': 1,
        'memory': '48g',
        'account': 'NChain',
        'walltime': '48:00:00'
    }

    spec = '''
    export PATH=$PATH:/home/marnit/miniconda3/envs/HyLiTE/bin/

    bash {subset}
    '''.format(subset=subset)

    return inputs, outputs, options, spec


for i, subset in enumerate(zip(subsets, subsets_out)):
    gwf.target_from_template("Subset"+str(i),
                             run_subset(subset[0], subset[1]))

def mergeHyLiTE(subset_out, directory, reference):
    inputs = subset_out
    outputs = []
    options = {
        'cores': 1,
        'memory': '12g',
        'account': 'NChain',
        'walltime': '48:00:00'
    }

    spec = '''
    cd {dir}

    /home/marnit/miniconda3/envs/HyLiTE/bin/HyLiTE-merge -r {ref}
    '''.format(dir=directory, ref=reference)

    return inputs, outputs, options, spec


gwf.target_from_template("mergeHyLiTE",
                         mergeHyLiTE(subsets_out, "./slicing_directories/", "."+reference_file))
