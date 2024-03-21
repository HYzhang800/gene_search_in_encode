configfile: "config.yml"
workdir: config['workdir']
import os
import pandas as pd
import polars as pl
from gtfparse import read_gtf
from os import path
import glob

SNAKEDIR = path.dirname(workflow.snakefile)


encode_dir = config['gtf_dir']

if not encode_dir.endswith("/"):
    encode_dir += "/"

files = os.listdir(encode_dir)

#encode_array = [file for file in files if file.endswith(".gtf")]
encode_array = [os.path.join(encode_dir, file) for file in files if file.endswith(".gtf")]

print("Files to be processed:", encode_array)


          
gene_coord_array = []

ref = config["reference_annotation"]
gene_file = config["target_genes"]

# Function to extract gene coordinates
def extract_gene_coordinates(reference_annotation, target_genes):
    with open(target_genes, 'r') as file:
       genes= [line.strip() for line in file.readlines()]

    ref_gencode = read_gtf(reference_annotation)
    ref_gencode = ref_gencode.to_pandas()
    ref_gencode = ref_gencode[ref_gencode["feature"] == "gene"]
    ref_gencode = ref_gencode[ref_gencode["gene_name"].isin(genes)]
    gffread_input_list = ref_gencode["strand"].astype(str) + ref_gencode["seqname"].astype(str) + ":" + ref_gencode["start"].astype(str) + ".." + ref_gencode["end"].astype(str)
    os.makedirs("Results", exist_ok=True)
    ref_gencode.to_csv("Results/gencode_genes.csv", index=False)
    return gffread_input_list.tolist(), genes


# Run the function to extract gene coordinates
gene_coord_array, genes = extract_gene_coordinates(ref, gene_file)

print("Gene coordinates extracted:", gene_coord_array)


rule all:
    input:
          #"Results/top_transcripts_sum.png",
          expand("Merged_gtf/merged_{gtf_name}", gtf_name = encode_array),
          trans_info = expand("Results/transcript_usage_{gene}.csv", gene = genes)

for p in encode_array:
  for g in gene_coord_array:

# Select transcripts from annotation files
   rule:
     name: f"gffread_select_{p}_for_{g}" 
     input:
        encode_gtf = f"{p}"
        
     output:
        selected_gtf = f"processing_gtf/selected_{g}_{p}"  
     params:        
        gene_coord = f"{g}"
     shell:
          """
          gffread -r {params.gene_coord} {input.encode_gtf} -T -o {output.selected_gtf}
          """


rule merge_gtf:
    input:
        expand("processing_gtf/selected_{g}_{p}", g=gene_coord_array, p=encode_array)
    output:
        #f"Merged_gtf/merged_{p}"
        merged_gtf = expand("Merged_gtf/merged_{p}", p=encode_array)
    run:
        # Iterate over each p
          for p in encode_array:
            # Get the list of GTF files containing the current p in their filename

            gtf_files = glob.glob('processing_gtf/selected*' + str(p))
            print("gtf_files for ", str(p), ": ", gtf_files)

            # Create the merged GTF file path

            #merged_gtf_file = output["merged_gtf"].join(f"Merged_gtf/merged_{p}.gtf")
            merged_gtf_file = ("Merged_gtf/merged_" + str(p))
            # Concatenate the GTF files into the merged GTF file
            with open(merged_gtf_file, "w") as merged_file:
               for gtf_file in gtf_files:
                  with open(gtf_file, "r") as individual_file:
                     merged_file.write(individual_file.read())


rule generate_file_list:
    input:
         merged_files = rules.merge_gtf.output.merged_gtf
    output:
         file_list = "Results/filtered_gtf_list.txt"
    run:
         with open(output.file_list, "w") as textfile:
             for filename in os.listdir("Merged_gtf/"):
                if filename.startswith("merged"):
                   textfile.write(filename + "\n")
          
 


rule data_process:
    input:
        selected_gtf_dir = rules.generate_file_list.output.file_list,
        gene_info = "Results/gencode_genes.csv"
    output:
        #top_transcripts = "Results/top_transcripts_sum.png",
        #dou = "Results/dou.png",
        #orf_usage_plot = "Results/orf_usage.png",
        #trans_info = "Results/transcript_usage.csv",
        #trans_sum = "Results/transcript_summary.csv",
        #gene_quanti = "Results/gene_quanti.csv",
        #orf_usage = "Results/orf_usage_sum.csv"
        trans_info = expand("Results/transcript_usage_{gene}.csv", gene = genes)
    params:
        #script = "/mnt/c/Users/Haoyu/Desktop/scripts/gene_test_ENCODE.R",
        script = SNAKEDIR + "/Scripts/gene_test_ENCODE.R",
        metadata = config["metadata"],
        quanti_file_path = config["quanti_file_path"],
        reference_genome = config["reference_genome"],
        reference_annotation = config["reference_annotation"],
        #organ_info = "/mnt/e/Encode_annotation/all/organ_info_full.csv",
        organ_info = SNAKEDIR + "/Data/organ_info_full.csv",
        workdir = config["workdir"]
    log: "Results/Rlog/Rlog.log"
    shell:
       """
       Rscript {params.script} {input.selected_gtf_dir} {input.gene_info} {params.metadata} {params.quanti_file_path} {params.reference_genome} {params.reference_annotation} {params.organ_info} {params.workdir}
       """


###
