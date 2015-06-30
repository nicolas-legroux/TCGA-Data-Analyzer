ROOT_DIR = "/home/nlegroux/TCGA-Data-Analyzer/"
EDGES = ROOT_DIR + "data/graph/biogrid-edges.txt"
NEGATIVE_WEIGHTS = ["0.5", "1.5", "2.5", "3.5", "4.5", "5.5"]
NEGATIVE_WEIGHTS_STRING = ",".join(NEGATIVE_WEIGHTS)
MAX_CONTROL = 100
MAX_TUMOR = 700
CANCERS = "BRCA,COAD,GBM,HNSC,KIRC,LGG,LUAD,LUSC,OV,PRAD,THCA,UCEC"
SAMPLES = []

try:
	SAMPLES = open(ROOT_DIR + "data/heinz/samples.list", "r").readlines()
	SAMPLES = map(str.strip, SAMPLES)
except (FileNotFoundError, IOError):
	print("samples.list was not found in the initialization. Please run 'snakemake init' before running 'snakemake'.")
	
rule all:
	input: ROOT_DIR + "data/heinz/samples.list",  \
expand(ROOT_DIR + "data/heinz/output/{weight}_{sample}.stats", weight=NEGATIVE_WEIGHTS, sample=SAMPLES), \
expand(ROOT_DIR + "data/heinz/output/{weight}_{sample}.nodes", weight=NEGATIVE_WEIGHTS, sample=SAMPLES)

rule graph:
	input: ROOT_DIR + "data/heinz/raw_output/{weight}_{sample}.txt"
	output: ROOT_DIR + "data/heinz/output/{weight}_{sample}.png"
	shell: "dot {input} -T png > {output}"

rule heinz_stats:
	input: ROOT_DIR + "data/heinz/raw_output/{weight}_{sample}.txt", ROOT_DIR + "data/heinz/output/{weight}_{sample}.png"
	output: ROOT_DIR + "data/heinz/output/{weight}_{sample}.stats", ROOT_DIR + "data/heinz/output/{weight}_{sample}.nodes"
	params: w="{weight}", s="{sample}"
	shell: "{ROOT_DIR}TCGA-Analyzer -mode 2 -f {params.w}_{params.s} >/dev/null"

rule heinz:
	input: ROOT_DIR + "data/heinz/input/{weight}_{sample}.txt"
	output: ROOT_DIR + "data/heinz/raw_output/{weight}_{sample}.txt"
	threads: 1
	shell: "/home/thume/heinz/build/heinz -e {EDGES} -n {input} >{output} 2>/dev/null"

rule init:
	output: ROOT_DIR + "data/heinz/samples.list", ROOT_DIR + "data/heinz/negative-weights.txt"
	shell: "{ROOT_DIR}TCGA-Analyzer -mode 1 -cancers {CANCERS} -weights {NEGATIVE_WEIGHTS_STRING} -maxcontrol {MAX_CONTROL} -maxtumor {MAX_TUMOR}"

rule clean:
	shell: "rm -f {ROOT_DIR}data/heinz/samples.list; rm -f {ROOT_DIR}data/heinz/input/*; rm -f {ROOT_DIR}data/heinz/output/*; rm -f {ROOT_DIR}data/heinz/raw_output/*;rm -f {ROOT_DIR}data/heinz/negative-weights.txt; rm -f {ROOT_DIR}export/*"
