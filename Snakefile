import os

configfile: "config/config.yaml"
workdir: config["workdir"]


CONS= config["consFile"]
Wind=int(config["windowSize"])
Kmer=int(config["kmerSize"])
maxGap=int(config["maxHORgap"])
maxHOR=int(config["maxHORlength"])
minrep=int(config["minHORrep"])
discTree=config["discreteTree"]
HOR_t=int(config["HORdet_threads"])

prefix, chrom= glob_wildcards("fasta/{PREFIX}.chr{CHR}.fasta")

rule all:
	input:
		tree=expand("results/HOR/{PREFIX}.{CHR}.tree.xml", PREFIX=prefix, CHR=chrom),
		rawF=expand("results/decomposition/{PREFIX}/{CHR}/final_decomposition_raw.tsv", PREFIX=prefix, CHR=chrom)


rule periodicity:
	input:
		seq="fasta/{PREFIX}.chr{CHR}.fasta"
	output:
		windSumm="results/wind_summary/{PREFIX}.chr{CHR}.txt",
		bed="results/wind2analize/{PREFIX}.chr{CHR}.windows.bed"
	params:
		plotFold="results/plots/{PREFIX}_chr{CHR}/",
		script=os.path.join(workflow.basedir, "scripts/periodicityScript_full_fin.py"),
		W=Wind,
		K=Kmer
	shell:
		'''
		python {params.script} --fasta {input.seq} --plot-fold {params.plotFold} --wind-summ {output.windSumm} --wind-bed {output.bed} --window-size {params.W} --kmer-size {params.K}  
		'''

rule concatenateWind:
        input:
                wind2conc=expand("results/wind2analize/{PREFIX}.chr{CHR}.windows.bed", CHR=chrom, allow_missing=True)
        output:
                fullWind="results/wind2analize/{PREFIX}.FULLchr.windows.bed"
        shell:
                '''
                cat {input.wind2conc} > {output.fullWind}
                '''

rule selectWind:
	input:
		windF=rules.concatenateWind.output.fullWind
	output:
		selWind="results/wind2analize/{PREFIX}.FULLchr.windows.filtered.bed"
	params:
		outfold="results/wind2analize/filtering/",
		script=os.path.join(workflow.basedir, "scripts/windowsFiltering.py")
	shell:
		'''
		python {params.script} --bed {input.windF} --out {params.outfold}
		cp {params.outfold}*_main.bed {output.selWind}
		'''

rule extractCons:
	input:
		seq="fasta/{PREFIX}.chr{CHR}.fasta",
		wind=rules.selectWind.output.selWind
	output:
		consExt="monomers/{PREFIX}_chr0cons.fasta"
	params:
		script=os.path.join(workflow.basedir, "scripts/Extract5mon_NOcons.py")
	shell:
		'''
		python {params.script} --fasta {input.seq} --bed {input.wind} --out-cons-file {output.consExt}		
		'''

rule extractMon:
	input:
		seq="fasta/{PREFIX}.chr{CHR}.fasta",
		wind=rules.selectWind.output.selWind,
		cons=CONS if os.path.isfile(CONS) else rules.extractCons.output.consExt
	output:
		mon="results/monomers/{PREFIX}.chr{CHR}_mons.fasta"	
	params:
		script=os.path.join(workflow.basedir, "scripts/Extract5mon.py")
	shell:
		'''
		python {params.script} --fasta {input.seq} --bed {input.wind} --cons-file {input.cons} --out {output.mon}
		'''

rule concatenateMons:
	input:
		mon=expand("results/monomers/{PREFIX}.chr{CHR}_mons.fasta", CHR=chrom,PREFIX=prefix),
		cons=CONS if os.path.isfile(CONS) else rules.extractCons.output.consExt
	output:
		fullMons="results/monomers/{PREFIX}.FULLchr_mons.fasta"
	shell:
		'''
		cat {input.mon} {input.cons} > {output.fullMons}
		'''

rule stringDec:
	input:
		seq="fasta/{PREFIX}.chr{CHR}.fasta",
		fullMons=rules.concatenateMons.output.fullMons
	params:
		outdir="results/decomposition/{PREFIX}/{CHR}"
	output:
		rawF="results/decomposition/{PREFIX}/{CHR}/final_decomposition_raw.tsv"
	threads: 4
	shell:
		'''
		stringdecomposer {input.seq} {input.fullMons} -o {params.outdir} -t {threads}
		'''

rule mon2bed:
	input:
		stringDec=rules.stringDec.output.rawF,
		wind_bed=rules.selectWind.output.selWind
	output:
		decBed="results/decomposition/{PREFIX}/{CHR}/final_decomposition.bed"
	params:
		script=os.path.join(workflow.basedir, "scripts/dec2bed_fin.R")
	shell:
		'''
		Rscript {params.script} {input.stringDec} {input.wind_bed} {output.decBed}
		'''

rule HORdet:
	input:
		seq="fasta/{PREFIX}.chr{CHR}.fasta",
		fullStringDec=rules.mon2bed.output.decBed
	params:
		matr="results/HOR/distMatr/chr{CHR}/",
		script=os.path.join(workflow.basedir, "scripts/findHORsFromMonomers_new.py"),
		MG=maxGap,
		MH=maxHOR,
		mL=minrep,
		DL=discTree
	output:
		tree="results/HOR/{PREFIX}.{CHR}.tree.xml"
	log:
		"results/{PREFIX}_{CHR}_HOR_clustering.log"
	threads: HOR_t
	shell:
		'''
		python {params.script} --fasta {input.seq} --bed {input.fullStringDec} --out-tree {output.tree} --out-distMatr {params.matr} --log {log} --t {threads} --max-gap {params.MG} --max-HOR {params.MH} --min-loops {params.mL} --dis-lev {params.DL}
		'''
