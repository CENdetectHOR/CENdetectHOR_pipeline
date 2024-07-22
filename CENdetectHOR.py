import os

#FASTA_FOL = os.path.abspath(config["fasta-folder"])
#FASTA_FOL = os.path.abspath("~/CENdetectHOR")

workdir: "/lustrehome/alessiadaponte/CENdetectHOR/"
#chrom=[15]
#prefix="HSA"
CONS="human_cons_alphasat.txt"
#chrom=[10,15]
#prefix="HSA"


prefix, chrom= glob_wildcards("fasta/{PREFIX}.chr{CHR}.fasta")

rule all:
	input:
		tree=expand("results/HOR/{PREFIX}.{CHR}.tree.xml", PREFIX=prefix, CHR=chrom),
		#disMatr=expand("/lustre/home/alessiadaponte/h2/results/HOR/distMatr/{PREFIX}.{CHR}/{CHR}.dist_matrix.npy", PREFIX=prefix, CHR=chrom),
		#hist=expand("/lustrehome/alessiadaponte/centromeres_fromGlennis/chimp/h1/results/HOR/plots/{PREFIX}.dist_hist.pdf", PREFIX=prefix),
                #matr=expand("/lustrehome/alessiadaponte/centromeres_fromGlennis/chimp/h1/results/HOR/plots/{PREFIX}.dist_matrix.pdf", PREFIX=prefix),
		#fullWind=expand("/lustrehome/alessiadaponte/centromeres_fromGlennis/chimp/h1/results/wind2analize/{PREFIX}.FULLchr.windows.filtered.bed", PREFIX=prefix),
		#fullMons=expand("/lustre/home/alessiadaponte/h2/results/monomers/{PREFIX}.FULLchr_mons.fasta", PREFIX=prefix),
		#seq=expand("/lustre/home/alessiadaponte/chimp/h2/fasta/{PREFIX}.chr{CHR}.fasta", PREFIX=prefix, CHR=chrom),
		rawF=expand("results/decomposition/{PREFIX}/{CHR}/final_decomposition_raw.tsv", PREFIX=prefix, CHR=chrom)
		#consExt=expand("results/monomers/{PREFIX}_c0_cons.fasta", PREFIX=prefix)

#rule stainedGlass:
#	input:
#		fasta=FASTA,
#		sample=""
#	threads: 8
#	shell:
#		'''
#		snakemake -use-conda -cores {threads}  -config fasta={input.fasta} sample={input.sample} make_figures
#		'''

rule periodicity:
	input:
		seq="fasta/{PREFIX}.chr{CHR}.fasta"
	output:
		windSumm="results/wind_summary/{PREFIX}.chr{CHR}.txt",
		bed="results/wind2analize/{PREFIX}.chr{CHR}.windows.bed"
	params:
		plotFold="results/plots/chr{CHR}/"
	shell:
		'''
		#mkdir {params.plotFold}
		python scripts/periodicityScript_full_fin.py --file {input.seq} --plot-fold {params.plotFold} --wind-summ {output.windSumm} --wind-bed {output.bed}
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
		outfold="results/wind2analize/filtering/"
	shell:
		'''
		python windowsFiltering.py --bed {input.windF} --out {params.outfold}
		cp {params.outfold}/*_main.bed {output.selWind}
		'''

rule extractCons:
	input:
		seq="fasta/{PREFIX}.chr{CHR}.fasta",
		wind=rules.selectWind.output.selWind
	output:
		consExt="monomers/{PREFIX}_chr0cons.fasta"
	shell:
		'''
		python scripts/Extract5mon_NOcons.py --fasta {input.seq} --bed {input.wind} --cons-file {output.consExt}		
		'''

#rule copyCons:
#	input:
#		consF=CONS
#	output:
#		consNewLoc="results/monomers/{prefix}.{CHR}mons.fasta"
#	shell:
#		'''
#		cp {input.consF} {output.consNewLoc}
#		'''

rule extractMon:
	input:
		seq="fasta/{PREFIX}.chr{CHR}.fasta",
		wind=rules.selectWind.output.selWind,
		cons=CONS if os.path.isfile(CONS) else rules.extractCons.output.consExt
	output:
		mon="results/monomers/{PREFIX}.chr{CHR}_mons.fasta"	
	shell:
		'''
		python scripts/Extract5mon.py --fasta {input.seq} --bed {input.wind} --cons-file {input.cons} --out {output.mon}
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

#rule concatenateSeq:
#	input:
#		seq2conc=expand("/lustre/home/alessiadaponte/centromeres_fromGlennis/chimp/h2/fasta/{PREFIX}.chr{CHR}.fasta", CHR=chrom, allow_missing=True)
#	output:
#		fullSeq="/lustre/home/alessiadaponte/centromeres_fromGlennis/chimp/h2/fasta/{PREFIX}.FULLchr.fasta"
#	shell:
#		'''
#		cat {input.seq2conc} > {output.fullSeq}
#		'''

#rule concatenateWind:
#        input:
#                wind2conc=expand("/lustre/home/alessiadaponte/centromeres_fromGlennis/chimp/h1/results/wind2analize/{PREFIX}.chr{CHR}.windows.filtered.bed", CHR=chrom, allow_missing=True)
#        output:
#                fullWind="/lustrehome/alessiadaponte/centromeres_fromGlennis/chimp/h1/results/wind2analize/{PREFIX}.FULLchr.windows.filtered.bed"
#        shell:
#                '''
#                cat {input.wind2conc} > {output.fullWind}
#                '''

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
	shell:
		'''
		Rscript scripts/dec2bed_fin.R {input.stringDec} {input.wind_bed} {output.decBed}
		'''

rule HORdet:
	input:
		seq="fasta/{PREFIX}.chr{CHR}.fasta",
		fullStringDec=rules.mon2bed.output.decBed
	params:
		matr="results/HOR/distMatr/chr{CHR}/"
	output:
		tree="results/HOR/{PREFIX}.{CHR}.tree.xml"
		#hist="/lustrehome/alessiadaponte/centromeres_fromGlennis/chimp/h1/results/HOR/plots/{PREFIX}.dist_hist.pdf",
		#disMatr="/lustre/home/alessiadaponte/chimp/h2/results/HOR/distMatr/{PREFIX}.{CHR}/{CHR}.dist_matrix.npy"
	log:
		"results/{PREFIX}_{CHR}_HOR_clustering.log"
	threads: 32
	shell:
		'''
		python scripts/findHORsFromMonomers_FullCen.py --fasta {input.seq} --bed {input.fullStringDec} --out-tree {output.tree} --out-distMatr {params.matr} --log {log} --t {threads}
		'''
