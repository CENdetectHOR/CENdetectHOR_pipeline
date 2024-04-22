import os

#FASTA_FOL = os.path.abspath(config["fasta-folder"])
#FASTA_FOL = os.path.abspath("~/CENdetectHOR")

#chrom=[15]
#prefix="HSA"
CONS="human_cons_alphasat.txt"
#chrom=[10,15]
#prefix="HSA"

prefix, chrom= glob_wildcards("fasta/{PREFIX}.chr{CHR}.fasta")

rule all:
	input:
		tree=expand("/lustrehome/alessiadaponte/CENdetectHOR/results/HOR/{PREFIX}.FULLchr.tree.xml", PREFIX=prefix),
		hist=expand("/lustrehome/alessiadaponte/CENdetectHOR/results/HOR/plots/{PREFIX}.dist_hist.pdf", PREFIX=prefix),
                matr=expand("/lustrehome/alessiadaponte/CENdetectHOR/results/HOR/plots/{PREFIX}.dist_matrix.pdf", PREFIX=prefix),
		fullWind=expand("/lustrehome/alessiadaponte/CENdetectHOR/results/wind2analize/{PREFIX}.FULLchr.windows.filtered.bed", PREFIX=prefix),
		fullMons=expand("/lustrehome/alessiadaponte/CENdetectHOR/results/monomers/{PREFIX}.FULLchr_mons.fasta", PREFIX=prefix),
		fullSeq=expand("/lustrehome/alessiadaponte/CENdetectHOR/fasta/{PREFIX}.FULLchr.fasta", PREFIX=prefix),
		rawF=expand("/lustrehome/alessiadaponte/CENdetectHOR/results/decomposition/{PREFIX}/final_decomposition_raw.tsv", PREFIX=prefix)
		#consExt=expand("results/monomers/{PREFIX}_chr0_cons.fasta", PREFIX=prefix)

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
		seq="/lustrehome/alessiadaponte/CENdetectHOR/fasta/{PREFIX}.chr{CHR}.fasta"
	output:
		windSumm="/lustrehome/alessiadaponte/CENdetectHOR/results/wind_summary/{PREFIX}.chr{CHR}.txt",
		bed="/lustrehome/alessiadaponte/CENdetectHOR/results/wind2analize/{PREFIX}.chr{CHR}.windows.bed"
	params:
		plotFold="/lustrehome/alessiadaponte/CENdetectHOR/results/plots/chr{CHR}/"
	shell:
		'''
		#mkdir {params.plotFold}
		python ./scripts/periodicityScript_full_fin.py --file {input.seq} --plot-fold {params.plotFold} --wind-summ {output.windSumm} --wind-bed {output.bed}
		'''

rule selectWind:
	input:
		windF=rules.periodicity.output.bed
	output:
		selWind="/lustrehome/alessiadaponte/CENdetectHOR/results/wind2analize/{PREFIX}.chr{CHR}.windows.filtered.bed"
	shell:
		'''
		num=$(cut -f 4 {input.windF} |sort|uniq |head -n 1); awk -v n="$num" '$4 == n' {input.windF} > {output.selWind}
		'''

rule extractCons:
	input:
		seq="/lustrehome/alessiadaponte/CENdetectHOR/fasta/{PREFIX}.chr{CHR}.fasta",
		wind=rules.selectWind.output.selWind
	output:
		consExt="/lustrehome/alessiadaponte/CENdetectHOR/results/monomers/{PREFIX}_chr0cons.fasta"
	shell:
		'''
		python ./scripts/Extract5mon_NOcons.py --fasta {input.seq} --bed {input.wind} --cons-file {output.consExt}		
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
		seq="/lustrehome/alessiadaponte/CENdetectHOR/fasta/{PREFIX}.chr{CHR}.fasta",
		wind=rules.selectWind.output.selWind,
		cons=CONS if os.path.isfile(CONS) else rules.extractCons.output.consExt
	output:
		mon="/lustrehome/alessiadaponte/CENdetectHOR/results/monomers/{PREFIX}.chr{CHR}_mons.fasta"	
	shell:
		'''
		python ./scripts/Extract5mon.py --fasta {input.seq} --bed {input.wind} --cons-file {input.cons} --out {output.mon}
		'''

rule concatenateMons:
	input:
		mon=expand("/lustrehome/alessiadaponte/CENdetectHOR/results/monomers/{PREFIX}.chr{CHR}_mons.fasta", CHR=chrom,PREFIX=prefix),
		cons=CONS if os.path.isfile(CONS) else rules.extractCons.output.consExt
	output:
		fullMons="/lustrehome/alessiadaponte/CENdetectHOR/results/monomers/{PREFIX}.FULLchr_mons.fasta"
	shell:
		'''
		cat {input.mon} {input.cons} > {output.fullMons}
		'''

rule concatenateSeq:
	input:
		seq2conc=expand("/lustrehome/alessiadaponte/CENdetectHOR/fasta/{PREFIX}.chr{CHR}.fasta", CHR=chrom, allow_missing=True)
	output:
		fullSeq="/lustrehome/alessiadaponte/CENdetectHOR/fasta/{PREFIX}.FULLchr.fasta"
	shell:
		'''
		cat {input.seq2conc} > {output.fullSeq}
		'''

rule concatenateWind:
        input:
                wind2conc=expand("/lustrehome/alessiadaponte/CENdetectHOR/results/wind2analize/{PREFIX}.chr{CHR}.windows.filtered.bed", CHR=chrom, allow_missing=True)
        output:
                fullWind="/lustrehome/alessiadaponte/CENdetectHOR/results/wind2analize/{PREFIX}.FULLchr.windows.filtered.bed"
        shell:
                '''
                cat {input.wind2conc} > {output.fullWind}
                '''

rule stringDec:
	input:
		fullSeq=rules.concatenateSeq.output.fullSeq,
		fullMons=rules.concatenateMons.output.fullMons
	params:
		outdir="/lustrehome/alessiadaponte/CENdetectHOR/results/decomposition/{PREFIX}"
	output:
		rawF="/lustrehome/alessiadaponte/CENdetectHOR/results/decomposition/{PREFIX}/final_decomposition_raw.tsv"
	threads: 16
	shell:
		'''
		stringdecomposer {input.fullSeq} {input.fullMons} -o {params.outdir} -t {threads}
		'''

rule mon2bed:
	input:
		stringDec=rules.stringDec.output.rawF,
		wind_bed=rules.concatenateWind.output.fullWind
	output:
		decBed="/lustrehome/alessiadaponte/CENdetectHOR/results/decomposition/{PREFIX}_final_decomposition.bed"
	shell:
		'''
		Rscript ./scripts/dec2bed_fin_fullChr.R {input.stringDec} {input.wind_bed} {output.decBed}
		'''

rule HORdet:
	input:
		fullSeq=rules.concatenateSeq.output.fullSeq,
		fullStringDec=rules.mon2bed.output.decBed
	#params:
		#plotdir="/lustrehome/alessiadaponte/CENdetectHOR/results/HOR/plots"
	output:
		tree="/lustrehome/alessiadaponte/CENdetectHOR/results/HOR/{PREFIX}.FULLchr.tree.xml",
		hist="/lustrehome/alessiadaponte/CENdetectHOR/results/HOR/plots/{PREFIX}.dist_hist.pdf",
		matr="/lustrehome/alessiadaponte/CENdetectHOR/results/HOR/plots/{PREFIX}.dist_matrix.pdf"
	log:
		"/lustrehome/alessiadaponte/CENdetectHOR/results/{PREFIX}_HOR_clustering.log"
	shell:
		'''
		python ./scripts/findHORsFromMonomersAD.py --fasta {input.fullSeq} --bed {input.fullStringDec} --out-tree {output.tree} --out-distHist {output.hist} --out-distMatr {output.matr} --log {log}
		'''

