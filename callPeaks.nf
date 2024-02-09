#!/usr/bin/env nextflow
/*
* Name:         callPeaks.nf
*
* Description:  Call peaks for ChIP-seq data, tries to use a methodology reltivly close to the
*               already established REMAP pipeline.
*
* Input:        A reference genome to align to in FASTA format and one or more directories containing
*               raw ChIP-seq reads in FASTQ format.
*
* Output:       Two results per directory, one with narrow peaks and one with broad peaks, same file layout
*               as a normal MACS2 run
*/

nextflow.enable.dsl = 2

// Genome and directories containing reads
params.genomeFasta = "genome.fa"

params.chipSeqFastq = ""
params.chipSeqFastqList = params.chipSeqFastq?.split(',') as List

params.controlFastq = null

// changeable parameters
params.threads = 8
params.outputDir = "output"

// Run in paired end mode, instead of treating each file separately
params.paired = false
// Use the deprecated tool 'samtools rmdup' instead of the newer tools to remove duplicates
params.useRmdup = false
// Multithread samtools, but beware this may use a lot of RAM
params.threadSamtools = false

/**
* Build an index using bowtie2 for the reference genome, used in downstream processes
*/
process buildBowtieIndex {
    conda 'bowtie2'
    cpus params.threads

    input:
        path genomeFile

    output:
        path 'genomeIndex', emit: bowtieIndex

    shell:
    '''
    mkdir genomeIndex
    bowtie2-build --threads !{task.cpus} !{genomeFile} genomeIndex/index
    '''
}

/**
* Trim the raw ChIP-seq reads using trim galore, relying on autodetecting the settings 
*/
process trimChipReplicates {
    conda 'trim-galore'
    cpus params.threads

    input:
        tuple file(readFile), val(readName)

    output:
        path '*trimmed_reads', emit: trimmedReadsDir

    shell:
    if (params.paired) {
        '''
        trim_galore --output_dir !{readName}trimmed_reads --length 30 --quality 20 --stringency 1 -e 0.1 !{readFile}/!{readName}* --paired
        '''
    } else {
        '''
        trim_galore --output_dir !{readName}trimmed_reads --length 30 --quality 20 --stringency 1 -e 0.1 !{readFile}/!{readName}*
        '''
    }
}

/**
* Align the reads against the genome using bowtie2
*/
process alignToGenome {
    conda 'bowtie2'
    cpus params.threads
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path bowtieIndex
        path trimmedReadsDir

    output:
        path "*alignedReads.sam", emit: alignedReads

    shell:
    if (params.paired) {
        '''
        bowtie2 --threads !{task.cpus} --end-to-end --sensitive -x !{bowtieIndex}/index -1 !{trimmedReadsDir}/*1_val_1.fq -2 !{trimmedReadsDir}/*2_val_2.fq -S !{trimmedReadsDir}alignedReads.sam
        '''
    } else {
        '''
        bowtie2 --threads !{task.cpus} --end-to-end --sensitive -x !{bowtieIndex}/index -U !{trimmedReadsDir}/*.fq -S !{trimmedReadsDir}alignedReads.sam
        '''
    }
}

/**
* Use samtools to remove PCR duplicates, may use the old rmdup or new markdup -r, depending on parameters
*/
process removePCRDuplicates {
    conda 'samtools'
    cpus params.threads

    input:
        path inputSam

    output:
        path 'deduplicated.sam', emit: outputSam

    shell:
    threadCommand = params.threadSamtools ? "-@ ${task.cpus}" : ''

    if (params.useRmdup) {
        if (params.paired) {
            '''
            samtools rmdup !{inputSam} deduplicated.sam
            '''
        } else {
            '''
            samtools rmdup -s !{inputSam} deduplicated.sam
            '''
        }
    } else {
        '''
        samtools sort -n !{inputSam} -o sorted.sam !{threadCommand}
        samtools fixmate -m sorted.sam scored.sam !{threadCommand}
        samtools sort scored.sam -o scored_sorted.sam !{threadCommand}
        samtools markdup -r scored_sorted.sam deduplicated.sam !{threadCommand}
        '''
    }
}

/**
* Trim the control raw ChIP-seq reads using trim galore, relying on autodetecting the settings 
*/
process trimControl {
    conda 'trim-galore'
    cpus params.threads

    input:
        file controlFiles

    output:
        path 'trimmed_reads', emit: trimmedReadsDir

    shell:
    if (params.paired) {
        '''
        trim_galore --output_dir trimmed_reads --length 30 --quality 20 --stringency 1 -e 0.1 !{controlFiles}/* --paired
        '''
    } else {
        '''
        trim_galore --output_dir trimmed_reads --length 30 --quality 20 --stringency 1 -e 0.1 !{controlFiles}/*
        '''
    }
}

/**
* Align the control reads against the genome using bowtie2
*/
process alignControlToGenome {
    conda 'bowtie2'
    cpus params.threads

    input:
        path bowtieIndex
        path trimmedReadsDir

    output:
        path 'alignedReads.sam', emit: alignedReads

    shell:
    if (params.paired) {
        '''
        bowtie2 --threads !{task.cpus} --end-to-end --sensitive -x !{bowtieIndex}/index -1 !{trimmedReadsDir}/*1_val_1.fq -2 !{trimmedReadsDir}/*2_val_2.fq -S alignedReads.sam
        '''
    } else {
        '''
        bowtie2 --threads !{task.cpus} --end-to-end --sensitive -x !{bowtieIndex}/index -U !{trimmedReadsDir}/*.fq -S alignedReads.sam
        '''
    }
}

/**
* Use samtools to remove PCR duplicates from the control, may use the old rmdup or new markdup -r, depending on parameters
*/
process removeControlPCRDuplicates {
    conda 'samtools'
    cpus params.threads

    input:
        path inputSam

    output:
        path 'deduplicatedControl.sam', emit: outputSam

    shell:
    threadCommand = params.threadSamtools ? "-@ ${task.cpus}" : ''

    if (params.useRmdup) {
        if (params.paired) {
            '''
            samtools rmdup !{inputSam} deduplicatedControl.sam
            '''
        } else {
            '''
            samtools rmdup -s !{inputSam} deduplicatedControl.sam
            '''
        }
    } else {
        '''
        samtools sort -n !{inputSam} -o sorted.sam !{threadCommand}
        samtools fixmate -m sorted.sam scored.sam !{threadCommand}
        samtools sort scored.sam -o scored_sorted.sam !{threadCommand}
        samtools markdup -r scored_sorted.sam deduplicatedControl.sam !{threadCommand}
        '''
    }
}

/**
* Call both broad and narrow peaks using the MACS2 peakcaller
*/
process macs2CallPeaksWithInput {
    conda 'macs2 biopython'
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path alignmentSam
        tuple file(inputFastqDir), val(inputFastq)
        path genomeFasta
        path controlAlignment

    output:
        path '*peakCallingResultsNarrow', emit: narrowPeaks
        path '*peakCallingResultsBroad', emit: broadPeaks

    shell:
    '''
    fastaLength=$(getFastaLength.py !{genomeFasta})
    macs2 callpeak -c !{controlAlignment} -g $fastaLength -q 0.00001 -t !{alignmentSam} --outdir !{inputFastq}_peakCallingResultsNarrow --name peaks
    macs2 callpeak -c !{controlAlignment} --broad -g $fastaLength -q 0.00001 -t !{alignmentSam} --outdir !{inputFastq}_peakCallingResultsBroad --name peaks
    '''
}

process macs2CallPeaks {
    conda 'macs2 biopython'
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path alignmentSam
        tuple file(inputFastqDir), val(inputFastq)
        path genomeFasta

    output:
        path '*peakCallingResultsNarrow', emit: narrowPeaks
        path '*peakCallingResultsBroad', emit: broadPeaks

    shell:
    '''
    fastaLength=$(getFastaLength.py !{genomeFasta})
    macs2 callpeak -g $fastaLength -q 0.00001 -t !{alignmentSam} --outdir !{inputFastq}_peakCallingResultsNarrow --name peaks
    macs2 callpeak --broad -g $fastaLength -q 0.00001 -t !{alignmentSam} --outdir !{inputFastq}_peakCallingResultsBroad --name peaks
    '''
}

process mergeReplicates {
    conda 'bedtools'
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path narrowPeakFiles
        path broadPeakFiles

    output:
        path 'mergedNarrowPeaks.bed'
        path 'mergedBroadPeaks.bed'

    shell:
    allNarrowPeakBEDs = findNarrowPeakBEDs(narrowPeakFiles)
    allBroadPeakBEDs = findBroadPeakBEDs(broadPeakFiles)

    '''
    cat !{allNarrowPeakBEDs.join(" ")} > unmergedNarrowPeaks.bed
    sort -k1,1 -k2,2n unmergedNarrowPeaks.bed > unmergedNarrowPeaks.sorted.bed
    bedtools merge -i unmergedNarrowPeaks.sorted.bed > mergedNarrowPeaks.bed
    cat !{allBroadPeakBEDs.join(" ")} > unmergedBroadPeaks.bed
    sort -k1,1 -k2,2n unmergedBroadPeaks.bed > unmergedBroadPeaks.sorted.bed
    bedtools merge -i unmergedBroadPeaks.sorted.bed > mergedBroadPeaks.bed
    '''
}

/**
* Detect the names of sets of paired end reads by looking for files with filenames *1.fastq and *2.fastq
*/
def intoFastqBasenames(inputFastqPath, isPaired) {
    def fastqs = [] as Set

    for (directory in inputFastqPath) {
        def baseDir = new File(directory)
        def files = baseDir.listFiles()

        for (fileObject in files) {
            def filename = fileObject.getName()
            def matches = null
            if (isPaired) {
                matches = (filename =~ /(.+)\d\.fastq/).findAll()
            } else {
                matches = (filename =~ /(.+)\.fastq/).findAll()
            }

            fastqs.add(new Tuple(file(directory), matches[0][1]))
        }
    }

    return Channel.fromList(fastqs)
}

def findNarrowPeakBEDs(narrowPeakDirs) {
    def narrowPeakBeds = []
    for (fileName in narrowPeakDirs) {
        narrowPeakBeds.add(fileName + '/peaks_peaks.narrowPeak')
    }

    return narrowPeakBeds
}

def findBroadPeakBEDs(broadPeakDirs) {
    def broadPeakBeds = []
    for (fileName in broadPeakDirs) {
        broadPeakBeds.add(fileName + '/peaks_peaks.broadPeak')
    }

    return broadPeakBeds
}

workflow {
    def genomeFile = Channel.value(file(params.genomeFasta))
    def fastqData = Channel.fromPath(params.chipSeqFastqList)
    def fastqBasenames = intoFastqBasenames(params.chipSeqFastqList, params.paired)

    buildBowtieIndex(genomeFile)

    def controlAlignment = null
    if (params.controlFastq != null) {
        def controlDir = file(params.controlFastq)

        trimControl(controlDir)
        alignControlToGenome(buildBowtieIndex.out.bowtieIndex, trimControl.out.trimmedReadsDir)
        removeControlPCRDuplicates(alignControlToGenome.out.alignedReads)
        
        controlAlignment = removeControlPCRDuplicates.out.outputSam
    }

    trimChipReplicates(fastqBasenames)
    alignToGenome(buildBowtieIndex.out.bowtieIndex, trimChipReplicates.out.trimmedReadsDir)
    removePCRDuplicates(alignToGenome.out.alignedReads)

    if (controlAlignment  == null) {
        macs2CallPeaks(removePCRDuplicates.out.outputSam, fastqBasenames, genomeFile)
        mergeReplicates(macs2CallPeaks.out.narrowPeaks.collect(), macs2CallPeaks.out.broadPeaks.collect())
    } else {
        macs2CallPeaksWithInput(removePCRDuplicates.out.outputSam, fastqBasenames, genomeFile, controlAlignment)
        mergeReplicates(macs2CallPeaksWithInput.out.narrowPeaks.collect(), macs2CallPeaksWithInput.out.broadPeaks.collect())
    }
}