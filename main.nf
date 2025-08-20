/*
vim: syntax=groovy
-*- mode: groovy;-*-

==============================
IPAW: a pipeline to search smallORF database
==============================
@Authors
Yafeng Zhu @yafeng

https://github.com/yafeng/smORF-workflow

- make sure conda install all excutive programs or scripts first
- create a bin folder to store all scripts
*/

nf_required_version = '0.26.0'
if( ! nextflow.version.matches(">= ${nf_required_version}") ){
  println("Nextflow version too old, ${nf_required_version} required")
  exit(1)
}


/* SET DEFAULT PARAMS */
params.isobaric = false
params.outdir = 'result'

mods = file(params.mods)
msgf_db = file(params.tdb)
msgf_jar = file(params.msgfjar)
converter = file(params.mzid2tsvConverter)

params.PrecursorMassTolerance = '10ppm'
params.FragmentMethodID = 3 // 3 for HCD spectra, 1 for CID, 2 for ETD
params.inst = 3 // 0: Low-res LCQ/LTQ (Default), 1: Orbitrap/FTICR/Lumos, 2: TOF, 3: Q-Exactive
plextype = params.isobaric ? params.isobaric.replaceFirst(/[0-9]+plex/, "") : false
msgfprotocol = params.isobaric ? [tmt:4, itraq:2][plextype] : 0

activationtype = [HCD:'High-energy collision-induced dissociation',CID:'Collision-induced dissociation',ETD:'Electron transfer dissociation',no:''][params.activationFilter]
massshifts = [tmt:0.0013, itraq:0.00125, false:0]
massshift = massshifts[plextype]

params.qval = 0.05 
params.MS2error = 0.02 // unit: Da

params.novheaders = '^SEP;^decoy_SEP' 
params.noclassfdr = false
novheaders = params.novheaders == true ? false : params.novheaders


///////////////////

println "PrecursorMassTolerance used in MSGF+:${params.PrecursorMassTolerance}"
println "q-value cutoff:${params.qval}"
println "MS2 error tolerance used in SpectrumAI:${params.MS2error}"
println "IsobaricAnalyzer filtered by activation:${params.activationFilter}"
println "InstrumentID used in MSGF+: ${params.inst}"
println "FragmentMethodID used in MSGF+: ${params.FragmentMethodID}"
println "ProtocolID used in MSGF+: $msgfprotocol"


/* PIPELINE START */

// Either feed an mzmldef file (tab separated lines with filepath\tsetname)
Channel
  .from(file("${params.mzmldef}").readLines())
  .map { it -> it.tokenize('\t') }
  .set { mzml_in }


mzml_in
  .tap { sets }
  .map { it -> [ it[1], (it[0]=~/.*\/(.+)\..*$/)[0][1], file(it[0])]} // create, set, samplename and file
  .tap{ mzmlfiles; mzml_msgf ; mzml_isobaric }
  .count()
  .set{ amount_mzml }

sets
  .map{ it -> it[1] }
  .unique()
  .tap { sets_for_denoms }
  .collect()
  .subscribe { println "Detected setnames: ${it.join(', ')}" }


// Get denominators if isobaric experiment
// passed in form --denoms 'set1:126:128N set2:131 set4:129N:130C:131'

if (params.isobaric && params.denoms) {
  setdenoms = [:]
  params.denoms.tokenize(' ').each{ it -> x=it.tokenize(':'); setdenoms.put(x[0], x[1..-1])}
  set_denoms = Channel.value(setdenoms)

}else if (params.isobaric) {
  setdenoms = [:]
  sets_for_denoms.reduce(setdenoms){ a, b -> a.put(b, ['_126']); return a }.set{ set_denoms }
}


process IsobaricQuant {

  when: params.isobaric

  input:
  set val(setname), val(sample), file(infile) from mzml_isobaric

  output:
  set val(sample), file("${infile}.consensusXML") into isobaricxml

  """
  source activate openms-2.4.0
  IsobaricAnalyzer  -type $params.isobaric -in $infile -out "${infile}.consensusXML" -extraction:select_activation "$activationtype" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true 
  """
}

isobaricamount = params.isobaric ? amount_mzml.value : 1

isobaricxml
  .ifEmpty(['NA', 'NA'])
  .buffer(size: isobaricamount)
  .flatMap { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> it[1] }
  .collect()
  .set { sorted_isoxml }

mzmlfiles
  .tap { groupset_mzmls }
  .buffer(size: amount_mzml.value)
  .map { it.sort( {a, b -> a[1] <=> b[1]}) }  //to sort by the first element of a tuple in descending order
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }] }
  .set{ mzmlfiles_all }

process createSpectraLookup {

  input:
  file(isobxmls) from sorted_isoxml 
  set val(setnames), file(mzmlfiles) from mzmlfiles_all
  
  output:
  file('mslookup_db.sqlite') into spec_lookup

  script:
  if(params.isobaric)
  """
  source activate nf-core-ddamsproteomics-1.0.0
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  msslookup isoquant --dbfile mslookup_db.sqlite -i ${isobxmls.join(' ')} --spectra ${mzmlfiles.join(' ')}
  """
  else
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  """
}

process msgfPlus {

  input:
  set val(setname), val(sample), file(x) from mzml_msgf

  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids
  
  """
  java -Xmx8000M -jar $msgf_jar -tasks 1 -thread 1 -d $msgf_db -s $x -o "${sample}.mzid" -mod $mods -tda 1 -t ${params.PrecursorMassTolerance} -ti -1,2 -m ${params.FragmentMethodID} -inst ${params.inst} -e 1 -protocol ${msgfprotocol} -ntt 2 -minLength 8 -maxLength 40 -minCharge 2 -maxCharge 4 -maxMissedCleavages 2 -n 1 -addFeatures 1
  """
}

process ConvertmzidTotsv {

  input:
  set val(setname), val(sample), file(x) from mzids
  
  output:
  set val(setname), file('out.mzid.tsv') into mzidtsvs 

  """
  source activate mono
  mono $converter -mzid:$x -tsv:out.mzid.tsv -showDecoy
  """
}

mzidtsvs
  .groupTuple()
  .set {mzidtsvs_byset}

process MergeTSVbyset {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(setname), file('tsv?') from mzidtsvs_byset

  output:
  file("${setname}_allpsms.txt") into merged_tsvs
  """
  head -1 tsv1 > psmheader  
  tail -q -n +2 tsv* > psmmerge
  cat psmheader psmmerge > ${setname}_allpsms.txt
  """
}

