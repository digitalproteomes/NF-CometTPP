if(params.help) {
    log.info ""
    log.info "Comet - TPP workflows"
    log.info "--------------------"
    log.info ""
    log.info "Pooled workflow (default)"
    log.info ""
    log.info ".    mzXML files"
    log.info ".         +"
    log.info ".         |"
    log.info ". +-------v--------+"
    log.info ". |DB engine search|"
    log.info ". +----------------+"
    log.info ". |"
    log.info ". |"
    log.info ". | +--------+"
    log.info ". +->search 1+--+                           ++        ++"
    log.info ". | +--------+  |                           |          |"
    log.info ". |             | Pool  +----------------+  | +------+ |"
    log.info ". +->  ...  +-----------> PeptideProphet +--> |XPRESS| |"
    log.info ". |             |       +----------------+  | |Libra | |"
    log.info ". | +--------+  |                           | +------+ |"
    log.info ". +->search 2+--+                           ++   +    ++"
    log.info ".   +--------+                                   |"
    log.info ".                                                |"
    log.info ".                              +----+            |"
    log.info ".                              |Mayu<-+          |"
    log.info ".                              +----+ | +--------v-----+"
    log.info ".                                     +-+ProteinProphet|"
    log.info ".                       +-----------+ | +--------------+"
    log.info ".                       |calctppstat<-+"
    log.info ".                       +-----------+ |"
    log.info ".                                     |"
    log.info ".                           +-------+ |"
    log.info ".                           |StPeter<-+"
    log.info ".                           +-------+"
    log.info ""
    log.info ""	   
    log.info "Separate workflow (--no_pool):"
    log.info ""
    log.info ".    mzXML files"
    log.info ".         +"
    log.info ".         |"
    log.info ". +-------v--------+"
    log.info ". |DB engine search|"
    log.info ". +----------------+"
    log.info ". |                                                                +-------+"
    log.info ". |                                                              +->StPeter|"
    log.info ". |                                                              | +-------+"
    log.info ". |                                                              |"    
    log.info ". |                                ++        ++                  | +----+"
    log.info ". |                                |          |                  +->Mayu|"
    log.info ". | +--------+ +----------------+  | +------+ | +--------------+ | +----+"
    log.info ". +->search 1+-> PeptideProphet +--> |XPRESS| +->ProteinProphet+-+ +-----------+"
    log.info ". | +--------+ +----------------+  | |Libra | | +--------------+ +->calctppstat|"
    log.info ". |                                | +------+ |                    +-----------+"
    log.info ". +->  ...           ...           |   ...    |       ...          +----+"
    log.info ". |                                |          |                  +->Mayu|"
    log.info ". | +--------+ +----------------+  | +------+ | +--------------+ | +----+"
    log.info ". +->search 2+-> PeptideProphet +--> |XPRESS| +->ProteinProphet+-+ +-----------+"
    log.info ".   +--------+ +----------------+  | |Libra | | +--------------+ +->calctppstat|"
    log.info ".                                  | +------+ |                  | +-----------+"
    log.info ".                                  ++        ++                  |"
    log.info ".                                                                | +-------+"
    log.info ".                                                                +->StPeter|"
    log.info ".                                                                  +-------+"
    log.info ""
    log.info ""    
    log.info ""
    log.info "Options:"
    log.info "*  --help:          show this message and exit"
    log.info "*  --dda_folder:    folder with DDA files to be searched (default: $params.dda_folder)"
    log.info "*  --comet_params:  comet parameter file (default: $params.comet_params)"
    log.info "*  --comet_threads: number of cores to be used in comet search (default: $params.comet_threads)"
    log.info "*  --protein_db:    fasta formatted sequence database file (default: $params.protein_db)"
    log.info "*  --tpp:           options to pass to TPP xinteract (default: $params.tpp)"
    log.info "*  --decoy:         decoy prefix used in protein_db (default: $params.decoy)"
    log.info "*  --no_pool:       do not pool results at the TPP step (default: $params.no_pool)"
    log.info "*  --libra_params:  libra parameter file (default: $params.libra_params)"
    log.info ""
    log.info "Results will be stored in Results/Comet"
    log.info "TSV export of pep and prot XML (filtered at 1% FDR) will be under Results/TppExport"
    log.info ""
    exit 1
}


Channel.fromPath("${params.dda_folder}/*.mzXML")
    .map{ file ->
         def key = file.name.toString()
         return tuple(key, file)
    }
    .groupTuple()
    .set{ keyed_mzxmls }


process cometSearch {
    // Search all mzXML files in $params.dda_folder with Comet
    cpus params.comet_threads
    // Human, 3 variable mods, semi, 2 missed cleavages and some margin for safety
    memory 30.GB
    
    publishDir 'Results/Comet', mode: 'link'

    input:
    set key, file(mzXML) from keyed_mzxmls
    file comet_params from file(params.comet_params)
    file protein_db from file(params.protein_db)

    output:
    set key, file('*.pep.xml') into cometOut
    file mzXML into cometMzXMLOut

    """
    # Set proteins DB
    sed -i s,db_path,$protein_db, $comet_params
    sed -i 's,num_threads = 0,num_threads = ${params.comet_threads},' $comet_params

    comet -P$comet_params $mzXML
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g ${mzXML.simpleName}.pep.xml
    sed -ri 's|<search_database local_path="|<search_database local_path="${workflow.launchDir}/Results/Comet/|'g ${mzXML.simpleName}.pep.xml
    """
}


// Overview of the files created by xinteract:
// ------------------------------------------
//
// For a standard pooled run tppProtOutPtm and tppProtOutIpro are
// going to be empty
// NOTE: these are also going to be empty if the xinteract parameters
// don't request ProteinProphet in combination with iProphet or PTMProphet

// For a pooled PTMProphet run we are going to have
//
// * PeptideProphet
// comet_merged.pep.xml
// comet_merged.pep-MODELS.html
// comet_merged.pep.xml.index
// comet_merged.pep.xml.pIstats

// * iProphet
// comet_merged.ipro.pep.xml
// comet_merged.ipro.pep-MODELS.html
// comet_merged.ipro.pep.xml.index

// * PTMProphet
// comet_merged.ptm.ipro.pep.xml
// comet_merged.ptm.ipro.pep-MODELS.html
// comet_merged.ptm.ipro.pep.xml.index

// * ProteinProphet
// comet_merged.ptm.ipro.prot.xml
// comet_merged.ptm.ipro.prot-MODELS.html
// comet_merged.ptm.ipro.prot.xml

// For a pooled iProphet run we are going to have
//
// * PeptideProphet
// comet_merged.pep.xml
// comet_merged.pep-MODELS.html
// comet_merged.pep.xml.index
// comet_merged.pep.xml.pIstats

// * iProphet
// comet_merged.ipro.pep.xml
// comet_merged.ipro.pep-MODELS.html
// comet_merged.ipro.pep.xml.index

// * ProteinProphet
// comet_merged.ipro.prot.xml
// comet_merged.ipro.prot-MODELS.html
// comet_merged.ipro.prot.xml

if(!params.no_pool) {
    // Aggregate individual search results into a merged TPP analysis
    
    // We need to handle simple, iProphet and PTMProphet cases
    // separately to account for the different output files that will
    // be generated
    if ( params.tpp.indexOf("-M") != -1 ) {
	// PTMProphet
	process pooledTppPtm {
	    publishDir 'Results/Comet', mode: 'link'
	    
	    input:
	    file pepxmls from cometOut.map{ it[1] }.collect()
            file protein_db from file(params.protein_db)
	    file mzXML from cometMzXMLOut.collect()
	    file libra_params from file(params.libra_params)

	    output:
	    file 'comet_merged.ptm.ipro.pep.xml' into tppPepOutRaw
	    file 'comet_merged.ptm.ipro.pep-MODELS.html' // NOTE we
							 // are not
							 // using this
							 // one
							 // because
							 // the Error
							 // Table has
							 // a bug and
							 // the
							 // extracted
							 // probabilities
							 // are wrong
	    file 'comet_merged.ipro.pep-MODELS.html' into tppPepModelOut
	    file 'comet_merged.pep-MODELS.html'
	    file 'comet_merged.ptm.pep.xml.index' optional true
	    file 'comet_merged.ptm.ipro.prot-MODELS.html' into tppProtModelOut
	    file 'comet_merged.ptm.ipro.prot.xml' into tppProtOutRaw
	    file(protein_db) // Required for ProteinProphet visualization
	    
	    // xinteract for PTMProphet
	    """
            xinteract $params.tpp -d$params.decoy -Ncomet_merged.pep.xml $pepxmls
            """
	}	
    }
    else if ( params.tpp.indexOf("-i") != -1 ) {
	// iProphet
	process pooledTppIpro {
	    publishDir 'Results/Comet', mode: 'link'
	    
	    input:
	    file pepxmls from cometOut.map{ it[1] }.collect()
            file protein_db from file(params.protein_db)
	    file mzXML from cometMzXMLOut.collect()
	    file libra_params from file(params.libra_params)

	    output:
	    file 'comet_merged.ipro.pep.xml' into tppPepOutRaw
	    file 'comet_merged.ipro.pep-MODELS.html' into tppPepModelOut
	    file 'comet_merged.pep-MODELS.html'
	    file 'comet_merged.ipro.pep.xml.index' optional true
	    file 'comet_merged.ipro.prot-MODELS.html' optional true into tppProtModelOut
	    file 'comet_merged.ipro.prot.xml' optional true into tppProtOutRaw
	    file(protein_db) // Required for ProteinProphet visualization
	    
	    // xinteract for iProphet
	    """
            xinteract $params.tpp -d$params.decoy -Ncomet_merged.pep.xml $pepxmls
            """
	}	
    } 
    else {
	// Simple
	process pooledTpp {
	    publishDir 'Results/Comet', mode: 'link'
	    
	    input:
	    file pepxmls from cometOut.map{ it[1] }.collect()
            file protein_db from file(params.protein_db)
	    file mzXML from cometMzXMLOut.collect()
	    file libra_params from file(params.libra_params)

	    output:
	    // Normal run
	    file 'comet_merged.pep.xml' into tppPepOutRaw
	    file 'comet_merged.pep-MODELS.html' into tppPepModelOut
	    file 'comet_merged.pep.xml.index'
	    file 'comet_merged.pep.xml.pIstats'
	    file 'comet_merged.prot-MODELS.html' optional true into tppProtModelOut
	    file 'comet_merged.prot.xml' optional true into tppProtOutRaw
	    file(protein_db) // Required for ProteinProphet visualization
	    
	    // xinteract
	    """
            xinteract $params.tpp -d$params.decoy -Ncomet_merged.pep.xml $pepxmls
            """
	}
    }
}
else {
    // Perform a separate TPP analysis for each search result

    if ( params.tpp.indexOf("-M") != -1 ) {
	//PTMProphet
	process splitTppPtm {
	    publishDir 'Results/Comet', mode: 'link'
	    
	    input:
	    set key, file(pepxml) from cometOut
	    file protein_db from file(params.protein_db)
	    file mzXML from cometMzXMLOut.collect() // IMPROVE: We don't actually need them all.
	    file libra_params from file(params.libra_params)
	    
	    output:
	    set key, file('*_sep.ptm.ipro.pep.xml') into tppPepOutRaw
	    set key, file('*_sep.ptm.ipro.prot.xml') optional true into tppProtOutRaw
	    set key, file('*_sep.ptm.ipro.pep-MODELS.html') into tppPepModelOut
	    file '*_sep.ipro.pep-MODELS.html'
	    file '*_sep.pep-MODELS.html'
	    file '*_sep.ptm.ipro.pep.xml.index' optional true
	    set key, file('*_sep.ptm.ipro.prot-MODELS.html') optional true into tppProtModelOut
	    file(protein_db) // Required for ProteinProphet visualization

	    // xinteract and refactor links in prot.xml 
	    """
            xinteract $params.tpp -d$params.decoy -N${pepxml}_sep.pep.xml $pepxml
            """
	}
    }
    else if ( params.tpp.indexOf("-i") != -1 ) {
	//iProphet
	process splitTppIpro {
	    publishDir 'Results/Comet', mode: 'link'
	    
	    input:
	    set key, file(pepxml) from cometOut
	    file protein_db from file(params.protein_db)
	    file mzXML from cometMzXMLOut.collect() // IMPROVE: We don't actually need them all.
	    file libra_params from file(params.libra_params)
	    
	    output:
	    set key, file('*_sep.ipro.pep.xml') into tppPepOutRaw
	    set key, file('*_sep.ipro.prot.xml') optional true into tppProtOutRaw
	    set key, file('*_sep.ipro.pep-MODELS.html') into tppPepModelOut
	    file '*_sep.pep-MODELS.html'
	    file '*_sep.ipro.pep.xml.index' optional true
	    set key, file('*_sep.ipro.prot-MODELS.html') optional true into tppProtModelOut
	    file(protein_db) // Required for ProteinProphet visualization

	    // xinteract and refactor links in prot.xml 
	    """
            xinteract $params.tpp -d$params.decoy -N${pepxml}_sep.pep.xml $pepxml
            """
	}
    }
    else {
	// Simple
	process splitTpp {
	    publishDir 'Results/Comet', mode: 'link'
	    
	    input:
	    set key, file pepxml from cometOut
	    file protein_db from file(params.protein_db)
	    file mzXML from cometMzXMLOut.collect() // IMPROVE: We don't actually need them all.
	    file libra_params from file(params.libra_params)
	    
	    output:
	    set key, file('*_sep.pep.xml') into tppPepOutRaw
	    set key, file('*_sep.prot.xml') optional true into tppProtOutRaw
	    set key, file('*_sep.pep-MODELS.html') into tppPepModelOut
	    file '*_sep.pep.xml.index'
	    file '*_sep.pep.xml.pIstats'
	    set key, file('*_sep.prot-MODELS.html') optional true into tppProtModelOut
	    file(protein_db) // Required for ProteinProphet visualization

	    // xinteract and refactor links in prot.xml 
	    """
            xinteract $params.tpp -d$params.decoy -N${pepxml}_sep.pep.xml $pepxml
            """
	}
    }
}


// pepXml uses absolute links. Fix them after moving them to the
// publishDir
process refactorPepXml {
    publishDir 'Results/Comet', mode: 'link'

    input:
    set key, file(pepxml) from tppPepOutRaw

    output:
    set key, file(pepxml) into tppPepOut

    """
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g $pepxml
    """
}


// protXml uses absolute links. Fix them after moving them to the
// publishDir
process refactorProtXml {
    publishDir 'Results/Comet', mode: 'link'

    input:
    set key, file(protxml) from tppProtOutRaw

    output:
    set key, file(protxml) into tppProtOut

    """
    sed -ri 's|/work/.{2}/.{30}|/Results/Comet|'g $protxml 
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g $protxml
    """
}

// Multiply tppPepOut channel so we can feed it to multiple processes
tppPepOut.into{ tppPepOut1; tppPepOut2; tppPepOut3; tppPepOut4 }

// Multiply tppProtOut channel so we can feed it to multiple processes
tppProtOut.into{ tppProtOut1; tppProtOut2 }


// Run label-free quantification
process stPeter {
    input:
    tppPepOut3.join(tppProtOut1).set{ joined_stPeter }
    set key, file(pepxml), file(protxml) from joined_stPeter
    file mzXML from file("${params.dda_folder}/*.mzXML")

    output:
    set key, file(protxml) into stPeterOut
    
    script:
    """
    StPeter $protxml
    """
}


process mayu {
    // For each TPP analysis run Mayu
    publishDir 'Results/Comet', mode: 'link'

    errorStrategy 'ignore'
    
    input:
    file pepxml from tppPepOut1.map{ it[1] }
    file comet_params from (params.comet_params)
    file db from file(params.protein_db)
    
    output:
    file("mayu_*")
    
    """
    Mayu.pl -A $pepxml \
    -C $db \
    -E $params.decoy \
    -M mayu_$pepxml \
    -P pep FDR=0.01:1
    """
}


process tppStat {
    // For each TPP analysis run calctppstat
    publishDir 'Results/Comet', mode: 'link'
    
    input:
    tppPepOut2.join(tppProtOut2).set{ joined_tppStat }
    set key, file(pepxml), file(protxml) from joined_tppStat

    output:
    file '*.summary.txt' into tppStatOut
    
    """
    /usr/local/tpp/cgi-bin/calctppstat.pl -i $pepxml -d $params.decoy --full > ${pepxml}.summary.txt
    """
}


process pepxml2tsv {
    tag "$pepXml"
    
    publishDir 'Results/TppExport', mode: 'link'

    input:
    tppPepOut4.join(tppPepModelOut).set{ joined_pepxml2tsv }
    set key, file(pepXml), file(pepXmlModels) from joined_pepxml2tsv
    file pepxsl from file("$baseDir/Xslt/pepxml2tsv.xsl")
    
    output:
    set key, file ('*.tsv') into pepTsvOut

    """
    PROB=\$(get_prophet_prob.py -i $pepXmlModels)
    xsltproc --param p_threshold \$PROB $pepxsl $pepXml > ${pepXml}.tsv
    """
}


process protxml2tsv {
    tag "$protXml"
    
    publishDir 'Results/TppExport', mode: 'link'

    input:
    //tppProtOut3.join(tppProtModelOut).set{ joined_protxml2tsv }
    stPeterOut.join(tppProtModelOut).set{ joined_protxml2tsv }
    set key, file(protXml), file(protXmlModels) from joined_protxml2tsv
    file protxsl from file("$baseDir/Xslt/protxml2tsv.xsl")
    
    output:
    set key, file('*.tsv') into protTsvOut

    """
    PROB=\$(get_prophet_prob.py -i $protXmlModels)
    xsltproc --param p_threshold \$PROB $protxsl $protXml > ${protXml}.tsv
    """
}


// Filters the peptide list created by pepxml2tsv by removing peptides
// not assigned to a protein included in the protxml2tsv generated
// protein list.
process filterPepTsv {
    tag "$pepTsv - $protTsv"
    
    publishDir 'Results/TppExport', mode: 'link'

    input:
    pepTsvOut.join(protTsvOut).set{ joined_filterPepTsv }
    set key, file(pepTsv), file(protTsv) from joined_filterPepTsv

    output:
    file '*.tsv'
    
    script:
    """
    filter_pep_pro.py -p $pepTsv -P $protTsv -o ${pepTsv.baseName}_filtered.tsv
    """
}


workflow.onComplete {
    // Make the Comet results folder writable for the www-data group
    // for TPP visualization.
    "chmod g+w Results/Comet".execute()
}
