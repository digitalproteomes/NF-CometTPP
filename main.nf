if(params.help) {
    log.info ""
    log.info "Comet - TPP workflows"
    log.info "--------------------"
    log.info ""
    log.info "Pooled workflow"
    log.info ""
    log.info "   mzXML files"
    log.info "        +"
    log.info "        |"
    log.info "+-------v--------+"
    log.info "|DB engine search|"
    log.info "+----------------+"
    log.info "|"
    log.info "|"
    log.info "| +--------+"
    log.info "+->search 1+--+                           ++        ++"
    log.info "| +--------+  |                           |          |"
    log.info "|             | Pool  +----------------+  | +------+ |"
    log.info "+->  ...  +-----------> PeptideProphet +--> |XPRESS| |"
    log.info "|             |       +----------------+  | +------+ |"
    log.info "| +--------+  |                           |          |"
    log.info "+->search 2+--+                           ++   +    ++"
    log.info "  +--------+                                   |"
    log.info "                                               |"
    log.info "                             +----+            |"
    log.info "                             |Mayu<-+          |"
    log.info "                             +----+ | +--------v-----+"
    log.info "                                    +-+ProteinProphet|"
    log.info "                      +-----------+ | +--------------+"
    log.info "                      |calctppstat<-+"
    log.info "                      +-----------+"
    log.info ""
    log.info ""	   
    log.info "Separate workflow (--no_pool):"
    log.info ""
    log.info "   mzXML files"
    log.info "        +"
    log.info "        |"
    log.info "+-------v--------+"
    log.info "|DB engine search|"
    log.info "+----------------+"
    log.info "|"
    log.info "|                                ++        ++                    +----+"
    log.info "|                                |          |                  +->Mayu|"
    log.info "| +--------+ +----------------+  | +------+ | +--------------+ | +----+"
    log.info "+->search 1+-> PeptideProphet +--> |XPRESS| +->ProteinProphet+-+ +-----------+"
    log.info "| +--------+ +----------------+  | +------+ | +--------------+ +->calctppstat|"
    log.info "|                                |          |                    +-----------+"
    log.info "+->  ...           ...           |   ...    |       ...          +----+"
    log.info "|                                |          |                  +->Mayu|"
    log.info "| +--------+ +----------------+  | +------+ | +--------------+ | +----+"
    log.info "+->search 2+-> PeptideProphet +--> |XPRESS| +->ProteinProphet+-+ +-----------+"
    log.info "  +--------+ +----------------+  | +------+ | +--------------+ +->calctppstat|"
    log.info "                                 |          |                    +-----------+"
    log.info "                                 ++        ++"
    log.info ""
    log.info ""
    log.info "Options:"
    log.info "  --help:         show this message and exit"
    log.info "  --dda_folder: files to be searched (default: $params.dda_folder)"
    log.info "  --comet_params: comet parameter file (default: $params.comet_params)"
    log.info "  --protein_db:   comet parameter file (default: $params.protein_db)"
    log.info "  --tpp:          TPP options (default: $params.tpp)"
    log.info "  --decoy:        decoy prefix (default: $params.decoy)"
    log.info "  --no_pool:      do not pool results at the TPP step (default: $params.no_pool)"
    log.info ""
    log.info "Results will be stored in Results/Comet"
    log.info "TSV export of pep and prot XML (filtered at 1% FDR) will be under Results/TppExport"
    log.info ""
    exit 1
}



process cometSearch {
    // Search all mzXML files in $params.dda_folder with Comet
    cpus params.comet_threads
    // Human, 3 variable mods, semi, 2 missed cleavages and some margin for safety
    memory 30.GB
    
    publishDir 'Results/Comet', mode: 'link'

    input:
    file mzXML from file("${params.dda_folder}/*.mzXML")
    file comet_params from file(params.comet_params)
    file protein_db from file(params.protein_db)

    output:
    file '*.pep.xml' into cometOut
    file mzXML into cometMzXMLOut

    """
    # Set proteins DB
    sed -i s,db_path,$protein_db, $comet_params
    sed -i 's,num_threads = 0,num_threads = ${params.comet_threads},' $comet_params

    comet $mzXML
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g ${mzXML.simpleName}.pep.xml
    sed -ri 's|<search_database local_path="|<search_database local_path="${workflow.launchDir}/Results/Comet/|'g ${mzXML.simpleName}.pep.xml
    """
}


if(!params.no_pool) {
    // Aggregate individual search results into a merged TPP analysis
    process pooledTpp {
	publishDir 'Results/Comet', mode: 'link'
	
	input:
	file pepxmls from cometOut.collect()
        file protein_db from file(params.protein_db)
	file mzXML from cometMzXMLOut.collect()
	file libra_params from file(params.libra_params)

	output:
	file 'comet_merged.pep.xml' into tppPepOut
	file 'comet_merged.pep-MODELS.html' into tppPepModelOut
	file 'comet_merged.pep.xml.index'
	file 'comet_merged.pep.xml.pIstats'
	file 'comet_merged.prot-MODELS.html' into tppProtModelOut
	file 'comet_merged.prot.xml' into tppProtOut
	file(protein_db) // Required for ProteinProphet visualization

	// xinteract and refactor links in prot.xml 
	"""
        xinteract $params.tpp -d$params.decoy -Ncomet_merged.pep.xml $pepxmls
	sed -ri 's|/work/.{2}/.{30}|/Results/Comet|'g comet_merged.prot.xml
        """
    }
}
else {
    // Perform a separate TPP analysis for each search result
    process splitTpp {
	publishDir 'Results/Comet', mode: 'link'
	
	input:
	file pepxml from cometOut
	file protein_db from file(params.protein_db)
	file mzXML from cometMzXMLOut.collect() // IMPROVE: We don't actually need them all.
	file libra_params from file(params.libra_params)
	
	output:
	file '*_sep.pep.xml' into tppPepOut
	file '*_sep.prot.xml' into tppProtOut
	file '*_sep.pep-MODELS.html' into tppPepModelOut
	file '*_sep.pep.xml.index'
	file '*_sep.pep.xml.pIstats'
	file '*_sep.prot-MODELS.html' into tppProtModelOut
	file(protein_db) // Required for ProteinProphet visualization

	// xinteract and refactor links in prot.xml 
	"""
        xinteract $params.tpp -d$params.decoy -N${pepxml}_sep.pep.xml $pepxml
	sed -ri 's|/work/.{2}/.{30}|/Results/Comet|'g ${pepxml}_sep.prot.xml
	sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g ${pepxml}_sep.pep.xml
	sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g ${pepxml}_sep.prot.xml
        """
    }
}

// Duplicate tppPepOut channel so we can feed it to two processes
tppPepOut.into{ tppPepOut1; tppPepOut2; tppPepOut3; tppPepOut4}

tppProtOut.into{ tppProtOut1; tppProtOut2; tppProtOut3}


process stPeter {
    input:
    file protxml from tppProtOut1
    file pepxml from tppPepOut3
    file mzXML from file("${params.dda_folder}/*.mzXML")

    output:

    
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
    file pepxml from tppPepOut1
    file comet_params from (params.comet_params)
    file db from file(params.protein_db)
    
    output:
    file("mayu_*")
    
    """
    Mayu.pl -A $pepxml \
    -C $db \
    -E $params.decoy \
    -M mayu_$pepxml \
    -P pepFDR=0.01:1
    """
}


process tppStat {
    // For each TPP analysis run calctppstat
    publishDir 'Results/Comet', mode: 'link'
    
    input:
    file pepxml from tppPepOut2
    file protxml from tppProtOut2

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
    file pepXmlModels from tppPepModelOut
    file pepXml from tppPepOut4
    file pepxsl from file("$baseDir/Xslt/pepxml2tsv.xsl")
    
    output:
    file '*.tsv'

    """
    PROB=\$(get_prophet_prob.py -i $pepXmlModels)
    xsltproc --param p_threshold \$PROB $pepxsl $pepXml > ${pepXml}.tsv
    """
}


process protxml2tsv {
    tag "$protXml"
    
    publishDir 'Results/TppExport', mode: 'link'

    input:
    file protXmlModels from tppProtModelOut
    file protXml from tppProtOut3
    file protxsl from file("$baseDir/Xslt/protxml2tsv.xsl")
    
    output:
    file '*.tsv'

    """
    PROB=\$(get_prophet_prob.py -i $protXmlModels)
    xsltproc --param p_threshold \$PROB $protxsl $protXml > ${protXml}.tsv
    """
}


workflow.onComplete {
    // Make the Comet results folder writable for the www-data group
    // for TPP visualization.
    "chmod g+w Results/Comet".execute()
}
