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
    log.info "  --mzxml_folder: files to be searched (default: $params.mzxml_folder)"
    log.info "  --comet_params: comet parameter file (default: $params.comet_params)"
    log.info "  --protein_db:   comet parameter file (default: $params.protein_db)"
    log.info "  --tpp:          TPP options (default: $params.tpp)"
    log.info "  --decoy:        decoy prefix (default: $params.decoy)"
    log.info "  --no_pool:      do not pool results at the TPP step (default: $params.no_pool)"
    log.info ""
    log.info "Results will be stored in Results/Comet"
    log.info ""
    exit 1
}



process cometSearch {
    // Search all mzXML files in $params.mzxml_folder with Comet
    publishDir 'Results/Comet'
    
    input:
    file mzXML from file(params.mzxml_folder)
    file comet_params from file(params.comet_params)
    file protein_db from file(params.protein_db)

    output:
    file '*.pep.xml' into cometOut
    file mzXML

    """
    comet $mzXML
    """
}


if(!params.no_pool) {
    // Aggregate individual search results into a merged TPP analysis
    process pooledTpp {
	publishDir 'Results/Comet'
	
	input:
	file pepxmls from cometOut.collect()
        file protein_db from file(params.protein_db)

	output:
	file 'comet_merged.pep.xml' into tppPepOut
	file 'comet_merged.pep-MODELS.html'
	file 'comet_merged.pep.xml.index'
	file 'comet_merged.pep.xml.pIstats'
	file 'comet_merged.prot-MODELS.html'
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
	publishDir 'Results/Comet'
	
	input:
	file pepxml from cometOut
	file protein_db from file(params.protein_db)

	output:
	file '*_sep.pep.xml' into tppPepOut
	file '*_sep.prot.xml' into tppProtOut
	file '*_sep.pep-MODELS.html'
	file '*_sep.pep.xml.index'
	file '*_sep.pep.xml.pIstats'
	file '*_sep.prot-MODELS.html'
	file(protein_db) // Required for ProteinProphet visualization

	// xinteract and refactor links in prot.xml 
	"""
        xinteract $params.tpp -d$params.decoy -N${pepxml}_sep.pep.xml $pepxml
	sed -ri 's|/work/.{2}/.{30}|/Results/Comet|'g ${pepxml}_sep.prot.xml
        """
    }
}

// Duplicate tppPepOut channel so we can feed it to two processes
tppPepOut.into{ tppPepOut1; tppPepOut2}


process mayu {
    // For each TPP analysis run Mayu
    publishDir 'Results/Comet'

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
    publishDir 'Results/Comet'
    
    input:
    file pepxml from tppPepOut2
    file protxml from tppProtOut

    output:
    file '*.summary.txt' into tppStatOut
    
    """
    /usr/local/tpp/cgi-bin/calctppstat.pl -i $pepxml -d $params.decoy --full > ${pepxml}.summary.txt
    """
}

workflow.onComplete {
    // Make the Comet results folder writable for the www-data group
    // for TPP visualization.
    "chmod g+w Results/Comet".execute()
}
