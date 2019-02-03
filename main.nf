if(params.help) {
    log.info ""
    log.info "Comet - TPP workflow"
    log.info "--------------------"
    log.info ""
    log.info "Options:"
    log.info "  --help:         show this message and exit"
    log.info "  --mzxml_folder: files to be searched (default: $params.mzxml_folder)"
    log.info "  --comet_params: comet parameter file (default: $params.comet_params)"
    log.info "  --protein_db:   comet parameter file (default: $params.comet_params)"
    log.info "  --tpp:          TPP options (default: $params.tpp_opt)"
    log.info "  --decoy:        decoy prefix (default: $params.decoy)"
    log.info "  --no_pool:      do not pool results at the TPP step (default: $params.no_pool)"
    log.info ""
    log.info "Results will be stored in Results/Comet"
    log.info ""
    exit 1
}



process cometSearch {
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
    // Aggregate individual search results into a merged pep.xml
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

	"""
        xinteract $params.tpp -d$params.decoy -Ncomet_merged.pep.xml $pepxmls
        """
    }
}
else {
    // Separate TPP analysis for each search result
    process splitTpp {
	publishDir 'Results/Comet'
	
	input:
	file pepxml from cometOut
	file protein_df from file(params.protein_db)

	output:
	file '*.pep.xml' into tppPepOut
	file '*.prot.xml' into tppProtOut

	"""
        xinteract $params.tpp -d$params.decoy $pepxml
        """
    }
}

// Duplicate tppPepOut channel so we can feed it to two processes
tppPepOut.into{ tppPepOut1; tppPepOut2}


process mayu {
    publishDir 'Results/Comet'

    input:
    file pepxml from tppPepOut1
    file comet_params from (params.comet_params)
    file params.protein_db
    
    output:
    file("mayu_*")
    
    """
    Mayu.pl -A $pepxml \
    -C $params.protein_db \
    -E $params.decoy \
    -M mayu_$pepxml \
    -P pepFDR=0.01:1
    """
}


process tppStat {
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
