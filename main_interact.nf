nextflow.enable.dsl=2

include {run_xinteract;
	 tpp_summaries;
	 tpp_exports;
	 run_proteinprophet;
	 tpp_peter;
	 apply_progenesis_patch	 
} from './nfmodule_tpp/tpp_workflows.nf'

workflow {
    main:
    log.info("++++++++++========================================")
    log.info("Executing xinteract workflow")
    log.info("")
    log.info("Parameters:")
    log.info(" Search folder:\t $params.search_folder")
    log.info(" DDA folder:\t $params.dda_folder")
    log.info(" Protein DB:\t $params.protein_db")
    log.info(" Decoy tag:\t $params.decoy")
    log.info(" TPP params:\t $params.tpp")
    if(params.no_pool) {
	log.info("  All searches analyzed separately in TPP")
    }
    else {
	log.info("  All searches merged for TPP analysis")
    }
    log.info(" xinteract threads:\t $params.xinteract_threads")
    log.info("++++++++++========================================")

    pepxmls = channel.fromPath("$params.search_folder/*.pep.xml")
    if(!params.no_pool) {
	pepxmls = pepxmls.collect()
    }
    mzxmls = channel.fromPath("$params.dda_folder/*.{mzXML,mzML}")
    if(!params.no_pool) {
	mzxmls = mzxmls.collect()
    }
    
    run_xinteract(pepxmls,
		  params.protein_db,
		  mzxmls,
		  params.tpp,
		  params.decoy,
		  params.no_pool,
		  params.xinteract_threads)
        
    tpp_summaries(run_xinteract.out.pepxmls_rf,
		  params.comet_params,
		  params.protein_db,
		  params.decoy)

    if(params.no_pool) {
	// StPeter only makes sense if did not pool pep.xmls
	//
	// Channels in nextflow are FIFO and tpp_main.out.pepxmls and
	// tpp_main.out.protxmls are loaded with the corresponding pep
	// and prot xml files in the corresponding positions

	tpp_peter(run_xinteract.out.pepxmls_rf,
		  run_xinteract.out.protxmls_rf,
		  channel.fromPath("$params.dda_folder/*").collect())
    }

    if(params.export_to_tsv) {
	if(params.no_pool) {
	    // StPeter changes the prot.xml files in place we need to wait for it to finish
	    tpp_exports(run_xinteract.out.pepxmls_rf,
			run_xinteract.out.pepxmlmodels,
			tpp_peter.out.protxmls,
			run_xinteract.out.protxmlmodels)
	}
	else {
	    tpp_exports(run_xinteract.out.pepxmls_rf,
			run_xinteract.out.pepxmlmodels,
			run_xinteract.out.protxmls_rf,
			run_xinteract.out.protxmlmodels)
	}
    }

    if(params.patch_for_progenesis) {
	apply_progenesis_patch(run_xinteract.out.pepxmls_rf)
    }
}
