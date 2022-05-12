# Nextflow workflow for searching DDA MS data through Comet and the Trans Proteomics Pipeline


## Overview:

This is a [Nextflow](https://www.nextflow.io/) based pipeline for the assignment and validation of DDA mass spectrometry data.

[Comet](http://comet-ms.sourceforge.net/) is used for the assignment step.

The [Trans Proteomics Pipeline](http://tools.proteomecenter.org/wiki/index.php?title=Software:TPP) is used for peptide validation, protein inference, and quantification.

The default workflow will pool identifications from multiple database searches into aggregated pep and prot xml files:

       mzXML files
            +
            |
    +-------v--------+
    |DB engine search|
    +----------------+
    |
    |
    | +--------+
    +->search 1+--+                           ++        ++
    | +--------+  |                           |          |
    |             | Pool  +----------------+  | +------+ |
    +->  ...  +-----------> PeptideProphet +--> |XPRESS| |
    |             |       +----------------+  | |Libra | |
    | +--------+  |                           | +------+ |
    +->search 2+--+                           ++   +    ++
      +--------+                                   |
                                                   |
                                 +----+            |
                                 |Mayu<-+          |
                                 +----+ | +--------v-----+
                                        +-+ProteinProphet|
                          +-----------+ | +--------------+
                          |calctppstat<-+
                          +-----------+ |
                                        |
                              +-------+ |
                              |StPeter<-+
                              +-------+


The `--no_pool` parameter can be used to bypass this behaviour and analyze all MS files independently:

       mzXML files
            +
            |
    +-------v--------+
    |DB engine search|
    +----------------+
    |                                                                +-------+
    |                                                              +->StPeter|
    |                                                              | +-------+
    |                                                              |
    |                                ++        ++                  | +----+
    |                                |          |                  +->Mayu|
    | +--------+ +----------------+  | +------+ | +--------------+ | +----+
    +->search 1+-> PeptideProphet +--> |XPRESS| +->ProteinProphet+-+ +-----------+
    | +--------+ +----------------+  | |Libra | | +--------------+ +->calctppstat|
    |                                | +------+ |                    +-----------+
    +->  ...           ...           |   ...    |       ...          +----+
    |                                |          |                  +->Mayu|
    | +--------+ +----------------+  | +------+ | +--------------+ | +----+
    +->search 2+-> PeptideProphet +--> |XPRESS| +->ProteinProphet+-+ +-----------+
      +--------+ +----------------+  | |Libra | | +--------------+ +->calctppstat|
                                     | +------+ |                  | +-----------+
                                     ++        ++                  |
                                                                   | +-------+
                                                                   +->StPeter|
                                                                     +-------+


## Workflow parameters:

*  `--dda_folder`:    folder with DDA files to be searched (default: Data/DDA)
*  `--comet_params`:  comet parameter file (default: Params/comet.params)
*  `--comet_threads`: number of cores to be used in comet search (default: 20)
*  `--protein_db`:    fasta formatted sequence database file (default: Results/Databases/proteome.fasta)
*  `--tpp`:           options to pass to TPP xinteract (default: -OAPdplIw -PPM)
*  `--decoy`:         decoy prefix used in protein_db (default: DECOY_)
*  `--no_pool`:       do not pool results at the TPP step (default: false)
*  `--xinteract_threads`:       number of cores to be used at xinteract step (default: 50)


## Workflow results:

* pep and prot xml files are stored under Results/Comet
* TSV exports of pep (pep.xml.tsv) and prot (prot.xml.tsv) XML (filtered at 1% FDR) are stored under Results/TppExport. A peptide level TSV export filtered at 1% FDR for both protein and peptides is also created in the same folder (pep.xml_filtered.tsv)


## Requirements:

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/)


## Examples:

Run using default settings:
```sh
nextflow run digitalproteomes/NF-CometTPP
```

Specify MS files location and analyze all them separately:
```sh
nextflow run digitalproteomes/NF-CometTPP --dda_folder /mnt/ms_data/ --no_pool
```
