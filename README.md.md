
# Table of Contents

1.  [&#x2013;help:          show an help message](#org25dfd6d)
2.  [&#x2013;dda<sub>folder</sub>:    folder with DDA files to be searched (default: Data/DDA)](#org23ffc58)
3.  [&#x2013;comet<sub>params</sub>:  comet parameter file (default: Params/comet.params)](#org463d475)
4.  [&#x2013;comet<sub>threads</sub>: number of cores to be used in comet search (default: 8)](#org0c24314)
5.  [&#x2013;protein<sub>db</sub>:    fasta formatted sequence database file (default: Results/Databases/proteome.fasta)](#org85aa573)
6.  [&#x2013;tpp:           options to pass to TPP xinteract (default: -OAPdplIw -PPM)](#org03f31ce)
7.  [&#x2013;decoy:         decoy prefix used in protein<sub>db</sub> (default: DECOY\_)](#org7bc4899)
8.  [&#x2013;no<sub>pool</sub>:       do not pool results at the TPP step (default: false)](#org6296479)
9.  [&#x2013;libra<sub>params</sub>:  libra parameter file (default: NO<sub>FILE</sub>)](#org06ef832)
10. [pep and prot xml files are stored under Results/Comet](#org03c6893)
11. [TSV exports of pep and prot XML (filtered at 1% FDR) are stored under Results/TppExport](#org4f00aa5)

\## Overview:

This is a [Nextflow](<https://www.nextflow.io/>) based pipeline for the assignment and validation of DDA mass spectrometry data.

[Comet](<http://comet-ms.sourceforge.net/>) is used for the assignment step.

The [Trans Proteomics Pipeline](<http://tools.proteomecenter.org/wiki/index.php?title=Software:TPP>) is used for peptide validation, protein inference, and quantification.

The default workflow will pool identifications from multiple database searches into aggregated pep and prot xml files:

.    mzXML files
.         +
.         |
. <del>--&#x2013;&#x2014;v--------</del>
. |DB engine search|
. <del>----------------</del>
. |
. |
. | <del>--------</del>
. <del>->search 1</del>&#x2013;+                           <del>+        +</del>
. | <del>--------</del>  |                           |          |
. |             | Pool  <del>----------------</del>  | <del>------</del> |
. <del>->  &#x2026;  +------&#x2013;&#x2014;> PeptideProphet +&#x2013;> |XPRESS| |
. |             |       +----------------</del>  | |Libra | |
. | <del>--------</del>  |                           | <del>------</del> |
. <del>->search 2</del>&#x2013;+                           <del>+   +    +</del>
.   <del>--------</del>                                   |
.                                                |
.                              <del>----</del>            |
.                              |Mayu<-+          |
.                              <del>----</del> | <del>---&#x2013;&#x2014;v-----</del>
.                                     <del>-+ProteinProphet|
.                       +-----------</del> | <del>--------------</del>
.                       |calctppstat<-+
.                       <del>-----------</del> |
.                                     |
.                           <del>-------</del> |
.                           |StPeter<-+
.                           <del>-------</del>

The **&#x2013;no<sub>pool</sub>** parameter can be used to bypass this behaviour and analyse all MS files independently:

.    mzXML files
.         +
.         |
. <del>--&#x2013;&#x2014;v--------</del>
. |DB engine search|
. <del>----------------</del>
. |                                                                <del>-------</del>
. |                                                              <del>->StPeter|
. |                                                              | +-------</del>
. |                                                              |
. |                                <del>+        +</del>                  | <del>----</del>
. |                                |          |                  <del>->Mayu|
. | +--------</del> <del>----------------</del>  | <del>------</del> | <del>--------------</del> | <del>----</del>
. <del>->search 1</del>-> PeptideProphet <del>&#x2013;> |XPRESS| +->ProteinProphet</del>-+ <del>-----------</del>
. | <del>--------</del> <del>----------------</del>  | |Libra | | <del>--------------</del> <del>->calctppstat|
. |                                | +------</del> |                    <del>-----------</del>
. <del>->  &#x2026;           &#x2026;           |   &#x2026;    |       &#x2026;          +----</del>
. |                                |          |                  <del>->Mayu|
. | +--------</del> <del>----------------</del>  | <del>------</del> | <del>--------------</del> | <del>----</del>
. <del>->search 2</del>-> PeptideProphet <del>&#x2013;> |XPRESS| +->ProteinProphet</del>-+ <del>-----------</del>
.   <del>--------</del> <del>----------------</del>  | |Libra | | <del>--------------</del> <del>->calctppstat|
.                                  | +------</del> |                  | <del>-----------</del>
.                                  <del>+        +</del>                  |
.                                                                | <del>-------</del>
.                                                                <del>->StPeter|
.                                                                  +-------</del>

\## Requirements:

\## Workflow parameters:


<a id="org25dfd6d"></a>

# &#x2013;help:          show an help message


<a id="org23ffc58"></a>

# &#x2013;dda<sub>folder</sub>:    folder with DDA files to be searched (default: Data/DDA)


<a id="org463d475"></a>

# &#x2013;comet<sub>params</sub>:  comet parameter file (default: Params/comet.params)


<a id="org0c24314"></a>

# &#x2013;comet<sub>threads</sub>: number of cores to be used in comet search (default: 8)


<a id="org85aa573"></a>

# &#x2013;protein<sub>db</sub>:    fasta formatted sequence database file (default: Results/Databases/proteome.fasta)


<a id="org03f31ce"></a>

# &#x2013;tpp:           options to pass to TPP xinteract (default: -OAPdplIw -PPM)


<a id="org7bc4899"></a>

# &#x2013;decoy:         decoy prefix used in protein<sub>db</sub> (default: DECOY\_)


<a id="org6296479"></a>

# &#x2013;no<sub>pool</sub>:       do not pool results at the TPP step (default: false)


<a id="org06ef832"></a>

# &#x2013;libra<sub>params</sub>:  libra parameter file (default: NO<sub>FILE</sub>)

\## Workflow results:


<a id="org03c6893"></a>

# pep and prot xml files are stored under Results/Comet


<a id="org4f00aa5"></a>

# TSV exports of pep and prot XML (filtered at 1% FDR) are stored under Results/TppExport

\## Examples:

