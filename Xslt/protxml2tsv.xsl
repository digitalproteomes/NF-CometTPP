<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
		xmlns:xs="http://www.w3.org/2001/XMLSchema" exclude-result-prefixes="xs"
		xmlns:xd="http://www.oxygenxml.com/ns/doc/xsl" version="1.0"
		xmlns:protx="http://regis-web.systemsbiology.net/protXML">
  <!-- PRODUCES A TAB SEPARATED FILE FROM A PROT.XML FILE
       
       Takes one parameter to the define the protein prophet probability cutoff.

       Example: xsltproc -\-param p_threshold 0.5 <XSL_FILE> <XML_FILE>


       06-09-19: - Initial commit
       Patrick Pedrioli-->

  <xsl:output method="text" indent="no"/>
  <xsl:strip-space elements="*"/>

  <xsl:param name="p_threshold" select="1"/>
  
  <xsl:variable name="iprophet" 
		select="contains(protx:protein_summary/protx:protein_summary_header/protx:program_details/protx:proteinprophet_details/@run_options , 'IPROPHET')"/>
  <xsl:variable name="asap"
		select="boolean(protx:protein_summary/protx:analysis_summary[@analysis='asapratio'])"/>
  <xsl:variable name="asap_p"
		select="boolean(protx:protein_summary/protx:analysis_summary[@analysis='asapratio_pvalue'])"/>
  <xsl:variable name="xpress"
		select="boolean(protx:protein_summary/protx:analysis_summary[@analysis='xpress'])"/>
  <xsl:variable name="libra" select="boolean(protx:analysis_result[@analysis='libra'])"/>
  <xsl:variable name="stpeter"
		select="boolean(protx:protein_summary/protx:analysis_summary[@analysis='stpeter'])"/>

  
  <!-- Headers -->
  <xsl:template match="/">protein_group&#9;group_sibling_id&#9;protein&#9;subsuming_protein_entry&#9;protein_probability&#9;percent_coverage&#9;num_unique_peps&#9;percent_share_of_spectrum_id<xsl:if test="$asap">&#9;ratio_mean&#9;ratio_stdev&#9;ratio_num_peptides</xsl:if><xsl:if test="$asap_p">&#9;adjusted_ratio_mean&#9;adjusted_ratio_stdev&#9;pvalue</xsl:if><xsl:if test="$xpress">&#9;xpress_ratio_mean&#9;xpress_stdev&#9;xpress_num_peptides</xsl:if><xsl:if test="$libra">libra_mz&#9;libra_ratio&#9;libra_error</xsl:if><xsl:if test="$stpeter">&#9;stpepter_sin&#9;stpepter_counts</xsl:if>&#9;peptide_sequence&#9;non-degenerate<xsl:text>
</xsl:text>
<xsl:apply-templates/>
  </xsl:template>

  <xsl:template match="protx:protein_group">
    <xsl:for-each select="protx:protein">
      <xsl:if test="@probability &gt;= $p_threshold">
	<xsl:for-each select="protx:peptide">
	  <!-- ProteinProphet -->
          <xsl:value-of select="../../@group_number"/>
          <xsl:text>&#9;</xsl:text>
          <xsl:value-of select="../@group_sibling_id"/>
          <xsl:text>&#9;</xsl:text>
          <xsl:value-of select="../@protein_name"/>
          <xsl:if test="../@n_indistinguishable_proteins &gt; '1'">
            <xsl:for-each select="../protx:indistinguishable_protein">,<xsl:value-of
            select="../@protein_name"/></xsl:for-each>
          </xsl:if>
          <xsl:text>&#9;</xsl:text>
          <xsl:if test="../@subsuming_protein_entry">
            <xsl:value-of select="../@subsuming_protein_entry"/>
          </xsl:if>
          <xsl:text>&#9;</xsl:text>
          <xsl:value-of select="../@probability"/>
          <xsl:text>&#9;</xsl:text>
          <xsl:value-of select="../@percent_coverage"/>
          <xsl:text>&#9;</xsl:text>
          <xsl:value-of select="../@total_number_peptides"/>
          <xsl:text>&#9;</xsl:text>
          <xsl:value-of select="../@pct_spectrum_ids"/>

          <!-- ASAPRatio -->
          <xsl:if test="$asap">
            <xsl:if test="../protx:analysis_result[@analysis='asapratio']">
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'asapratio']/protx:ASAPRatio/@ratio_mean"/>
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'asapratio']/protx:ASAPRatio/@ratio_standard_dev"/>
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'asapratio']/protx:ASAPRatio/@ratio_number_peptides"
                  />
            </xsl:if>
            <xsl:if test="not(../protx:analysis_result[@analysis='asapratio'])">
              <xsl:text>&#9;</xsl:text>
              <xsl:text>&#9;</xsl:text>
              <xsl:text>&#9;</xsl:text>
            </xsl:if>
          </xsl:if>

          <xsl:if test="$asap_p">
            <xsl:if test="../protx:analysis_result[@analysis='asapratio_pvalue']">
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'asapratio_pvalue']/protx:ASAPRatio_pvalue/@adj_ratio_mean"/>
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'asapratio_pvalue']/protx:ASAPRatio_pvalue/@adj_ratio_standard_dev"/>
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'asapratio_pvalue']/protx:ASAPRatio_pvalue/@pvalue"
                  />
            </xsl:if>
            <xsl:if test="not(../protx:analysis_result[@analysis='asapratio_pvalue'])">
              <xsl:text>&#9;</xsl:text>
              <xsl:text>&#9;</xsl:text>
              <xsl:text>&#9;</xsl:text>
            </xsl:if>
          </xsl:if>

          <!-- XPRESS -->
          <xsl:if test="$xpress">
            <xsl:if test="../protx:analysis_result[@analysis='xpress']">
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'xpress']/protx:XPressRatio/@ratio_mean"/>
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'xpress']/protx:XPressRatio/@ratio_standard_dev"/>
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'xpress']/protx:XPressRatio/@ratio_number_peptides"
                  />
            </xsl:if>
            <xsl:if test="not(../protx:analysis_result[@analysis='xpress'])">
              <xsl:text>&#9;</xsl:text>
              <xsl:text>&#9;</xsl:text>
              <xsl:text>&#9;</xsl:text>
            </xsl:if>
          </xsl:if>

          <!-- Libra -->
          <xsl:if test="$libra">
            <xsl:if test="../protx:analysis_result[@analysis='libra']">
              <xsl:for-each
                  select="../protx:analysis_result[@analysis = 'libra']/protx:libra_result/protx:intensity">
		<xsl:text>&#9;</xsl:text>
		<xsl:value-of select="../@mz"/>
		<xsl:text>&#9;</xsl:text>
		<xsl:value-of select="../@ratio"/>
		<xsl:text>&#9;</xsl:text>
		<xsl:value-of select="../@error"/>
              </xsl:for-each>
            </xsl:if>
            <xsl:if test="not(../protx:analysis_result[@analysis='libra'])">
              <xsl:for-each
                  select="../protx:analysis_result[@analysis = 'libra']/protx:libra_result/protx:intensity">
		<xsl:text>&#9;</xsl:text>
		<xsl:text>&#9;</xsl:text>
		<xsl:text>&#9;</xsl:text>
              </xsl:for-each>
            </xsl:if>
          </xsl:if>

	  <!-- StPeter -->
          <xsl:if test="$stpeter">
            <xsl:if test="../protx:analysis_result[@analysis='stpeter']">
              <xsl:text>&#9;</xsl:text>
              <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'stpeter']/protx:StPeterQuant/@SIn"/>
              <xsl:text>&#9;</xsl:text>
	      <xsl:value-of
                  select="../protx:analysis_result[@analysis = 'stpeter']/protx:StPeterQuant/@counts"/>
            </xsl:if>
            <xsl:if test="not(../protx:analysis_result[@analysis='stpeter'])">
              <xsl:text>&#9;</xsl:text>
            </xsl:if>
          </xsl:if>       

	  
          <!-- Peptides for this protein -->

          <xsl:if test="not($iprophet)">
            <xsl:text>&#9;</xsl:text>
	    <xsl:if test="not(protx:modification_info)">
	      <xsl:value-of select="@peptide_sequence"/>
	    </xsl:if>
	    <xsl:if test="protx:modification_info">
	      <xsl:value-of select="protx:modification_info/@modified_peptide"/>
	    </xsl:if>
            <xsl:text>&#9;</xsl:text>
            <xsl:value-of select="@is_nondegenerate_evidence"/>
	  </xsl:if>
          <xsl:if test="$iprophet">
	    <xsl:for-each select="protx:indistinguishable_peptide">
	      <xsl:text>&#9;</xsl:text>
	      <xsl:if test="not(protx:modification_info)">
		<xsl:value-of select="@peptide_sequence"/>
	      </xsl:if>
	      <xsl:if test="protx:modification_info">
		<xsl:value-of select="protx:modification_info/@modified_peptide"/>
	      </xsl:if>
	      <xsl:text>&#9;</xsl:text>
	      <xsl:value-of select="../@is_nondegenerate_evidence"/>
            </xsl:for-each>
          </xsl:if>
	  <xsl:text>
</xsl:text>
	  <!-- For each peptide in the protein -->
	</xsl:for-each> 
	<!-- Only if protein probability is greater equal  threshold -->
      </xsl:if>
      	<!-- For each protein in the group -->
      </xsl:for-each>
    </xsl:template>

  </xsl:stylesheet>
