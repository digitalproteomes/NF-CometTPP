<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
		xmlns:xs="http://www.w3.org/2001/XMLSchema" exclude-result-prefixes="xs"
		xmlns:xd="http://www.oxygenxml.com/ns/doc/xsl" version="1.0"
		xmlns:pepx="http://regis-web.systemsbiology.net/pepXML">
  <!-- PRODUCES A TAB SEPARATED FILE FROM A PEP.XML FILE
       
       Takes one parameter to the define the peptide prophet probability cutoff.

Example: xsltproc -\-param p_threshold 0.5 <XSL_FILE> <XML_FILE>


06-09-19: - Initial commit
Patrick Pedrioli-->


  <xsl:output method="text" indent="no"/>
  <xsl:strip-space elements="*"/>

  <xsl:param name="p_threshold" select="1"/>
  
  <xsl:variable name="asap" select="boolean(pepx:msms_pipeline_analysis/pepx:msms_run_summary/pepx:spectrum_query/pepx:search_result/pepx:search_hit/pepx:analysis_result[@analysis='asapratio'])"/>
  <xsl:variable name="xpress" select="boolean(pepx:msms_pipeline_analysis/pepx:msms_run_summary/pepx:spectrum_query/pepx:search_result/pepx:search_hit/pepx:analysis_result[@analysis='xpress'])"/>
  <xsl:variable name="iprophet" select="boolean(pepx:msms_pipeline_analysis/pepx:msms_run_summary/pepx:spectrum_query/pepx:search_result/pepx:search_hit/pepx:analysis_result[@analysis='interprophet'])"/>

  <xsl:template match="/">peptide&#9;charge&#9;neutral_mass&#9;prior_aa&#9;post_aa&#9;modified_seq&#9;probability<xsl:if test="$iprophet">&#9;iprophet_probability</xsl:if><xsl:if test="$asap">&#9;asapratio_ratio&#9;asapratio_error</xsl:if><xsl:if test="$xpress">&#9;xpress_ratio&#9;xpress_light_area&#9;xpress_heavy_area</xsl:if>&#9;spectrum&#9;protein<xsl:text>
</xsl:text>
<xsl:apply-templates/>
  </xsl:template>

  <xsl:template match="pepx:search_hit[@hit_rank='1']">
    <!-- Make sure the probability passes the threshold criteria -->
    <xsl:variable name="pep_p" select="pepx:analysis_result[@analysis='peptideprophet']/pepx:peptideprophet_result/@probability"/>
    <xsl:if test="$pep_p &gt;= $p_threshold">
      <xsl:value-of select="@peptide"/>
      <xsl:text>&#9;</xsl:text><xsl:value-of select="../../@assumed_charge"/>
      <xsl:text>&#9;</xsl:text><xsl:value-of select="@calc_neutral_pep_mass"/>
      <xsl:text>&#9;</xsl:text><xsl:value-of select="@peptide_prev_aa"/>
      <xsl:text>&#9;</xsl:text><xsl:value-of select="@peptide_next_aa"/>
      <xsl:text>&#9;</xsl:text>
      <xsl:if test="not(pepx:modification_info)">
	<xsl:value-of select="@peptide"/>
      </xsl:if>
      <xsl:if test="pepx:modification_info">
	<xsl:value-of select="pepx:modification_info/@modified_peptide"/>
      </xsl:if>

      <!-- Peptide Prophet -->
      <xsl:text>&#9;</xsl:text><xsl:value-of select="$pep_p"></xsl:value-of>

      <!-- Peptide Prophet -->
      <xsl:if test="$iprophet">
	<xsl:if test="pepx:analysis_result[@analysis='interprophet']">
	  <xsl:text>&#9;</xsl:text><xsl:value-of select="pepx:analysis_result[@analysis='interprophet']/pepx:interprophet_result/@probability"></xsl:value-of>
	</xsl:if>
	<xsl:if test="not(pepx:analysis_result[@analysis='interprophet'])">
	  <xsl:text>&#9;</xsl:text>
	</xsl:if>
      </xsl:if>

      <!-- ASAPRAtio -->
      <xsl:if test="$asap">
	<xsl:if test="pepx:analysis_result[@analysis='asapratio']">
	  <xsl:text>&#9;</xsl:text><xsl:value-of select="pepx:analysis_result[@analysis='asapratio']/pepx:asapratio_result/@mean"></xsl:value-of>
	  <xsl:text>&#9;</xsl:text><xsl:value-of select="pepx:analysis_result[@analysis='asapratio']/pepx:asapratio_result/@error"></xsl:value-of> 
	</xsl:if>
	<xsl:if test="not(pepx:analysis_result[@analysis='asapratio'])">
	  <xsl:text>&#9;</xsl:text>
	  <xsl:text>&#9;</xsl:text>
	</xsl:if>
      </xsl:if>

      <!-- XPRESS -->
      <xsl:if test="$xpress">
	<xsl:if test="pepx:analysis_result[@analysis='xpress']">
	  <xsl:text>&#9;</xsl:text><xsl:value-of select="pepx:analysis_result[@analysis='xpress']/pepx:xpressratio_result/@decimal_ratio"></xsl:value-of>
	  <xsl:text>&#9;</xsl:text><xsl:value-of select="pepx:analysis_result[@analysis='xpress']/pepx:xpressratio_result/@light_area"></xsl:value-of>
	  <xsl:text>&#9;</xsl:text><xsl:value-of select="pepx:analysis_result[@analysis='xpress']/pepx:xpressratio_result/@heavy_area"></xsl:value-of>
	</xsl:if>
	<xsl:if test="not(pepx:analysis_result[@analysis='xpress'])">
	  <xsl:text>&#9;</xsl:text>
	  <xsl:text>&#9;</xsl:text>
	  <xsl:text>&#9;</xsl:text>
	</xsl:if>
      </xsl:if>
      <xsl:text>&#9;</xsl:text><xsl:value-of select="../../@spectrum"/>

      <xsl:text>&#9;</xsl:text><xsl:value-of select="@protein"/>
      <xsl:for-each select="pepx:alternative_protein">
	<xsl:text>,</xsl:text><xsl:value-of select="@protein"/>
      </xsl:for-each>
      
      <xsl:text>
</xsl:text>
    </xsl:if>
  </xsl:template>

</xsl:stylesheet>

