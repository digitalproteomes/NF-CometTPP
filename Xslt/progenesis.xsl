<xsl:stylesheet
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"
    xmlns:pepx="http://regis-web.systemsbiology.net/pepXML">

  <!-- Copy everything not matching more specific templates-->
  <xsl:template match="node()|@*">
    <xsl:copy>
      <xsl:apply-templates select="node()|@*"/>
    </xsl:copy>
  </xsl:template>
  

  <!-- Remove any search_hit that is not top scoring one -->
  <xsl:template match="pepx:search_hit[@hit_rank!='1']"/>

  
  <!-- Remove spectrum_query with more than 1 search_hit -->
  <xsl:template match="pepx:spectrum_query[pepx:search_result[count(pepx:search_hit) > 1]]"/>

  
  <!-- Replace value of xcorr with value of PeptideProphet p-value (since Progenesis doesn't know how to properly parse pep.xml files) -->
  <xsl:template match="pepx:search_score[@name='xcorr']">
    <xsl:copy>
      <!-- Set new value of xcorr -->
      <xsl:attribute name="value">
	<xsl:value-of select="../pepx:analysis_result[@analysis='peptideprophet']/pepx:peptideprophet_result/@probability"/>
      </xsl:attribute>
      <!-- Copy other attributes unchanged -->
      <xsl:apply-templates select="@*[name(.)!='value']" />
    </xsl:copy>
  </xsl:template>

</xsl:stylesheet>
