#!/usr/bin/env python3

import pandas as pd
import argparse


def get_arguments(parser):
  """Set up command line parameters
  """
  parser.add_argument("-p", "--pepfile",
                      help="""The tsv exported peptide list.""",
                      metavar="FILE",
                      required=True)
  parser.add_argument("-P", "--protfile",
                      help="""The tsv exported protein list.""",
                      metavar="FILE",
                      required=True)

  parser.add_argument("-o", "--outfile",
                      help="""The output file name.""",
                      metavar="FILE",
                      required=True)
    
  return parser.parse_args()


def main():
  """Main
  """
  parser = argparse.ArgumentParser(description="""Filter the 1% FDR filtered peptide list by removing those peptides
  that have not been assigned to a protein in the 1% FDR list.""")
  args = get_arguments(parser)

  pepfile = args.pepfile
  protfile = args.protfile

  # Get the list of 1%FDR proteins
  prot_df = pd.read_csv(protfile, sep='\t', usecols=['protein'])
  protein_list =  prot_df['protein'].str.split(',').apply(pd.Series, 1).stack().unique()

  # Filter 1%FDR peptide list to only those that come from a protein
  # in the 1%FDR list
  pep_df = pd.read_csv(pepfile, sep='\t')
  filtered_pep_list = pd.DataFrame()
  first = True
  for protein in protein_list:
    if first:
      filtered_pep_list = pep_df[pep_df['protein'].str.contains(protein)]
      first = False
    filtered_pep_list = pd.concat([filtered_pep_list, pep_df[pep_df['protein'].str.contains(protein)]])
  filtered_pep_list = filtered_pep_list.drop_duplicates()

  filtered_pep_list.to_csv(args.outfile , sep='\t')


if __name__ == "__main__":
  main()
  
