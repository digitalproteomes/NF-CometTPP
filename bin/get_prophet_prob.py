#!/usr/bin/env python3

from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
import argparse


def get_arguments(parser):
  """Set up command line parameters
  """
  parser.add_argument("-i", "--infile",
                      help="""The MODELS.html input file.""",
                      metavar="FILE",
                      required=True)
  parser.add_argument("-e", "--error_rate",
                      help="""The desired error rate (0 to 0.1).""",
                      type=float,
                      metavar=("error_rate"),
                      default=0.01)
  parser.add_argument("-t", "--print_table",
                      help="""Print error table to stdout""",
                      action='store_true')
  parser.add_argument("-v", "--verbose",
                      help="""Print verbose messages""",
                      action='store_true')

  return parser.parse_args()


def get_table(soup, text='Error Table'):
  """Returns the table with text in the first td element
  """
  for table in soup.findAll('table'):
    if table.find('td').getText() == 'Error Table':
      return table
    

def get_generic_error_table(soup):
  """Convert content of pepxml error table to a pandas dataframe.
  Only the summary stats are extracted (not the ones for individual
  charge states.

  Returns:
  - DataFrame with error_rate; min_prob; num_correct; num_incorrect
  """
  table = get_table(soup)
  
  count = 0
  m_table = []
  for row in table.findAll('tr'):
    if row.has_attr('class'):
      if 'theader' in row.get('class'):
        count = count + 1
      if count > 1:
        break
    else:
      m_row = []
      for cell in row.findAll('td'):
        m_row.append(float(cell.getText()))
      m_table.append(m_row)

  return pd.DataFrame(m_table,
                      columns=['error_rate', 'min_prob', 'num_correct', 'num_incorrect'])


def main():
  """Main
  """    

  parser = argparse.ArgumentParser(description="""Get probability threshold for a given FDR from TPP peptide and
  protein prophet html MODELS file.""")
  args = get_arguments(parser)

  infile = args.infile
  error_rate = args.error_rate

  with open(infile) as fin:
    soup = BeautifulSoup(fin, 'html.parser')
    error_table = get_generic_error_table(soup)

    # The table does not contains all possible error_rates, let's find the closest one to the requested one.
  adjusted_error_rate_idx = error_table['error_rate'].sub(error_rate).abs().idxmin()
  adjusted_error_rate = float(error_table.loc[[adjusted_error_rate_idx]]['error_rate'])

  if args.verbose:
      if error_rate != adjusted_error_rate:
        print("Using adjusted error rate {}".format(adjusted_error_rate))
  
  print(float(error_table[error_table['error_rate'] == adjusted_error_rate]['min_prob'].min()))

  if args.print_table:
    print(error_table)
  
if __name__ == "__main__":
  main()
