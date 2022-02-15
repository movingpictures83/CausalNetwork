# CausalNetwork
# Language: R
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 1.1, R 4.0.0
# Dependency: bnlearn 4.7

PluMA plugin that computes a causal network (Sazal et al, 2019).

The plugin takes as input a tab-delimited TXT file of keyword-value pairs:

undirected: Undirected network (CSV)
directed: Directed network (CSV)
datafile: Input data (CSV)
bootfile: Bootstrap (CSV)
group1: List of Group 1 data (line-by-line)
group2: List of Group 2 data (line-by-line)
graphid: Integer identifier for network


Two files are produced using the user-specified output PREFIX:
prefix.correlation.csv: Correlation network
prefix.xgmml: Causal network 
