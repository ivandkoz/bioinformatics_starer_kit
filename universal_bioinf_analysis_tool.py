from bioinf_tools.dna_rna_analisys_tools import dna_rna_analysis
from bioinf_tools.fasta_analisys_tool import fastaq_analysis
from bioinf_tools.protein_analysis_tool import protein_analysis


class Analyze:
"""
This class includes three methods. 
Each allows you to work with its own data type, similar to the name of the method"""

    def fasta(*task):
        return fasta_analysis(*task)

    def protein(*task):
        return protein_analysis(*task)

    def dna_rna(*task):
        return dna_rna_analysis(*task)
