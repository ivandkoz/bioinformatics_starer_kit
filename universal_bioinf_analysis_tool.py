from bioinf_tools.dna_rna_analisys_tools import dna_rna_analysis
from bioinf_tools.fastq_analisys_tool import fastq_analysis
from bioinf_tools.protein_analysis_tool import protein_analysis
from bioinf_tools.bio_files_processor  import parse_blast_output, convert_multiline_fasta_to_oneline, select_genes_from_gbk_to_fasta, change_fasta_start_pos

class Analyze:
    """
    Class for performing analysis on different types of biological data.
    """

    def fasta(*args):
        """
        Performs analysis on FASTA data.

        Args:
            *args: Arguments for the fasta_analysis() function.

        Returns:
            The return value of the fasta_analysis() function.
        """

        return fasta_analysis(*args)

    def protein(*args):
        """
        Performs analysis on protein data.

        Args:
            *args: Arguments for the protein_analysis() function.

        Returns:
            The return value of the protein_analysis() function.
        """

        return protein_analysis(*args)

    def dna_rna(*args):
        """
        Performs analysis on DNA and RNA data.

        Args:
            *args: Arguments for the dna_rna_analysis() function.

        Returns:
            The return value of the dna_rna_analysis() function.
        """

        return dna_rna_analysis(*args)


class BioFile:
    """
    Class for working with biological files.
    """

    def parse_blast_output(*args):
        """
        Extracts information from BLAST output.

        Args:
            input_file (str): Path of the input file.
            output_file (str): Name of the output file. Default is None.

        Returns:
            The return value of the parse_blast_output() function.
        """

        return parse_blast_output(*args)

    def convert_multiline_fasta_to_oneline(*args):
        """
        Converts a FASTA file with multiline sequences to a one-line format.

        Args:
            *args: Arguments for the convert_multiline_fasta_to_oneline() function.

        Returns:
            The return value of the convert_multiline_fasta_to_oneline() function.
        """

        return convert_multiline_fasta_to_oneline(*args)

    def select_genes_from_gbk_to_fasta(*args):
        """
        Selects genes from a GenBank (GBK) file and converts them to a FASTA format.

        Args:
            *args: Arguments for the select_genes_from_gbk_to_fasta() function.

        Returns:
            The return value of the select_genes_from_gbk_to_fasta() function.
        """

        return select_genes_from_gbk_to_fasta(*args)

    def change_fasta_start_pos(*args):
        """
        Changes the starting position of a FASTA sequence.

        Args:
            *args: Arguments for the change_fasta_start_pos() function.

        Returns:
            The return value of the change_fasta_start_pos() function.
        """

        return change_fasta_start_pos(*args)

        
# Я понимаю, что это не совсем правильное написание класса, но тем не менее, работает же!
# Ну, и в моей голове идея вот такая вот, которая +- работает