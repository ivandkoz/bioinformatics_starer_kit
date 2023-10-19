N_BEFORE = 1
N_AFTER  = 1
OUTPUT_FASTA_NAME = 'selected_genes_from_gbk'


def read_input_file(input_path: str) -> list:
    """
    Reads the input file and returns a list of fasta lines.

    Args:
    input_path (str): Path to the input file.

    Returns:
    list: List of fasta lines.
    """
    
    fasta_dict = {}
    with open(input_path, 'r', encoding = 'utf-8') as f:
        fasta_lines = f.readlines()
    return fasta_lines


def make_fasta_dict(fasta_lines: list) -> dict:
    """
    Creates a dictionary of fasta sequences from a list of fasta lines.

    Args:
    fasta_lines (list): List of fasta lines.

    Returns:
    dict: Dictionary of fasta sequences.
    """
    
    fasta_dict = {}
    for seq in fasta_lines:
        if seq[0] == ">":
            seq_name = seq[:-1]
            fasta_dict[seq_name] = ""
        else:
            fasta_dict[seq_name] = fasta_dict[seq_name] + seq[:-1]
          # fastq_dict[seq_name] = (fastq_lines[index+1][0:-1], fastq_lines[index+3][0:-1])
    return fasta_dict


def write_output_file(
    input_path: str,
    output_filename: str or None,
    fasta_dict: dict) -> None:
    """
    Writes the fasta dictionary to an output file.
    
    Args:
    input_path (str): Path to the input file.
    output_filename (str or None): Name of the output file. If None, the output file will overwrite the input file.
    fasta_dict (dict): Dictionary of fasta sequences.
    
    Returns:
    None
    """
    
    if output_filename is None:
        output_path = input_path
    else:
        output_path = input_path[0: input_path.rfind('/')]  + '/' + output_filename + '.fasta'
    with open(output_path, 'w') as f:
        for seq_name, seq in fasta_dict.items():
            output_list = [seq_name, seq, '']    #наверное лучше сделать отдельную переменную для третьей строеки
            f.write('\n'.join(output_list))
            #print(*output_list, sep='\n', end='\n')  #Можно и так в принципе
    return


def convert_multiline_fasta_to_oneline(
    input_fasta_path: str,
    output_fasta: str) -> None:
    """
    Converts a multiline fasta file to a single line fasta file.
    
    Args:
    input_fasta_path (str): Path to the input fasta file.
    output_fasta (str): Path to the output fasta file.
    
    Returns:
    None
    """
    
    fasta_file = read_input_file(input_fasta_path)
    fasta_dict = make_fasta_dict(fasta_file)
    write_output_file(input_fasta_path, output_fasta, fasta_dict)
    return

def find_genes_and_seqs(fasta_lines: list) -> dict:
    """
    Finds genes and their corresponding sequences from a list of fasta lines.
    
    Args:
    fasta_lines (list): List of fasta lines.
    
    Returns:
    dict: Dictionary of genes and their sequences.
    """
    
    desired_seq = '/gene='
    gene_dict = {}
    for line in fasta_lines:
        if desired_seq not in line:
            continue        
        if desired_seq == '/gene=':
            gene = line[line.find(desired_seq)+7:-2]
            gene_dict[gene] = ""
            desired_seq = '/translation='
            continue
        if desired_seq == '/translation=':
            seq = line[line.find(desired_seq)+14:-1]
            if seq[-1] == '"':
                desired_seq = '/gene='
                seq = seq[:-1]
            else:
                desired_seq == ''
            gene_dict[gene] = gene_dict[gene] + seq
    return gene_dict
# Не очень мне нравится такой ститль с многим количеством условий,
# но я как-то не придумал более оптимального решения


def makes_target_genes_list(
    genes: str|list,
    n_before: int,
    n_after: int,
    gene_dict:dict) -> list:
    """
    Creates a target genes list based on a list of genes and their surrounding genes.
    
    Args:
    genes (str or list): List of genes or a single gene.
    n_before (int): Number of genes before the target gene.
    n_after (int): Number of genes after the target gene.
    gene_dict (dict): Dictionary of genes and their sequences.
    
    Returns:
    list: List of target genes and their sequences.
    """
    
    genes_in_file = list(gene_dict.keys())
    genes_target = []
    if isinstance(genes, str):    # Здесь проверяем если ген один подан в виде строки, то переделываем в список, чтобы итерироваться по нему
        genes = [genes]
    for gene in genes:
        if gene not in genes_in_file:
            raise ValueError(f'Gene {gene} not found in file')
        genes_target.extend([genes_in_file[genes_in_file.index(gene)+ num] for num in range(-n_before, n_after+1)]) # Просто берем гены с индексами от
                                                                                                                    # n_before до n_after 
                                                                                                                    # Т.к. ДНК может быть кольцевой, 
                                                                                                                    # то можно спокойно представить список в виде кольца
    return genes_target


def makes_fasta_with_target_genes(
    genes: str|list,
    n_before: int,
    n_after: int,
    gene_dict:dict) -> dict:
    """
    Constructs a fasta file containing target genes and their sequences.

    Args:
        genes (str or list): Genes to include in the fasta file.
        n_before (int): Number of characters to include before the gene sequence.
        n_after (int): Number of characters to include after the gene sequence.
        gene_dict (dict): Dictionary of genes and their sequences.

    Returns:
        dict: Dictionary of genes and their sequences in fasta format.
    """
    
    target_genes_list = makes_target_genes_list(genes, n_before,
                                                n_after, gene_dict)
    fasta_output_dict = {}
    for gene in target_genes_list:
        name = '>' + gene
        fasta_output_dict[name] = gene_dict[gene]
    return fasta_output_dict
    

def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: list|str,
    n_before: int = N_BEFORE,
    n_after: int = N_AFTER,
    output_fasta: str = OUTPUT_FASTA_NAME):
    """
    Converts selected genes from a GenBank file to a fasta file.

    Args:
        input_gbk (str): Path to the input GenBank file.
        genes (list or str): Genes to be selected.
        n_before (int, optional): The number of characters to include before the gene sequence.
        n_after (int, optional): The number of characters to include after the gene sequence.
        output_fasta (str, optional): Path for the output fasta file.

    Returns:
        None
    """
    
    file = read_input_file(input_gbk)
    gene_dict = find_genes_and_seqs(file)
    fasta_output_dict = makes_fasta_with_target_genes(genes, n_before,
                                                      n_after, gene_dict)
    write_output_file(input_gbk, output_fasta, fasta_output_dict)

    return 


def make_shift_in_each_seq(
    fasta_dict: dict,
    shift: int) -> dict:
    """
    Shifts each sequence in the fasta dictionary by a given number of characters.

    Args:
        fasta_dict (dict): Dictionary of fasta sequences.
        shift (int): The number of characters to shift each sequence.

    Returns:
        dict: Dictionary of shifted fasta sequences.
    """
    
    fasta_shifted_dict = {}
    for name, seq in fasta_dict.items():
        seq_shifted = seq[shift:] + seq[:shift]
        fasta_shifted_dict[name] = seq_shifted
    return fasta_shifted_dict


def change_fasta_start_pos(
    input_fasta: str,
    shift: int,
    output_fasta: str|None  = None) -> None:
    """
    Changes the starting position of each sequence in the fasta file.

    Args:
        input_fasta (str): Path to the input fasta file.
        shift (int): The number of characters to shift the starting position.
        output_fasta (str or None, optional): Path for the output fasta file.
            If None, the input fasta file will be overwritten.

    Returns:
        None
    """
    
    fasta_file = read_input_file(input_fasta)
    fasta_dict = make_fasta_dict(fasta_file)
    fasta_output_dict = make_shift_in_each_seq(fasta_dict, shift)
    write_output_file(input_fasta, output_fasta, fasta_output_dict)
    return


def find_names_for_queries(blast_result_file: list) -> dict:
    """
    Extracts protein names from a BLAST result file.

    Args:
        blast_result_file (list): List of lines from the BLAST result file.

    Returns:
        dict: List of protein names extracted from the BLAST result file.
    """
    
    desired_seq = 'Query #'
    protein_list = []
    for line in blast_result_file:
        if desired_seq not in line:
            continue        
        if desired_seq == 'Query #':
            desired_seq = '>'
            continue
        if desired_seq == '>':
            if '>MULTISPECIES: ' in line:
                protein_name = line[15:line.find('[')-1]
            else:
                protein_name = line[1:line.find('[')-1]
            protein_list.append(protein_name)
            desired_seq = 'Query #'
    protein_list.sort()
    return protein_list


def write_output_blast_file(
    input_path:str,
    output_filename:str,
    protein_list:list, 
    ) -> None:
    """
    Writes a list of protein names to an output file.

    Args:
        input_path (str): Path of the input file.
        output_filename (str): Name of the output file.
        protein_list (list): List of protein names.

    Returns:
        None
    """
    
    if output_filename is None:
        output_path = input_path
    else:
        output_path = input_path[0: input_path.rfind('/')]  + '/' + output_filename + '.txt'
    with open(output_path, 'w') as f:
        for protein in protein_list:
            f.write(protein+'\n')
    return
    

def parse_blast_output(
    input_file: str,
    output_file: str|None = None) -> None:
    """
    Parses a BLAST output file and writes the protein names to an output file.

    Args:
        input_file (str): Path of the input file.
        output_file (str): Name of the output file. Default is None.

    Returns:
        None
    """
    
    blast_result_file = read_input_file(input_file)
    protein_list = find_names_for_queries(blast_result_file)
    write_output_blast_file(input_file, output_file, protein_list)
    return

