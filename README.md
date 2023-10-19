# bioinformatics_starer_kit
![тык](https://github.com/ivandkoz/bioinformatics_starer_kit/blob/test_input_branch/img/kto_prochital.jpg)

> А если серьезно, я заболел и еле сделал задание :( <br>
> Еще и опоздал из-за желания сделать все функции <br>
> Мне правда стыдно, я постараюсь допилить ридми

![вжух](https://github.com/ivandkoz/bioinformatics_starer_kit/blob/test_input_branch/img/vzhuh.jpg)


# Bioinformatics Starter Kit

This repository contains a collection of tools and utilities for bioinformatics analysis. It provides ready-to-use classes and functions for various tasks such as analyzing biological data, working with different file formats, and performing common bioinformatics operations.

## Installation

To use the Bioinformatics Starter Kit, you need to have Python installed on your system. You can install it by following the official Python installation guide: [Install Python](https://www.python.org/downloads/)

Once Python is installed, you can clone this repository using the following command:

```
git clone https://github.com/ivandkoz/bioinformatics_starer_kit.git
```

After cloning the repository, navigate to the project directory:

```
cd bioinformatics_starer_kit
```

## Usage

The Bioinformatics Starter Kit provides two main classes: `Analyze` and `BioFile`. These classes offer various methods for working with biological data and files.

### Analyze Class

The `Analyze` class allows you to perform analysis on different types of biological data. Here are the available methods:

1. `fasta(*args)`: Performs analysis on FASTA data.
2. `protein(*args)`: Performs analysis on protein data.
3. `dna_rna(*args)`: Performs analysis on DNA and RNA data.

To use these methods, create an instance of the `Analyze` class and call the desired method, passing any required arguments.

Example:

```python
from bioinformatics_starter_kit import Analyze

analyzer = Analyze()

# Perform FASTA analysis
result = analyzer.fasta(data)

# Perform protein analysis
result = analyzer.protein(data)

# Perform DNA/RNA analysis
result = analyzer.dna_rna(data)
```

### BioFile Class

The `BioFile` class provides functionality for working with biological files. Here are the available methods:

1. `parse_blast_output(*args)`: Extracts information from BLAST output.
2. `convert_multiline_fasta_to_oneline(*args)`: Converts a multiline FASTA file to a one-line format.
3. `select_genes_from_gbk_to_fasta(*args)`: Selects genes from a GenBank (GBK) file and converts them to FASTA format.

To use these methods, create an instance of the `BioFile` class and call the desired method, passing any required arguments.

Example:

```python
from bioinformatics_starter_kit import BioFile

file_handler = BioFile()

# Parse BLAST output
result = file_handler.parse_blast_output(file)

# Convert multiline FASTA to one-line format
result = file_handler.convert_multiline_fasta_to_oneline(file)

# Select genes from GBK file and convert to FASTA format
result = file_handler.select_genes_from_gbk_to_fasta(file)
```

## Contributing

If you find any issues or have suggestions for improvements, feel free to open an issue or submit a pull request. Contributions are welcome!

## Contacts

This project is licensed under the [MIT License](LICENSE).

```

Please, if you have any suggestions for improvement or find a bug, email ivan.d.kozin@gmail.com
