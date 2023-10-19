GC_RANGE = (0, 100)
STANDART_LENGHT = (0, 2**32)
THRESHOLD_VALUE = (0, 41)


def check_conditions(*conditions) -> list:
    """
    Checks if the given conditions are met and returns a list of condition tuples.

    Args:
        *conditions: Variable number of conditions to check.

    Returns:
        list: List of condition tuples, where each tuple contains the condition
        and a value representing its success.
    """
    condition_list = []
    for condition in conditions:
        if type(condition) is int:
            if condition == conditions[-1]:
                condition_list.append((condition, 41))
            else:
                condition_list.append((0, condition))
        else:
            condition_list.append(condition)
    return condition_list


def is_length_in_borders(seq_length: int, length_bounds: tuple) -> bool:
    """
    Checks if the sequence length falls within the specified range.

    Args:
        seq_length (int): Length of the sequence.
        length_bounds (tuple): Minimum and maximum boundaries for the sequence length.

    Returns:
        bool: True if the length falls within the range, False otherwise.
    """
    
    if length_bounds[0] <= seq_length <= length_bounds[1]:
        return True
    return False


def gc_test(seq: str, seq_length: int, gc_bounds: tuple) -> bool:
    """
    Checks if the GC ratio of the sequence falls within the specified range.

    Args:
        seq (str): Sequence string.
        seq_length (int): Length of the sequence.
        gc_bounds (tuple): Minimum and maximum boundaries for the GC ratio.

    Returns:
        bool: True if the GC ratio falls within the range, False otherwise.
    """
    
    gc_count = seq.count("G") + seq.count("C")
    gc_ratio = gc_count / seq_length * 100
    if gc_bounds[0] <= gc_ratio <= gc_bounds[1]:
        return True
    return False


def quality_test(quality_info: str, quality_threshold: tuple) -> bool:
    """
    Checks if the quality value of the sequence falls within the specified range.

    Args:
        quality_info (str): String representing the quality values of each nucleotide in the sequence.
        quality_threshold (tuple): Minimum and maximum boundaries for the quality value.

    Returns:
        bool: True if the quality value falls within the range, False otherwise.
    """
    
    quality_each_nucl = [ord(i) - 33 for i in list(quality_info)]
    quality_value = sum(quality_each_nucl) / len(quality_each_nucl)
    if quality_threshold[0] <= quality_value <= quality_threshold[1]:
        return True
    return False


def read_input_file(input_path: str) -> dict:
    """
    Reads the input file and returns a dictionary containing the sequences and their corresponding qualities.

    Args:
        input_path (str): Path to the input file.

    Returns:
        dict: Dictionary containing the sequences as keys and their corresponding qualities as values.
    """
    
    fastq_dict = {}
    with open(input_path, 'r') as f:
        fastq_lines = f.readlines()
    for index, element in enumerate(fastq_lines, 0):
        if index % 4 == 0:
            seq_name = element[0:-1]
            fastq_dict[seq_name] = (fastq_lines[index+1][0:-1], fastq_lines[index+3][0:-1])
    return fastq_dict


def write_output_file(
    input_path: str,
    output_filename: str or None,
    filtered_seqs_dict: dict) -> None:
    """
    Writes the filtered sequences to an output file.

    Args:
        input_path (str): The path to the input file.
        output_filename (str or None): The name of the output file. If None, overwrites the input file.
        filtered_seqs_dict (dict): A dictionary containing the filtered sequences.

    Returns:
        None
    """
    
    if output_filename is None:
        output_path = input_path
    else:
        output_path = input_path[0: input_path.rfind('/')]  + '/' + output_filename
    with open(output_path, 'w') as f:
        for name, (seq, quality_info) in filtered_seqs_dict.items():
            output_list = [name, seq, '+' + name[1:], quality_info]    #наверное лучше сделать отдельную переменную для третьей строеки
            f.write('\n'.join(output_list))
    return


def fastq_analysis(
    input_path: str,
    output_filename: str|None = None,
    gc_bounds: float|tuple = GC_RANGE,
    length_bounds: float|tuple = STANDART_LENGHT,
    quality_threshold: int|tuple = THRESHOLD_VALUE,
) -> None:
    """
    Performs fastq analysis on a given input file.

    Args:
        input_path (str): The path to the input file.
        output_filename (str or None): The name of the output file. If None, overwrites the input file.
        gc_bounds (float or tuple): A tuple specifying the GC content range for filtering sequences. Default is GC_RANGE.
        length_bounds (float or tuple): A tuple specifying the length range for filtering sequences. Default is STANDART_LENGHT.
        quality_threshold (int or tuple): A tuple specifying the quality threshold for filtering sequences. Default is THRESHOLD_VALUE.

    Returns:
        None
    """

    seqs = read_input_file(input_path)
    gc_bounds, length_bounds, quality_threshold = check_conditions(
        gc_bounds, length_bounds, quality_threshold
    )
    filtered_seqs_dict = {}
    for name, (seq, quality_info) in seqs.items():
        seq_length = len(seq)
        test_borders_result = is_length_in_borders(seq_length, length_bounds)
        test_gc_result = gc_test(seq, seq_length, gc_bounds)
        test_quality_result = quality_test(quality_info, quality_threshold)
        if all([test_borders_result, test_gc_result, test_quality_result]):
            filtered_seqs_dict[name] = (seq, quality_info)
    write_output_file(input_path, output_filename, filtered_seqs_dict)
    return filtered_seqs_dict