GC_RANGE = (0, 100)
STANDART_LENGHT = (0, 2**32)
THRESHOLD_VALUE = (0, 41)


def fasta_analysis(
    seqs: dict,
    gc_bounds: tuple = GC_RANGE,
    length_bounds: tuple = STANDART_LENGHT,
    quality_threshold: tuple = THRESHOLD_VALUE,
) -> dict:
    """
    Filters values in the fast dictionary according to specified parameters

    Arguments:
      - seqs (dict), conditions values

      Return:
      - dictionary with values filtered by parameters"""

    gc_bounds, length_bounds, quality_threshold = check_conditions(
        gc_bounds, length_bounds, quality_threshold
    )
    filtered_seqs = {}
    for name in seqs.keys():
        seq, quality_info = seqs.get(name)
        seq_length = len(seq)
        test_1 = length_test(seq_length, length_bounds)
        test_2 = gc_test(seq, seq_length, gc_bounds)
        test_3 = quality_test(quality_info, quality_threshold)
        if all([test_1, test_2, test_3]):
            filtered_seqs.update({name: (seq, quality_info)})
    return filtered_seqs


def check_conditions(*conditions) -> list:
    """
    Checks if the condition is specified by a single number, then makes the condition a tuple with two boundaries

    Arguments:
      - list of conditions

      Return:
      - list of transformed condition values"""
    
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


def length_test(seq_length: int, length_bounds: tuple) -> bool:
    """
    Сhecks if the length of the sequence is within the specified bounds

    Arguments:
      - Length Conditions, sequence length

      Return:
      - True if test passed, else False"""
    
    if length_bounds[0] <= seq_length <= length_bounds[1]:
        return True
    return False


def gc_test(seq: str, seq_length: int, gc_bounds: tuple) -> bool:
    """
    Сhecks if the GC percentage of the sequence is within the specified bounds

    Arguments:
      - GC Conditions, sequence, sequence length

      Return:
      - True if test passed, else False"""
    
    gc_count = seq.count("G") + seq.count("C")
    gc_ratio = gc_count / seq_length * 100
    if gc_bounds[0] <= gc_ratio <= gc_bounds[1]:
        return True
    return False


def quality_test(quality_info: str, quality_threshold: tuple) -> bool:
    """
    Сhecks if the quality of the sequence is within the specified bounds

    Arguments:
      - Quality Conditions, sequence quality info

      Return:
      - True if test passed, else False"""
    
    quality_each_nucl = [ord(i) - 33 for i in list(quality_info)]
    quality_value = sum(quality_each_nucl) / len(quality_each_nucl)
    if quality_threshold[0] <= quality_value <= quality_threshold[1]:
        return True
    return False
  
