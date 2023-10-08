COMPLEMENT_DICT = {
    "a": "a",
    "c": "c",
    "t": "u",
    "g": "g",
    "A": "A",
    "C": "C",
    "T": "U",
    "G": "G",
}
COMPLEMENT_DICT_DNA = {
    "a": "t",
    "c": "g",
    "t": "a",
    "g": "c",
    "A": "T",
    "C": "G",
    "T": "A",
    "G": "C",
}
COMPLEMENT_DICT_RNA = {
    "a": "u",
    "c": "g",
    "u": "a",
    "g": "c",
    "A": "U",
    "C": "G",
    "U": "A",
    "G": "C",
}
NUCLEOTIDS = ["A", "a", "T", "t", "C", "c", "G", "g", "U", "u"]
GC_LETTERS = ["G", "g", "C", "c"]


def dna_rna_analysis(*task):
    """
    Main function. Distributes tasks among others
    Also checks for "valid" sequences,
    if missing, displays an error message"""

    seqs = task[:-1]
    operation = task[-1]
    tested_seqs = check_dna_rna(*seqs)
    operations = {
        "reverse": reverse,
        "transcribe": transcribe,
        "complement": complement,
        "reverse_complement": reverse_complement,
        "lenght": lenght,
        "gc_content": gc_content,
        "double_stranded": double_stranded
    }

    if not tested_seqs:
        raise ValueError(
            "All yours sequences aren't DNA or RNA, please check your input"
        )
    result = operations.get(operation)(*tested_seqs)
    if len(result) == 1:
        return result[0]
    elif not result:
        raise ValueError("Nothing to modify, please check your input")
    return result


def reverse(*seqs):
    """
    Performs a sequence reversal function"""

    reversed_seqs = [seq[::-1] for seq in seqs]
    return reversed_seqs


def transcribe(*seqs):
    """
    Transcribes the sequence"""

    transcribed_seqs = []

    for seq in seqs:
        if check_dna(seq) is True:
            seq = "".join([COMPLEMENT_DICT.get(letter) for letter in seq])
            transcribed_seqs.append(seq)
        else:
            print(f"{seq} is RNA, transcription failed")
    return transcribed_seqs


def complement(*seqs):
    """
    Creates a sequence complementary to a given one"""

    complement_seqs = []

    for seq in seqs:
        if check_dna(seq) is True:
            COMPLEMENT_DICT = COMPLEMENT_DICT_DNA
        else:
            COMPLEMENT_DICT = COMPLEMENT_DICT_RNA
        seq = "".join([COMPLEMENT_DICT.get(letter) for letter in seq])
        complement_seqs.append(seq)
    return complement_seqs


def reverse_complement(*seqs):
    """
    Combines two functions reverse and complement"""

    reversed_complement_seqs = complement(*reverse(*seqs))
    return reversed_complement_seqs


def check_dna_rna(*seqs):
    """
    Checking whether the sequence is DNA or RNA.
    If not, removes it from the queue and displays an error message."""

    checked_seqs = []

    for seq in seqs:
        test_nucl = all(letter in NUCLEOTIDS for letter in seq)
        test_t_u = any(letter in ["u", "U"] for letter in seq) and any(
            letter in ["t", "T"] for letter in seq
        )
        if test_t_u is True or test_nucl is False:
            print(f"Error {seq} is not DNA or RNA")
            continue
        checked_seqs.append(seq)
    return checked_seqs


def check_dna(seq):
    """
    Recognizes DNA or RNA. Needed for transcription"""

    result = not any(letter in seq for letter in ["u", "U"])
    return result


def lenght(*seqs):
    """
    Counts the length of each sequence"""

    result = [len(seq) for seq in seqs]
    return result


def gc_content(*seqs):
    """
    Counts the G and C contents in each sequence. Output as a percentage"""
    gc_index_list = []
    for seq in seqs:
        gc_count = sum(seq.count(letter) for letter in GC_LETTERS)
        len_seq = len(seq)
        index = round(gc_count / len_seq * 100, 2)
        gc_index_list.append(index)
    return gc_index_list


def double_stranded(*seqs):
    """
    Draws a double-stranded sequence"""

    seqs_2 = complement(*seqs)
    i = 0
    for seq in seqs:
        print("", seq, "\n", "".join(["|" for i in range(len(seq))]), "\n", seqs_2[i])
        i += 1
        if check_dna(seq) is False:
            print(f"{seq} is RNA, be careful, double stranded RNA is rare")
    return seqs_2
