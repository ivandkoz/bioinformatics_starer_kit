# Основная функция. Распределяет задачи по остальным
# Также проверяет наличие "верных" последовательностей,
# в случае отсутсвия выводит сообщение об ошибке

def run_dna_rna_tools(*task):
    seqs = task[:-1]
    operation = task[-1]
    tested_seqs = check_dna_rna(*seqs)
    if not tested_seqs:
        print("All yours sequences aren't DNA or RNA, please check your input")
        return
    result = eval(operation + '(*tested_seqs)')
    if len(result) == 1:
        return result[0]
    elif not result:
        print("Nothing to modify, please check your input")
        return
    else:
        return result

# Выполняет функцию разворота последовательности


def reverse(*seqs):
    reversed_seqs = [seq[::-1] for seq in seqs]
    return reversed_seqs

# Транскрибирует последовательность


def transcribe(*seqs):
    transcribed_seqs = []
    complement_dict = {'a': 'a', 'c': 'c',
                       't': 'u', 'g': 'g',
                       'A': 'A', 'C': 'C',
                       'T': 'U', 'G': 'G'}
    for seq in seqs:
        if check_dna(seq) is True:
            seq = ''.join([complement_dict.get(letter) for letter in seq])
            transcribed_seqs.append(seq)
        else:
            print(f'{seq} is RNA, transcription failed')
    return transcribed_seqs

# Создает последовательность комплементарную заданной


def complement(*seqs):
    complement_seqs = []
    complement_dict_dna = {'a': 't', 'c': 'g',
                           't': 'a', 'g': 'c',
                           'A': 'T', 'C': 'G',
                           'T': 'A', 'G': 'C'}
    complement_dict_rna = {'a': 'u', 'c': 'g',
                           'u': 'a', 'g': 'c',
                           'A': 'U', 'C': 'G',
                           'U': 'A', 'G': 'C'}
    for seq in seqs:
        if check_dna(seq) is True:
            complement_dict = complement_dict_dna
        else:
            complement_dict = complement_dict_rna
        seq = ''.join([complement_dict.get(letter) for letter in seq])
        complement_seqs.append(seq)
    return complement_seqs

# Объединяет две функции reverse и complement


def reverse_complement(*seqs):
    reversed_complement_seqs = complement(*reverse(*seqs))
    return reversed_complement_seqs

# Проверка является ли последовательность ДНК или РНК
# Если нет, удаляет из очереди и выводит сообщение об ошибке


def check_dna_rna(*seqs):
    checked_seqs = []
    nucleotids = ['A', 'a', 'T', 't',
                  'C', 'c', 'G', 'g',
                  'U', 'u']
    for seq in seqs:
        test_nucl = all(letter in nucleotids for letter in seq)
        test_t_u = any(letter in ['u', 'U']
                       for letter in seq) and any(letter in ['t', 'T']
                                                  for letter in seq)
        if test_t_u is True or test_nucl is False:
            print(f'Error {seq} is not DNA or RNA')
            continue
        checked_seqs.append(seq)
    return checked_seqs

# Узнает ДНК или РНК. Нужно для транскрипции


def check_dna(seq):
    result = not any(letter in seq for letter in ['u', 'U'])
    return result

# Считает длину каждой последовательности


def lenght(*seqs):
    result = [len(seq) for seq in seqs]
    return result

# Считает содержание G и С в каждой последовательности. Вывод в процентах


def gc_content(*seqs):
    gc_index_list = []
    gc_letters = ['G', 'g', 'C', 'c']
    for seq in seqs:
        gc_count = sum(seq.count(letter) for letter in gc_letters)
        len_seq = len(seq)
        index = round(gc_count / len_seq * 100, 2)
        gc_index_list.append(index)
    return gc_index_list

# Считает длину каждой последовательности


def double_stranded(*seqs):
    seqs_2 = complement(*seqs)
    i = 0
    for seq in seqs:
        print('', seq, '\n',
              ''.join(['|' for i in range(len(seq))]), '\n',
              seqs_2[i])
        i += 1
        if check_dna(seq) is False:
            print(f'{seq} is RNA, be careful, double stranded RNA is rare')
    return seqs_2
