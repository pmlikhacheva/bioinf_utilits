def correct_nucleic_acid(seq):
    """
    Проверяет, является ли последовательность корректной нуклеиновой кислотой.
    Args:
        seq (str): Последовательность для проверки
    Returns:
        bool: True если последовательность валидна, иначе False
    """
    if not seq:
        return False
    correct_dna_bases = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
    correct_rna_bases = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}
    seq_set = set(seq)
    """
    Проверяем, содержит ли последовательность только ДНК
    или только РНК основания
    """
    is_dna = seq_set.issubset(correct_dna_bases)
    is_rna = seq_set.issubset(correct_rna_bases)
    return bool(is_dna or is_rna)


def transcribe(seq):
    """
    Транскрибирует ДНК в РНК (заменяет T на U, t на u).
    Args:
        seq (str): ДНК последовательность
    Returns:
        str: РНК последовательность
    """
    transcription_dict = {'T': 'U', 't': 'u'}
    return ''.join(transcription_dict.get(base, base) for base in seq)


def reverse_sequence(seq):
    """
    Разворачивает последовательность.
    Args:
        seq (str): Последовательность для разворота
    Returns:
        str: Развернутая последовательность
    """
    return seq[::-1]


def complement_sequence(seq):
    """
    Возвращает комплементарную последовательность.
    Args:
        seq (str): Последовательность для комплементации
    Returns:
        str: Комплементарная последовательность
    """
    dna_complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g'
    }
    rna_complement = {
        'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
        'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'
    }
    # Определяем тип последовательности (ДНК или РНК)
    if 'U' in seq or 'u' in seq:
        complement_dict = rna_complement
    else:
        complement_dict = dna_complement
    return ''.join(complement_dict.get(base, base) for base in seq)


def reverse_complement_sequence(seq):
    """
    Возвращает обратную комплементарную последовательность.
    Args:
        seq (str): Последовательность для обработки
    Returns:
        str: Обратная комплементарная последовательность
    """
    complement = complement_sequence(seq)
    return reverse_sequence(complement)

