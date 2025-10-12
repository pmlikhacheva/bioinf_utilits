from modules.dna_tools import correct_nucleic_acid, transcribe, reverse_sequence, complement_sequence, reverse_complement_sequence
from modules.filter_fastq_sequences import filter_fastq_sequences, calculate_gc_content, calculate_quality_score

def run_dna_rna_tools(*args):
    """
    Основная функция для работы с последовательностями нуклеиновых кислот.
    Args:
        *args: Произвольное количество аргументов.
        Последний аргумент - название операции,
        остальные - последовательности ДНК/РНК.
    Returns:
        str or list: Результат обработки последовательностей
    """
    """
    Для того, чтобы функция отработала правильно нужно,
    чтобы на вход подавалось хотя бы два аргумента -
    одна последовательность и название опреации.
    """
    if len(args) < 2:
        return print("Нужна хотя бы одна последовательность и процедура")
    sequences = list(args[:-1])
    operation = args[-1]
    """
    Для операции is_nucleic_acid проверяем каждую последовательность
    и возвращаем результат.
    """
    if operation == 'is_nucleic_acid':
        results = [correct_nucleic_acid(seq) for seq in sequences]
        return results[0] if len(results) == 1 else results
    # Проверяем все последовательности на корректность
    for seq in sequences:
        if not correct_nucleic_acid(seq):
            return False
    # Выполняем запрошенную операцию
    if operation == 'transcribe':
        result = [transcribe(seq) for seq in sequences]
    elif operation == 'reverse':
        result = [reverse_sequence(seq) for seq in sequences]
    elif operation == 'complement':
        result = [complement_sequence(seq) for seq in sequences]
    elif operation == 'reverse_complement':
        result = [reverse_complement_sequence(seq) for seq in sequences]
    else:
        return print(f"Неизвестная операция: {operation}")
    return result[0] if len(result) == 1 else result

def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: tuple[float, float] | float = (0, 100),
    length_bounds: tuple[int, int] | int = (0, 2**32),
    quality_threshold: int = 0
) -> dict[str, tuple[str, str]]:
    """
    Фильтрует FASTQ последовательности по GC составу, длине и качеству
    
    Args:
        seqs: Словарь FASTQ последовательностей {имя: (последовательность, качество)}
        gc_bounds: Интервал GC состава в % или верхняя граница
        length_bounds: Интервал длины или верхняя граница
        quality_threshold: Пороговое значение среднего качества (phred33)
    
    Returns:
        Отфильтрованный словарь FASTQ последовательностей
    """
    return filter_fastq_sequences(seqs, gc_bounds, length_bounds, quality_threshold)
