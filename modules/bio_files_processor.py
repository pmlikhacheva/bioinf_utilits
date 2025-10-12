"""
BioFiles Processor - Скрипт для обработки биоинформатических файлов.
"""

import os

def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    """
    Конвертирует multiline FASTA в oneline FASTA.
    
    Args:
        input_fasta: Путь к входному FASTA файлу.
        output_fasta: Путь для сохранения результата (опционально).
    
    Returns:
        Путь к созданному файлу.
    
    """
    
    if output_fasta is None:
        base_name = os.path.splitext(input_fasta)[0]
        output_fasta = f"{base_name}_oneline.fasta"
    
    seqs = []
    curr_head = ""
    curr_seq = ""
    
    with open(input_fasta, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # Сохраняем предыдущую последовательность если есть
                if curr_header and curr_seq:
                    seqs.append((curr_header, curr_sequence))
                # Начинаем новую последовательность
                curr_header = line
                curr_seq = ""
            else:
                # Добавляем к текущей последовательности
                curr_seq += line
    
    # Добавляем последнюю последовательность
    if curr_header and curr_seq:
        seqs.append((curr_header, curr_seq))
    
    # Записываем в выходной файл
    with open(output_fasta, 'w') as file:
        for header, sequ in seqs:
            file.write(f"{header}\n")
            file.write(f"{seq}\n")

    # Выводим комментарий, как метку того, что программа корректно отработала
    print(f"Конвертировано {len(seqs)} последовательностей в {output_fasta}")
    return output_fasta

