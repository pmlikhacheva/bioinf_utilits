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


def parse_blast_output(input_file, output_file):
    """
    Парсит BLAST вывод и извлекает лучшие совпадения
    
    Args:
        input_file: Путь к файлу с результатами BLAST
        output_file: Путь для сохранения списка белков
    
    Returns:
        Количество найденных белков
    """
    
    proteins = set()
    
    with open(input_file, 'r') as file:
        content = file.read()
    
    # Разделяем на секции по каждому запросу
    query_sections = re.split(r'Query= ', content)
    
    for section in query_sections[1:]:  # Пропускаем первую пустую секцию
        # Ищем секцию с значимыми совпадениями
        if 'Sequences producing significant alignments:' in section:
            # Извлекаем часть после заголовка
            alignments_part = section.split('Sequences producing significant alignments:')[1]
            # Берем первую строку с описанием белка
            lines = alignments_part.split('\n')
            for line in lines[1:]:  # Пропускаем строку с заголовками
                line = line.strip()
                if line and not line.startswith('>') and line:
                    # Извлекаем описание белка (все после e-value)
                    parts = line.split()
                    if len(parts) >= 3:
                        # Объединяем все части кроме первых двух (score и e-value)
                        protein_name = ' '.join(parts[2:])
                        if protein_name and protein_name != 'Description':
                            proteins.add(protein_name)
                            break  # Берем только первый белок для каждого запроса
    
    # Сортируем белки по алфавиту
    sorted_proteins = sorted(proteins)
    
    # Записываем в файл
    with open(output_file, 'w') as file:
        for protein in sorted_proteins:
            file.write(f"{protein}\n")
    
    print(f"Найдено {len(sorted_proteins)} уникальных белков в {output_file}")
    return len(sorted_proteins)


