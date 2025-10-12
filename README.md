---
editor_options: 
  markdown: 
    wrap: 72
---

# BioUtils

Биоинформатические утилиты для работы с ДНК/РНК последовательностями и
FASTQ файлами.

## Установка

\`\`\`bash git clone
<https://github.com/pmlikhacheva/bioinf_utilits/bioutils.git> cd
bioutils

## Описание

### Основные функции

#### DNA/RNA Tools

run_dna_rna_tools(\*args)

Основная функция для работы с последовательностями нуклеиновых кислот.

Параметры:

\*args: Произвольное количество аргументов. Последний аргумент -
название операции, остальные - последовательности ДНК/РНК. Возвращает:

str или list: Результат обработки последовательностей

#### Отдельные DNA/RNA функции

python from modules.dna_tools import correct_nucleic_acid, transcribe,
reverse_sequence, complement_sequence

# Проверка валидности последовательности

is_valid = correct_nucleic_acid("ATCG") \# True

# Транскрипция ДНК в РНК

rna_seq = transcribe("ATCG") \# "AUCG"

# Разворот последовательности

reversed_seq = reverse_sequence("ATCG") \# "GCTA"

# Комплементарная последовательность

comp_seq = complement_sequence("ATCG") \# "TAGC"

#### FASTQ обработка

filter_fastq(input_fastq, output_fastq, gc_bounds, length_bounds,
quality_threshold)

Фильтрует FASTQ файл по GC составу, длине и качеству.

Параметры:

input_fastq: Путь к входному FASTQ файлу output_fastq: Путь для
сохранения отфильтрованного файла gc_bounds: Интервал GC состава в % или
верхняя граница (по умолчанию: (0, 100)) length_bounds: Интервал длины
или верхняя граница (по умолчанию: (0, 2\*\*32)) quality_threshold:
Пороговое значение среднего качества (phred33) (по умолчанию: 0)
Возвращает:

dict: Отфильтрованный словарь FASTQ последовательностей Примеры
использования:

python from bioutils import filter_fastq

# Базовая фильтрация

result = filter_fastq("input.fastq", "filtered/output.fastq")

# Фильтрация по GC составу (30-70%)

result = filter_fastq("input.fastq", "filtered/output.fastq",
gc_bounds=(30, 70))

# Фильтрация по качеству (\>20)

result = filter_fastq("input.fastq", "filtered/output.fastq",
quality_threshold=20)

# Комбинированная фильтрация

result = filter_fastq("input.fastq", "filtered/output.fastq",
gc_bounds=(30, 70), length_bounds=(50, 150), quality_threshold=20) \####
Вспомогательные FASTQ функции

python from modules.fastq_tools import calculate_gc_content,
calculate_quality_score

# Расчет GC содержания

gc_content = calculate_gc_content("ATGC") \# 50.0

# Расчет среднего качества

quality_score = calculate_quality_score("IIII") \# 10.0 Обработка
биоинформатических файлов

convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None)

Конвертирует multiline FASTA в oneline FASTA.

Параметры:

input_fasta: Путь к входному FASTA файлу output_fasta: Путь для
сохранения результата (опционально) Возвращает:

str: Путь к созданному файлу Пример использования:

python from bio_files_processor import
convert_multiline_fasta_to_oneline

# С автоматическим именем файла

output = convert_multiline_fasta_to_oneline("input.fasta")

# С указанием выходного файла

output = convert_multiline_fasta_to_oneline("input.fasta",
"output_oneline.fasta") parse_blast_output(input_file, output_file)

Парсит BLAST вывод и извлекает лучшие совпадения.

Параметры:

input_file: Путь к файлу с результатами BLAST output_file: Путь для
сохранения списка белков Возвращает:

int: Количество найденных белков Пример использования:

python from bio_files_processor import parse_blast_output


