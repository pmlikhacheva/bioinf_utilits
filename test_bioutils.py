"""
Тестирование всех функций пакета BioUtils
"""

import os
import sys

import tempfile
sys.path.insert(0, os.path.dirname(__file__))

from bioutils import run_dna_rna_tools, filter_fastq
from modules.bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output
from modules.dna_tools import (
    correct_nucleic_acid, transcribe, reverse_sequence, 
    complement_sequence, reverse_complement_sequence, run_sequence_operation
)
from modules.fastq_tools import (
    calculate_gc_content, calculate_quality_score, read_fastq_dict,
    write_fastq_from_dict, filter_fastq_sequences, filter_fastq_on_the_fly
)


def test_dna_tools():
    """Тестирование DNA/RNA инструментов"""
    print("=== Тестирование DNA/RNA инструментов ===")
    
    # Тест correct_nucleic_acid
    assert correct_nucleic_acid("ATCG") == True, "Валидная ДНК должна возвращать True"
    assert correct_nucleic_acid("AUCG") == True, "Валидная РНК должна возвращать True"
    assert correct_nucleic_acid("ATCX") == False, "Невалидная последовательность должна возвращать False"
    assert correct_nucleic_acid("") == False, "Пустая последовательность должна возвращать False"
    print("✓ correct_nucleic_acid: OK")
    
    # Тест transcribe
    assert transcribe("ATCG") == "AUCG", "Транскрипция ДНК в РНК"
    assert transcribe("atcg") == "aucg", "Транскрипция с маленькими буквами"
    print("✓ transcribe: OK")
    
    # Тест reverse_sequence
    assert reverse_sequence("ATCG") == "GCTA", "Разворот последовательности"
    print("✓ reverse_sequence: OK")
    
    # Тест complement_sequence
    assert complement_sequence("ATCG") == "TAGC", "Комплементарная ДНК"
    assert complement_sequence("AUCG") == "UAGC", "Комплементарная РНК"
    print("✓ complement_sequence: OK")
    
    # Тест reverse_complement_sequence
    assert reverse_complement_sequence("ATCG") == "CGAT", "Реверс-комплемент ДНК"
    print("✓ reverse_complement_sequence: OK")
    
    # Тест run_sequence_operation
    assert run_sequence_operation(["ATCG"], "reverse_complement") == "CGAT", "run_sequence_operation с одной последовательностью"
    assert run_sequence_operation(["ATCG", "GCTA"], "reverse") == ["GCTA", "ATCG"], "run_sequence_operation с несколькими последовательностями"
    assert run_sequence_operation(["ATCG", "AUCG"], "is_nucleic_acid") == [True, True], "run_sequence_operation с is_nucleic_acid"
    print("✓ run_sequence_operation: OK")
    
    print("Все DNA/RNA тесты пройдены! ✅")


def test_fastq_tools():
    """Тестирование FASTQ инструментов"""
    print("\n=== Тестирование FASTQ инструментов ===")
    
    # Тест calculate_gc_content
    assert calculate_gc_content("GGCC") == 100.0, "100% GC состав"
    assert calculate_gc_content("ATAT") == 0.0, "0% GC состав"
    assert calculate_gc_content("ATGC") == 50.0, "50% GC состав"
    assert calculate_gc_content("") == 0.0, "GC состав пустой строки"
    print("✓ calculate_gc_content: OK")
    
    # Тест calculate_quality_score
    # "IIII" в phred33: (73-33)/4 = 10.0, "!!!!" в phred33: (33-33)/4 = 0.0
    assert calculate_quality_score("IIII") == 10.0, "Качество для 'IIII'"
    assert calculate_quality_score("!!!!") == 0.0, "Качество для '!!!!'"
    print("✓ calculate_quality_score: OK")
    
    # Создаем тестовый FASTQ файл
    test_fastq_content = """@read1
ATCGATCG
+
IIIIIIII
@read2
GGCCGGCC
+
AAAAAAAA
@read3
TTTTTTTT
+
!!!!!!!!
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(test_fastq_content)
        test_fastq_path = f.name
    
    try:
        # Тест read_fastq_dict
        sequences = read_fastq_dict(test_fastq_path)
        assert len(sequences) == 3, "Должно быть 3 последовательности"
        assert "read1" in sequences, "read1 должен быть в словаре"
        assert sequences["read1"] == ("ATCGATCG", "IIIIIIII"), "Корректное чтение последовательности и качества"
        print("✓ read_fastq_dict: OK")
        
        # Тест filter_fastq_sequences
        filtered = filter_fastq_sequences(sequences, gc_bounds=(30, 70), quality_threshold=5)
        assert len(filtered) == 1, "Должна остаться только read1"
        assert "read1" in filtered, "read1 должен пройти фильтрацию"
        print("✓ filter_fastq_sequences: OK")
        
        # Тест write_fastq_from_dict
        output_path = "test_output.fastq"
        write_fastq_from_dict(filtered, output_path)
        assert os.path.exists(output_path), "Выходной файл должен быть создан"
        
        # Проверяем содержимое записанного файла
        with open(output_path, 'r') as f:
            content = f.read()
            assert "@read1" in content, "read1 должен быть в выходном файле"
            assert "ATCGATCG" in content, "Последовательность должна быть в файле"
        print("✓ write_fastq_from_dict: OK")
        
        # Тест filter_fastq_on_the_fly
        on_the_fly_output = "test_on_the_fly.fastq"
        count = filter_fastq_on_the_fly(test_fastq_path, on_the_fly_output, 
                                      gc_bounds=(30, 70), quality_threshold=5)
        assert count == 1, "Должна быть сохранена 1 последовательность"
        assert os.path.exists(on_the_fly_output), "Файл должен быть создан"
        print("✓ filter_fastq_on_the_fly: OK")
        
        # Очистка
        os.unlink(output_path)
        os.unlink(on_the_fly_output)
        
    finally:
        os.unlink(test_fastq_path)
    
    print("Все FASTQ тесты пройдены! ✅")


def test_main_functions():
    """Тестирование главных функций"""
    print("\n=== Тестирование главных функций ===")
    
    # Тест run_dna_rna_tools
    result = run_dna_rna_tools("ATCG", "reverse_complement")
    assert result == "CGAT", "run_dna_rna_tools с reverse_complement"
    
    result = run_dna_rna_tools("ATCG", "AUCG", "is_nucleic_acid")
    assert result == [True, True], "run_dna_rna_tools с is_nucleic_acid для нескольких последовательностей"
    print("✓ run_dna_rna_tools: OK")
    
    # Создаем тестовый FASTQ для filter_fastq
    test_fastq_content = """@test_read1
ATCGATCG
+
IIIIIIII
@test_read2
GGCCGGCC
+
AAAAAAAA
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(test_fastq_content)
        test_input = f.name
    
    try:
        # Тест filter_fastq
        output_file = "test_filtered.fastq"
        result_dict = filter_fastq(test_input, output_file, gc_bounds=(30, 70))
        
        assert isinstance(result_dict, dict), "Должен возвращаться словарь"
        assert len(result_dict) == 1, "Должна остаться 1 последовательность"
        assert "test_read1" in result_dict, "test_read1 должен пройти фильтрацию"
        assert os.path.exists(output_file), "Выходной файл должен быть создан"
        print("✓ filter_fastq: OK")
        
        # Очистка
        os.unlink(output_file)
        
    finally:
        os.unlink(test_input)
    
    print("Все главные функции работают! ✅")


def test_bio_files_processor():
    """Тестирование функций обработки файлов"""
    print("\n=== Тестирование BioFiles Processor ===")
    
    # Тест convert_multiline_fasta_to_oneline
    multiline_fasta = """>seq1
ATCG
ATCG
>seq2
GGCC
GGCC
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(multiline_fasta)
        test_fasta = f.name
    
    try:
        output_fasta = convert_multiline_fasta_to_oneline(test_fasta)
        assert os.path.exists(output_fasta), "Выходной файл должен быть создан"
        
        # Проверяем содержимое
        with open(output_fasta, 'r') as f:
            content = f.read()
            lines = content.strip().split('\n')
            assert len(lines) == 4, "Должно быть 4 строки (2 заголовка + 2 последовательности)"
            assert lines[1] == "ATCGATCG", "Последовательность должна быть в одной строке"
            assert lines[3] == "GGCCGGCC", "Вторая последовательность должна быть в одной строке"
        print("✓ convert_multiline_fasta_to_oneline: OK")
        
        # Очистка
        os.unlink(output_fasta)
        
    finally:
        os.unlink(test_fasta)
    
    # Тест parse_blast_output (создаем тестовый файл)
    test_blast_content = """Query= sequence1

Sequences producing significant alignments:
  Description         Score    E-value
  Protein A           100      1e-50
  Protein B           90       1e-45

Query= sequence2

Sequences producing significant alignments:
  Description         Score    E-value  
  Protein C           95       1e-48
  Protein D           85       1e-40
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write(test_blast_content)
        test_blast = f.name
    
    try:
        output_proteins = "test_proteins.txt"
        count = parse_blast_output(test_blast, output_proteins)
        
        assert count == 2, "Должно быть найдено 2 белка"
        assert os.path.exists(output_proteins), "Выходной файл должен быть создан"
        
        # Проверяем содержимое
        with open(output_proteins, 'r') as f:
            proteins = [line.strip() for line in f if line.strip()]
            assert len(proteins) == 2, "Должно быть 2 белка в файле"
            assert "Protein A" in proteins, "Protein A должен быть в списке"
            assert "Protein C" in proteins, "Protein C должен быть в списке"
        print("✓ parse_blast_output: OK")
        
        # Очистка
        os.unlink(output_proteins)
        
    finally:
        os.unlink(test_blast)
    
    print("Все BioFiles Processor тесты пройдены! ✅")


def test_error_handling():
    """Тестирование обработки ошибок"""
    print("\n=== Тестирование обработки ошибок ===")
    
    # Тест на несуществующий файл
    try:
        read_fastq_dict("nonexistent_file.fastq")
        assert False, "Должна быть вызвана ошибка FileNotFoundError"
    except FileNotFoundError:
        print("✓ Обработка несуществующего файла: OK")
    
    # Тест на некорректный FASTQ формат
    bad_fastq_content = """@read1
ATCG
+
II
"""  # Неправильная длина качества
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(bad_fastq_content)
        bad_fastq = f.name
    
    try:
        try:
            read_fastq_dict(bad_fastq)
            assert False, "Должна быть вызвана ошибка ValueError"
        except ValueError:
            print("✓ Обработка некорректного FASTQ: OK")
    finally:
        os.unlink(bad_fastq)
    
    # Тест на перезапись файла
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write("@test\nATCG\n+\nIIII\n")
        existing_file = f.name
    
    try:
        test_seqs = {"test": ("ATCG", "IIII")}
        try:
            write_fastq_from_dict(test_seqs, existing_file)
            assert False, "Должна быть вызвана ошибка при перезаписи"
        except ValueError:
            print("✓ Защита от перезаписи файлов: OK")
    finally:
        os.unlink(existing_file)
    
    print("Все тесты обработки ошибок пройдены! ✅")


def run_all_tests():
    """Запуск всех тестов"""
    print("🚀 Запуск комплексного тестирования BioUtils...\n")
    
    try:
        test_dna_tools()
        test_fastq_tools()
        test_main_functions()
        test_bio_files_processor()
        test_error_handling()
        
        print("\n🎉 Все тесты успешно пройдены! Биоинформатические утилиты работают корректно.")
        print("\n📋 Итог:")
        print("  • DNA/RNA tools: ✅")
        print("  • FASTQ processing: ✅") 
        print("  • File processors: ✅")
        print("  • Error handling: ✅")
        print("  • Main functions: ✅")
        
    except Exception as e:
        print(f"\n❌ Тесты не пройдены: {e}")
        raise


if __name__ == "__main__":
    run_all_tests()
