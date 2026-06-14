#!/bin/bash

# 1. Проверка входящих аргументов
if [ $# -eq 0 ]; then
    echo "Ошибка: Не указаны датасеты для обработки!"
    echo "Использование: $0 DATASET1 [DATASET2 DATASET3 ...]"
    exit 1
fi

# ---- НАСТРОЙКА ПАРАЛЛЕЛИЗМА ----
# Сколько датасетов обрабатывать ОДНОВРЕМЕННО (каждый берет до 10 CPU)
CONCURRENT_JOBS=2 
# --------------------------------

# Исправленные пути через $HOME (безопасно для кавычек)
export BASE_DIR="$HOME/storage/AnTrD_personal_task"
export REF_DIR="$BASE_DIR/reference"
export LOGS_DIR="$BASE_DIR/logs"

# Создаем центральную папку для логов заранее
mkdir -p "$LOGS_DIR"

# Инициализируем массив датасетов из аргументов командной строки
DATASETS=("$@")

# 2. Функция обработки ОДНОГО датасета
process_dataset() {
    local DATASET=$1
    echo "[$(date +%T)] >>> НАЧАЛО: $DATASET"

    # Шаг 1: Скачивание данных
        # Шаг 1: Скачивание данных в папку raw_data
    fastq-dump --split-3 --gzip --skip-technical -O "$BASE_DIR/raw_data" "$DATASET" > "$LOGS_DIR/${DATASET}_download.log" 2>&1

    # Шаг 2: Контроль качества FastQC
    mkdir -p "$BASE_DIR/fastqc_reports"
    fastqc -t 4 -o "$BASE_DIR/fastqc_reports" "$BASE_DIR/raw_data/${DATASET}.fastq.gz" > "$LOGS_DIR/${DATASET}_fastqc.log" 2>&1

    # Шаг 3: Выравнивание STAR
    mkdir -p "$BASE_DIR/star_alignment/${DATASET}_Aligned/"
    STAR --runThreadN 10 \
        --genomeDir "$REF_DIR/star_index" \
        --readFilesIn "$BASE_DIR/raw_data/${DATASET}.fastq.gz" \
        --readFilesCommand zcat \
        --genomeLoad LoadAndKeep \
        --outFileNamePrefix "$BASE_DIR/star_alignment/${DATASET}_Aligned/${DATASET}_" \
        --outSAMtype BAM Unsorted \
        --quantMode GeneCounts > "$LOGS_DIR/${DATASET}_star_mapping.log" 2>&1

    # Шаг 4: Сортировка и индексация BAM
    samtools sort -@ 10 \
        "$BASE_DIR/star_alignment/${DATASET}_Aligned/${DATASET}_Aligned.out.bam" \
        -o "$BASE_DIR/star_alignment/${DATASET}_Aligned/${DATASET}_Aligned.out_sorted.bam" > "$LOGS_DIR/${DATASET}_samtools_sort.log" 2>&1
        
    samtools index -@ 10 \
        "$BASE_DIR/star_alignment/${DATASET}_Aligned/${DATASET}_Aligned.out_sorted.bam" > "$LOGS_DIR/${DATASET}_samtools_index.log" 2>&1

    # Шаг 5: Контроль качества QoRTs
    mkdir -p "$BASE_DIR/qorts_reports/${DATASET}"
    java -Xmx4G -jar "$BASE_DIR/QoRTs-STABLE.jar" QC \
        --generatePlots \
        --singleEnded \
        "$BASE_DIR/star_alignment/${DATASET}_Aligned/${DATASET}_Aligned.out_sorted.bam" \
        "$REF_DIR/gencode.v49.primary_assembly.annotation.gtf" \
        "$BASE_DIR/qorts_reports/${DATASET}" > "$LOGS_DIR/${DATASET}_qorts.log" 2>&1

    # Шаг 6: Подсчет чтений featureCounts
    mkdir -p "$BASE_DIR/feature_counts/gene_level"
    
    # 6.1 Уровень экзонов
    featureCounts -s 0 -T 4 \
        -a "$REF_DIR/gencode.v49.primary_assembly.annotation.gtf" \
        -t exon -g gene_id \
        -o "$BASE_DIR/feature_counts/${DATASET}_counts.txt" \
        "$BASE_DIR/star_alignment/${DATASET}_Aligned/${DATASET}_Aligned.out_sorted.bam" > "$LOGS_DIR/${DATASET}_exon_featureCounts.log" 2>&1
    
    # 6.2 Уровень генов
    featureCounts -s 0 -M -T 4 \
        -a "$REF_DIR/gencode.v49.primary_assembly.annotation.gtf" \
        -t gene -g gene_id \
        -o "$BASE_DIR/feature_counts/gene_level/${DATASET}_counts.txt" \
        "$BASE_DIR/star_alignment/${DATASET}_Aligned/${DATASET}_Aligned.out_sorted.bam" > "$LOGS_DIR/${DATASET}_gene_featureCounts.log" 2>&1
    
    # Шаг 7: Квантификация Salmon
    mkdir -p "$BASE_DIR/salmon/transcripts_quant"
    salmon quant -i "$BASE_DIR/salmon/salmon_index" -l A \
        -r "$BASE_DIR/raw_data/${DATASET}.fastq.gz" \
        -o "$BASE_DIR/salmon/transcripts_quant/${DATASET}" \
        --validateMappings --useVBOpt --seqBias -p 4 > "$LOGS_DIR/${DATASET}_salmon_quant.log" 2>&1

    echo "[$(date +%T)] <<< ФИНИШ: $DATASET"
}

# Экспортируем функцию для GNU Parallel
export -f process_dataset

# 3. Загрузка общего индекса STAR в оперативную память
echo "Загрузка индекса STAR в оперативную память..."
STAR --genomeLoad LoadAndExit --genomeDir "$REF_DIR/star_index" > "$LOGS_DIR/star_genomeLoad_init.log" 2>&1

# 4. Параллельный запуск функции через GNU Parallel
echo "Запуск параллельной обработки через GNU Parallel (одновременно задач: $CONCURRENT_JOBS)..."
parallel -j "$CONCURRENT_JOBS" process_dataset ::: "${DATASETS[@]}"

# 5. Выгрузка индекса STAR из памяти
echo "Выгрузка индекса STAR из оперативной памяти..."
STAR --genomeLoad Remove --genomeDir "$REF_DIR/star_index" > "$LOGS_DIR/star_genomeLoad_remove.log" 2>&1

# 6. Сбор раздельных отчетов для каждого этапа (из папок с результатами)
echo "Запуск раздельного сбора MultiQC..."

# Создаем отдельные папки для отчетов MultiQC
mkdir -p "$BASE_DIR/multiqc_reports/fastqc"
mkdir -p "$BASE_DIR/multiqc_reports/star"
mkdir -p "$BASE_DIR/multiqc_reports/salmon"
mkdir -p "$BASE_DIR/multiqc_reports/feature_counts_exon"
mkdir -p "$BASE_DIR/multiqc_reports/feature_counts_gene"
mkdir -p "$BASE_DIR/multiqc_reports/qorts"

# MultiQC сканирует оригинальные файлы результатов (.html, .summary, .out, .json)
/opt/conda/bin/multiqc "$BASE_DIR/fastqc_reports" -o "$BASE_DIR/multiqc_reports/fastqc" -n fastqc_report
/opt/conda/bin/multiqc "$BASE_DIR/star_alignment" -o "$BASE_DIR/multiqc_reports/star" -n star_report
/opt/conda/bin/multiqc "$BASE_DIR/salmon"         -o "$BASE_DIR/multiqc_reports/salmon" -n salmon_report
/opt/conda/bin/multiqc "$BASE_DIR/qorts_reports"  -o "$BASE_DIR/multiqc_reports/qorts" -n qorts_report

# Разделяем сбор featureCounts: игнорируем gene_level для экзонов, и смотрим строго в gene_level для генов
/opt/conda/bin/multiqc "$BASE_DIR/feature_counts" --ignore gene_level -o "$BASE_DIR/multiqc_reports/feature_counts_exon" -n feature_counts_exon_report
/opt/conda/bin/multiqc "$BASE_DIR/feature_counts/gene_level" -o "$BASE_DIR/multiqc_reports/feature_counts_gene" -n feature_counts_gene_report

echo "Конвейер успешно завершен!"
echo "Все консольные логи (stdout/stderr) сохранены в: $LOGS_DIR"