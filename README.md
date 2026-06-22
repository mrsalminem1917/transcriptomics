# Bulk Transcriptomics Pipeline: Young vs Old Analysis

Индивидуальное задание по анализу данных bulk RNA-Seq. Проект включает в себя автоматизированный Bash-пайплайн первичной обработки данных, а также комплексный статистический анализ дифференциальной экспрессии генов (DEG) и функционального обогащения (GSEA).

---

## 📌 Содержание
- [Структура репозитория](#-структура-репозитория)
- [Используемый стек технологий](#-используемый-стек-технологий)
- [Архитектура пайплайна (pipeline.sh)](#-архитектура-пайплайна-pipelinesh)

---
## 📂 Структура репозитория

В репозитории содержатся ключевые скрипты для каждого этапа биоинформатического анализа:

*   [`pipeline.sh`](https://github.com/mrsalminem1917/transcriptomics/blob/main/pipeline.sh) — Bash-скрипт для сквозного параллельного процессинга сырых чтений.
*   [`FeatureCounts_DESeq2`](https://github.com/mrsalminem1917/transcriptomics/blob/main/FeatureCounts_DESeq2) — скрипт анализа дифференциальной экспрессии генов на уровне экзонов после STAR-выравнивания с фильтрацией батч-эффектов.
*   [`salmon_deg_script.R`](https://github.com/mrsalminem1917/transcriptomics/blob/main/salmon_deg_script.R) — альтернативный пайплайн квантификации транскриптов с помощью псевдовыравнивания Salmon и DESeq2.
*   [`gsea.R`](https://github.com/mrsalminem1917/transcriptomics/blob/main/gsea.R) — функциональный анализ обогащения генных сетей (GSEA) на основе Hallmark-сигнатур из базы MSigDB.

---

## 🛠 Используемый стек технологий

### Инструменты процессинга (Bash/Linux):
*   **Загрузка данных:** `SRA Toolkit (fastq-dump)`
*   **Контроль качества:** `FastQC`, `QoRTs`, `MultiQC`
*   **Картирование и квантификация:** `STAR` (с общим разделяемым индексом в RAM), `featureCounts`, `Salmon`
*   **Параллелизация:** `GNU Parallel` (оптимальное распределение CPU между датасетами)

### Статистический анализ (R / Bioconductor):
*   **Нормализация и DEG:** `DESeq2` (метод RLE/MRN), `edgeR` (фильтрация низкоэкспрессируемых генов)
*   **Коррекция батч-эффектов:** `sva (svaseq)` для поиска скрытых переменных, `limma (removeBatchEffect)` для визуализации
*   **Обогащение:** `fgsea`, `msigdbr`, `org.Hs.eg.db`, `decoupler`
*   **Визуализация:** `ggplot2`, `pheatmap`, `EnhancedVolcano`, `ggvenn`, `patchwork`

---

