# Maize RNA-seq Analysis Pipeline

A collaborative RNA-seq analysis workflow for maize genetics.

## Workflow Overview

```mermaid
graph TD
    subgraph "Indexing"
    A[STAR Index Genome] --> |indexed genome|B
    C[Kallisto Index Transcripts] --> |indexed transcriptome|D
    end
    
    subgraph "Mapping"
    B[STAR Align Reads] --> |aligned BAM files|E[STAR Quantification]
    end
    
    subgraph "Quantification"
    D[Kallisto Quantification] --> |transcript counts|F
    E --> |gene counts|F[Expression Analysis]
    end
    
    subgraph "Analysis"
    F --> G[Kallisto Stats]
    end
    
    subgraph "Input Files"
    H[Reference Genome] --> A
    I[Annotation GTF] --> A
    J[Reference Transcriptome] --> C
    K[RNA-seq Reads] --> B
    K --> D
    end