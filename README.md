# HiRGC – Reference-Based Genome Compression

HiRGC is a C++ implementation of a reference-based genome compression algorithm for FASTA files.

---

## Algorithm Background

This implementation is based on the HiRGC reference-based genome compression algorithm described in the paper:

**High-speed and high-ratio referential genome compression**  
Yuansheng Liu, Hui Peng, Limsoon Wong, Jinyan Li  
*Bioinformatics*, Volume 33, Issue 21, 2017  
https://academic.oup.com/bioinformatics/article/33/21/3364/3885699

---

## Implemented Features

The implementation follows the main ideas presented in the paper, including:

- preprocessing of FASTA auxiliary information
- separation of non-ACGT symbols
- 2-bit genome encoding
- greedy matching using hash tables and k-tuples
- delta encoding
- Run-Length Encoding (RLE)
- reference-based genome reconstruction

---

## Purpose

This project was developed for research and experimentation in genomic data compression.

---

# Requirements

## Linux / Ubuntu / WSL

Install:

- `g++`
- `7zip`

## Install dependencies

```bash
sudo apt update
sudo apt install g++
sudo apt install p7zip-full
```

---

# Download

## Clone the repository

```bash 
git clone https://github.com/StellaKis/HiRGC.git
cd HiRGC
```

## Or download manually

You can also download the repository as a ZIP archive directly from GitHub and extract it locally.

---

# Compilation

Compile the program using:

```bash 
g++ HiRGC.cpp -o HiRGC
```

---

# Usage

## 1. Full compression + decompression with detailed output

Runs:
- preprocessing  
- compression  
- decompression  
- genome reconstruction  
- validation  

### Command

```bash 
./HiRGC output <target.fna> <reference.fna>
```
## 2. Compression time tracking 

Runs compression only and measures total compression time.

### Command

```bash 
./HiRGC totaltime <target.fna> <reference.fna>
```

