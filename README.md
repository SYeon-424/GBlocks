# Gblocks GUI — Python-based Gblocks-like Trimmer

A lightweight and transparent **Gblocks-style alignment trimmer**.  
Includes both:
1. **`gblock.py`** — Command-line (CLI) version  
2. **`app.py`** — Graphical (GUI) version with tooltips, plots, and official Gblocks integration

---

## OVERVIEW

This project reimplements the core idea of **Gblocks**:  
retain only well-aligned, conserved blocks from a multiple sequence alignment (MSA).

### Supported Modes
- **Internal Python trimmer** (pure Python, Gblocks-like algorithm)
- **External Gblocks executable** (calls official Gblocks binary if installed)

### Features
- Parse FASTA/TXT input  
- Auto-align unaligned sequences (MAFFT / MUSCLE / built-in center-star)  
- Trim by conservation/gap ratios  
- Export FASTA/TXT/MEGA formats  
- Visualize sequence length reduction in the GUI  

---

## I. COMMAND-LINE VERSION (`gblock.py`)

### Usage Example

```bash
python gblock.py --in input.fasta --out trimmed.fasta   --min_block_len 10 --max_gap 0.5 --min_cons 0.7 --flank_cons 0.6
```

If input is unaligned:
> MAFFT → MUSCLE → built-in aligner (in that order)

### Key Flags

| Flag | Description | Default |
|------|--------------|----------|
| `--auto_align` | Auto-align if unaligned | true |
| `--aligner` | auto / mafft / muscle / none | auto |
| `--fallback_align` | Use built-in aligner if others missing | true |
| `--min_block_len` | Minimum conserved block length | 10 |
| `--max_gap` | Max allowed gap ratio per column (0–1) | 0.5 |
| `--min_cons` | Min required conservation ratio (0–1) | 0.7 |
| `--flank_cons` | Optional relaxed edge conservation | None |
| `--max_noncons_run` | Max short nonconserved run to absorb | 0 |
| `--drop_all_gap_columns` | Remove all-gap columns | flag |

---

## Gblocks-style Integer Parameters

You can emulate **official Gblocks integer parameters**:

```bash
python gblock.py --in msa.fasta --out trimmed.fasta   --use_gblocks_params true   --gb_min_cons_count 9 --gb_flank_cons_count 14   --gb_max_noncons_run 8 --gb_min_block_len 10 --gb_allowed_gap Half
```

Automatic conversion (for N sequences):

| Parameter | Formula | Example (N=20) |
|------------|----------|----------------|
| `min_cons` | `gb_min_cons_count / N` | 9/20 = 0.45 |
| `flank_cons` | `gb_flank_cons_count / N` | 14/20 = 0.7 |
| `min_block_len` | same as Gblocks | 10 |
| `max_noncons_run` | same as Gblocks | 8 |
| `max_gap` | None=0.0, Half=0.5, All=1.0 | 0.5 |

---

## II. GUI VERSION (`app.py`)

### Launch

```bash
python app.py
```

### Features
- Scrollable **Tkinter GUI**  
- Choose between internal trimmer / external Gblocks executable  
- Auto-detect data type (DNA / Protein)  
- Tooltip explanations for all parameters  
- Real-time result graph (Before vs After length)  
- Automatic cleanup of temporary files  

### UI Layout

1. Input File (FASTA / TXT)  
2. [Optional] Gblocks Executable Path  
3. **Internal Parameters**  
   - Minimum block length  
   - Maximum gap ratio  
   - Conservation threshold  
   - Flanking edge threshold  
4. **External Parameters (Official Gblocks)**  
   - Data type (-t): p (protein), d (DNA), or auto  
   - b1–b5 parameters  
   - Gap mode: None / Half / All  
5. **Run button** — executes trimming  
6. **Result area** — output summary, file paths, and bar plot  

### Scrollable Window
The entire window supports vertical scrolling for comfortable viewing, even with long parameter lists.

---

## FILE NAMING RULES

| Mode | Example |
|------|----------|
| Internal | `input_trimmed_m10_g0.5_c0.7_f0.4.fasta` |
| External | `input_trimmed_GB_t{type}_b1{b1}_b2{b2}_b3{b3}_b4{b4}_b5{b5}.fasta` |

---

## TOOLTIP SUMMARY

### Internal Parameters
| Name | Description |
|------|-------------|
| 연속된 보존 구간 최소 길이 | 최소 블록 길이 |
| gap 허용 정도 (0~1) | 칼럼 내 갭 허용 비율 |
| 보존 비율 (0~1) | 동일 문자 비율 기준 |
| 양쪽 끝부분 품질 (0~1) | 블록 경계 품질 |

### External (Official Gblocks)
| Name | Description |
|------|-------------|
| b1 | conserved position 최소 동일 서열 수 |
| b2 | flanking position 최소 동일 서열 수 |
| b3 | 연속된 비보존 허용 수 |
| b4 | 블록 최소 길이 |
| b5 | 갭 허용 정도 (None/Half/All) |

---

## OUTPUT FILES

| File | Description |
|------|--------------|
| `.fasta` | Trimmed alignment (original labels) |
| `.txt` | Plain text equivalent |
| `_mega.fasta` | Sanitized labels for MEGA import |
| `_label_map.txt` | Original ↔ Sanitized label mapping |
| (optional) | Raw Gblocks output (external mode) |

---

## REQUIREMENTS

- Python 3.8+  
- (Optional) [MAFFT](https://mafft.cbrc.jp/alignment/software/) or [MUSCLE v5](https://www.drive5.com/muscle5/)  
- (Optional) Official Gblocks binary  
- GUI dependencies:  
  ```bash
  pip install matplotlib
  ```

---

## LICENSE

- **Code:** MIT License  
- **MAFFT, MUSCLE, and Gblocks** retain their original licenses.  
- Educational / research use encouraged with citation.

---

**End of File**
