# GBlocks

Gblocks GUI & CLI — Python-based Gblocks-like Trimmer
=====================================================

A lightweight and transparent Gblocks-style alignment trimmer.
Includes both:
  1) gblock.py — Command-line (CLI) version
  2) app.py    — Graphical (GUI) version with tooltips, plots, and official Gblocks integration

───────────────────────────────────────────────

OVERVIEW
───────────────────────────────────────────────

This project reimplements the core idea of Gblocks: 
retain only well-aligned, conserved blocks from a multiple sequence alignment (MSA).

Two modes are supported:
- Internal Python trimmer (Gblocks-like algorithm)
- External Gblocks executable (calls the real Gblocks binary if installed)

It can:
• Parse FASTA/TXT input
• Auto-align unaligned sequences (MAFFT / MUSCLE / built-in center-star)
• Trim by conservation/gap ratios
• Export FASTA/TXT/MEGA formats
• Visualize sequence length reduction in the GUI

───────────────────────────────────────────────
I. COMMAND-LINE VERSION (gblock.py)
───────────────────────────────────────────────

USAGE EXAMPLES
--------------
python gblock.py --in input.fasta --out trimmed.fasta \
  --min_block_len 10 --max_gap 0.5 --min_cons 0.7 --flank_cons 0.6

If input is unaligned:
→ MAFFT → MUSCLE → built-in aligner (in that order)

Key Flags
----------
--auto_align true|false        (default: true)
--aligner auto|mafft|muscle|none
--fallback_align true|false    (default: true)
--min_block_len INT            (default: 10)
--max_gap FLOAT                (default: 0.5)
--min_cons FLOAT               (default: 0.7)
--flank_cons FLOAT|None        (default: None)
--max_noncons_run INT          (default: 0)
--drop_all_gap_columns         (flag)

───────────────────────────────────────────────
Gblocks-style Integer Parameters
───────────────────────────────────────────────
Instead of raw ratios, use integers similar to official Gblocks:

python gblock.py --in msa.fasta --out trimmed.fasta \
  --use_gblocks_params true \
  --gb_min_cons_count 9 --gb_flank_cons_count 14 \
  --gb_max_noncons_run 8 --gb_min_block_len 10 --gb_allowed_gap Half

Automatic conversion (for N sequences):
  min_cons        = gb_min_cons_count / N
  flank_cons      = gb_flank_cons_count / N
  min_block_len   = gb_min_block_len
  max_noncons_run = gb_max_noncons_run
  max_gap         = {None:0.0, Half:0.5, All:1.0}

───────────────────────────────────────────────
II. GUI VERSION (app.py)
───────────────────────────────────────────────

Launch
-------
python app.py

Features
---------
• Scrollable Tkinter GUI  
• Choose mode: Internal trimmer or official Gblocks executable  
• Auto-detect data type (DNA/protein)  
• Tooltip explanations for all parameters  
• Real-time result graph (Before vs After length)  
• Automatic cleanup of temp files  

UI Layout
----------
1. Input File (FASTA/TXT)
2. [Optional] Gblocks Executable Path
3. Internal Parameters:
   - Minimum block length
   - Maximum gap ratio
   - Conservation threshold
   - Flanking edge threshold
4. External (Gblocks) Parameters:
   - Data type (-t): p (protein), d (DNA), or auto-detect
   - b1–b5 (official Gblocks integers)
   - Gap mode: None / Half / All
5. Run button → executes trimming
6. Result area shows:
   - Output summary
   - Saved file paths (FASTA, TXT, MEGA, label map)
   - Bar plot comparing before/after sequence length

Scrollable Window
-----------------
The entire window scrolls vertically for comfortable viewing, even with long parameter lists.

Output
-------
Automatically saves:
  • _trimmed.fasta
  • _trimmed.txt
  • _trimmed_mega.fasta  (sanitized labels)
  • _label_map.txt
  • Optional: raw Gblocks output

───────────────────────────────────────────────
III. FILE NAMING RULES
───────────────────────────────────────────────
• Internal mode:
  input_trimmed_m10_g0.5_c0.7_f0.4.fasta

• External (official Gblocks) mode:
  A parameter folder named:
    GB_t{type}_b1{b1}_b2{b2}_b3{b3}_b4{b4}_b5{b5}
  Inside it: trimmed files (fasta, txt, mega, etc.)

───────────────────────────────────────────────
IV. PARAMETER TOOLTIP SUMMARY
───────────────────────────────────────────────
Internal parameters:
  - 연속된 보존 구간 최소 길이: 얼마나 긴 구간을 블록으로 인정할지
  - gap 허용 정도 (0~1): 칼럼 내 갭 허용 비율
  - 보존 비율 (0~1): 동일 문자 비율 기준
  - 양쪽 끝부분 품질 (0~1): 블록 경계 품질

External (official Gblocks):
  - b1: conserved position 최소 동일 서열 수
  - b2: flanking position 최소 동일 서열 수
  - b3: 연속된 비보존 허용 수
  - b4: 블록 최소 길이
  - b5: 갭 허용 정도 (None/Half/All)

───────────────────────────────────────────────
V. OUTPUT FILES
───────────────────────────────────────────────
- FASTA: trimmed alignment (original labels)
- TXT: plain text equivalent
- MEGA FASTA: sanitized labels for MEGA import
- Label map: original ↔ sanitized
- (external only) Gblocks output for reference

───────────────────────────────────────────────
VI. REQUIREMENTS
───────────────────────────────────────────────
- Python 3.8+
- (Optional) MAFFT or MUSCLE v5 for alignment
- (Optional) Official Gblocks executable (for external mode)
- GUI dependencies: tkinter, matplotlib

───────────────────────────────────────────────
VII. LICENSE
───────────────────────────────────────────────
- Code: MIT License
- MAFFT, MUSCLE, Gblocks retain their original licenses.
- Educational / research use encouraged with citation.

───────────────────────────────────────────────
END OF FILE
───────────────────────────────────────────────
