# 🧬 Gblocks GUI & CLI — Python-based Gblocks-like Trimmer
Lightweight, transparent reimplementation of **Gblocks** with both CLI and GUI modes.

## 🧩 Included Components
| File | Description |
|------|--------------|
| `gblock.py` | Command-line (CLI) version |
| `app.py` | GUI front-end for trimming & visualization |
| `fetch_pane.py` | NCBI Protein Fetch panel (GUI addon) |
| `main_gui.exe` | Standalone executable (GUI only, bundled with all code) |
| `mafft-win/` | MAFFT binaries (used for auto-alignment) |
| `Gblocks.exe` | Official Gblocks executable (optional) |

## 🚀 Overview
This project reimplements the core idea of **Gblocks**: retain only well-aligned, conserved blocks from a multiple sequence alignment (MSA).

Two modes are supported:
- 🧩 **Internal Python trimmer** (Gblocks-like algorithm)
- ⚙️ **External Gblocks** (calls the official Gblocks binary if present)

It can:
- Parse FASTA/TXT input  
- Auto-align unaligned sequences (MAFFT / MUSCLE / built-in aligner)  
- Trim by conservation and gap ratios  
- Export FASTA/TXT/MEGA formats  
- Visualize sequence reduction in GUI  

## 📂 Folder Structure (Final Distribution)
📁 **Program Folder**  
├─ `main_gui.exe` ← GUI executable (entry point)  
├─ `Gblocks.exe` ← Optional: official Gblocks binary  
└─ `mafft-win/`  
　　├─ `mafft.bat`  
　　├─ `mafft.exe`  
　　└─ supporting files  

💡 You can create a shortcut to `main_gui.exe` anywhere (e.g. Desktop).  
As long as `mafft-win` and `Gblocks.exe` remain in the same folder, everything works fine.

## 🧭 Quick Start (GUI)
1. **Run `main_gui.exe`** → Opens the graphical trimming interface.  
2. **Select Input File** (FASTA or TXT). Data type (DNA/protein) auto-detected.  
3. **Choose Mode**  
   - Internal Trimmer (Python)  
   - External Gblocks (requires `Gblocks.exe`)  
4. **(Optional)** Click “NCBI Fetch (Protein)” to open the integrated protein-fetch panel.  
5. **Run → Trim → Done!** Results saved automatically and summarized in the GUI.

## 🧪 Output Files
| File | Description |
|-------|-------------|
| `trimmed.fasta` | Main trimmed alignment |
| `trimmed.txt` | Plain text version |
| `trimmed_mega.fasta` | Sanitized labels for MEGA import |
| `label_map.txt` | Original ↔ sanitized label mapping |
| (External) Gblocks output | Original raw `.htm` or `.fasta` copy |

## ⚙️ Internal Parameters
| Parameter | Description |
|------------|--------------|
| **Min block length** | Minimum number of consecutive conserved positions |
| **Max gap ratio** | Allowed proportion of gaps per column |
| **Min conservation** | Fraction of identical characters required |
| **Flank conservation** | Edge quality threshold for each block |

## ⚙️ External (Official Gblocks) Parameters
| Parameter | Meaning |
|------------|----------|
| **b1–b4** | Correspond to official integer parameters |
| **b5** | Gap tolerance (None / Half / All) |
| **-t** | Data type: p = protein, d = DNA |
| **Executable** | Path auto-detected if `Gblocks.exe` exists nearby |

## 💻 CLI Version (for developers)
`python gblock.py --in input.fasta --out trimmed.fasta --min_block_len 10 --max_gap 0.5 --min_cons 0.7 --flank_cons 0.6`  

If input is unaligned:  
> MAFFT → MUSCLE → built-in aligner (in that order)

Additional options mimic Gblocks integer style:  
`python gblock.py --in msa.fasta --out trimmed.fasta --use_gblocks_params true --gb_min_cons_count 9 --gb_flank_cons_count 14 --gb_max_noncons_run 8 --gb_min_block_len 10 --gb_allowed_gap Half`

## 🧰 Requirements
- **Windows 10+**  
- **No installation required** (portable exe)  
- Python 3.8+ (for source version)  
- Optional: MAFFT / MUSCLE / Gblocks binaries  

## 🔒 License
- Code: MIT License  
- MAFFT, MUSCLE, Gblocks: their respective original licenses  
- Educational & research use encouraged (with citation)

## 🧠 Credits
- GUI Integration & Python Reimplementation: **Your build**  
- MAFFT: © Katoh & Standley  
- Gblocks: © Castresana Lab  
- Icons & UI styling inspired by modern bioinformatics suites  

---

> _Built with curiosity, caffeine, and too many test FASTAs._
