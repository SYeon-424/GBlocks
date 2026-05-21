"""
fetch.py
--------
NCBI 시퀀스 수집 기능 (로직 + GUI).
- 종 + 유전자명으로 NCBI에서 단백질/뉴클레오타이드 FASTA 자동 수집 → 멀티-FASTA(.txt) 저장
- NCBI 검색(db="protein" 또는 "nucleotide") → efetch로 FASTA 저장, 재시도/요청 간격 제한
- 입력: (1) species, genes / (2) 텍스트 파일(species<TAB>gene 또는 species만)
- Biopython 필요: pip install biopython

사용 예시 (코드에서)
-------------------
from fetch import fetch_from_txt, fetch_many, FetchPane

fetch_from_txt("pairs.txt", out_dir="protein_results", email="you@example.com", db="protein")
fetch_many(species, genes, out_dir="protein_results", email="you@example.com", db="protein")
"""

from __future__ import annotations

import os
import sys
import threading
import time
import ssl
from collections import defaultdict
from typing import List, Optional, Tuple

import tkinter as tk
from tkinter import filedialog, messagebox

from Bio import Entrez, SeqIO

ssl._create_default_https_context = ssl._create_unverified_context

DEFAULT_DELAY = 0.5  # 요청 간격(초)


def _ensure_email(email: Optional[str]):
    if not email or "@" not in email:
        raise ValueError("Entrez.email이 설정되지 않았습니다. 유효한 이메일을 입력하세요.")
    Entrez.email = email


def entrez_search_id(
    query: str, db: str = "protein", retmax: int = 10, retries: int = 3, delay: float = 2.0
) -> List[str]:
    for attempt in range(retries):
        try:
            with Entrez.esearch(db=db, term=query, retmax=retmax) as handle:
                rec = Entrez.read(handle)
            return rec.get("IdList", [])
        except Exception as e:
            print(f"⚠️ [검색 실패 {attempt+1}/{retries}] {query} → {e}")
            time.sleep(delay)
    print(f"❌ 검색 완전 실패: {query}")
    return []


def entrez_fetch_fasta(ncbi_id: str, db: str = "protein", retries: int = 3, delay: float = 2.0):
    """
    NCBI에서 FASTA 시퀀스를 다운로드합니다.
    """
    for attempt in range(retries):
        try:
            with Entrez.efetch(db=db, id=ncbi_id, rettype="fasta", retmode="text") as handle:
                seq_record = SeqIO.read(handle, "fasta")
            return seq_record
        except Exception as e:
            print(f"⚠️ [다운로드 실패 {attempt+1}/{retries}] ID={ncbi_id} → {e}")
            time.sleep(delay)
    print(f"❌ 다운로드 완전 실패: {ncbi_id}")
    return None


def _sanitize_filename(s: str) -> str:
    return "".join(ch for ch in s.replace(" ", "_") if ch.isalnum() or ch in "._-")


def fetch_one(
    species: str, gene: str, out_dir: str, db: str = "protein", delay: float = DEFAULT_DELAY
) -> Tuple[bool, Optional[str]]:
    """
    단일 종-유전자 조합에 대한 FASTA를 수집합니다.
    """
    os.makedirs(out_dir, exist_ok=True)
    # 미토콘드리아 등 별칭이 많은 유전자만 검색어 확장. 그 외 유전자(LEP, UCP1, BRCA1 등)는 입력한 이름으로 그대로 검색됨.
    gene_alias = {
        "COX1": ["COX1", "COI", "cytochrome c oxidase subunit I"],
        "COI": ["COI", "COX1", "cytochrome c oxidase subunit I"],
        "COX2": ["COX2", "cytochrome c oxidase subunit II"],
        "COX3": ["COX3", "cytochrome c oxidase subunit III"],
        "CYTB": ["CYTB", "cytochrome b"],
        "ND1": ["ND1", "NADH dehydrogenase subunit 1"],
        "ND2": ["ND2", "NADH dehydrogenase subunit 2"],
        "ND3": ["ND3", "NADH dehydrogenase subunit 3"],
        "ND4": ["ND4", "NADH dehydrogenase subunit 4"],
        "ND4L": ["ND4L", "NADH dehydrogenase subunit 4L"],
        "ND5": ["ND5", "NADH dehydrogenase subunit 5"],
        "ND6": ["ND6", "NADH dehydrogenase subunit 6"],
        "ATP6": ["ATP6"],
        "ATP8": ["ATP8"],
    }
    aliases = gene_alias.get(gene.upper(), [gene])
    alias_part = " OR ".join([f'"{a}"' for a in aliases])
    query = f"{species}[organism] AND ({alias_part})"
    ids = entrez_search_id(query, db=db)
    if not ids:
        print(f"🚫 {species} ({gene}) 검색 결과 없음")
        return False, None
    rec = entrez_fetch_fasta(ids[0], db=db)
    if rec is None:
        print(f"🚫 {species} ({gene}) 다운로드 실패")
        return False, None
    filename = f"{_sanitize_filename(species)}_{_sanitize_filename(gene)}.fasta"
    path = os.path.join(out_dir, filename)
    SeqIO.write(rec, path, "fasta")
    print(f"✅ 저장 완료: {path}")
    time.sleep(delay)
    return True, path


def _unique(seq: List[str]) -> List[str]:
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def fetch_many(
    species_list: List[str],
    genes: List[str],
    out_dir: str,
    email: str,
    db: str = "protein",
    combine_filename: str = "multi_fasta.txt",
    delay: float = DEFAULT_DELAY,
) -> str:
    """
    species_list의 각 종에 대해 genes 목록의 첫 번째 히트 FASTA를 받아 개별 파일로 저장하고,
    성공한 것들을 하나의 멀티-FASTA 텍스트로 병합하여 반환.
    """
    _ensure_email(email)
    os.makedirs(out_dir, exist_ok=True)
    success_paths: List[str] = []
    species_list = [s.strip() for s in species_list if s and s.strip()]
    genes = _unique([g.strip() for g in genes if g and g.strip()])
    for sp in species_list:
        for g in genes:
            ok, p = fetch_one(sp, g, out_dir=out_dir, db=db, delay=delay)
            if ok and p:
                success_paths.append(p)
    grouped = defaultdict(list)
    for fp in success_paths:
        gene = os.path.splitext(os.path.basename(fp))[0].split("_")[-1]
        grouped[gene].append(fp)
    combined = os.path.join(out_dir, combine_filename)
    for gene, paths in grouped.items():
        if len(paths) < 2:
            print(f"⚠️ {gene}: 서열이 {len(paths)}개라 Gblocks용 병합을 건너뜁니다.")
            continue
        combined = os.path.join(out_dir, f"multi_fasta_{gene}.txt")
        with open(combined, "w", encoding="utf-8") as out:
            for fp in paths:
                base = os.path.basename(fp)
                prefix = os.path.splitext(base)[0]
                with open(fp, "r", encoding="utf-8") as f:
                    lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith(">"):
                        lines[i] = f">{prefix} " + line[1:]
                        break
                out.write("".join(lines).rstrip() + "\n")
        print(f"📦 병합 저장: {combined} (entries: {len(paths)})")
    return combined


def fetch_from_txt(
    txt_path: str,
    out_dir: str,
    email: str,
    db: str = "protein",
    combine_filename: str = "multi_fasta.txt",
    delay: float = DEFAULT_DELAY,
) -> str:
    """
    txt 포맷: 'species<TAB>gene' 또는 'species' 단독 라인.
    반환: 병합된 multi_fasta.txt 경로
    """
    _ensure_email(email)
    os.makedirs(out_dir, exist_ok=True)
    pairs: List[Tuple[str, str]] = []
    only_species: List[str] = []
    with open(txt_path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            if "\t" in s:
                sp, ge = s.split("\t", 1)
                sp = sp.strip()
                ge = ge.strip()
                if sp and ge:
                    pairs.append((sp, ge))
            else:
                only_species.append(s)
    success_paths: List[str] = []
    for sp, ge in pairs:
        ok, p = fetch_one(sp, ge, out_dir=out_dir, db=db, delay=delay)
        if ok and p:
            success_paths.append(p)
    combined = os.path.join(out_dir, combine_filename)
    with open(combined, "w", encoding="utf-8") as out:
        for fp in success_paths:
            with open(fp, "r", encoding="utf-8") as f:
                out.write(f.read().rstrip() + "\n")
    print(f"📦 병합 저장: {combined}  (entries: {len(success_paths)})")
    if only_species:
        print(f"ℹ️ 종만 있는 라인 {len(only_species)}개 발견. 별도 genes 인자와 함께 fetch_many() 호출을 고려하세요.")
    return combined


# ---------------------------------------------------------------------------
# FetchPane GUI (NCBI 시퀀스 수집기)
# ---------------------------------------------------------------------------

class FetchPane(tk.Frame):
    def __init__(self, master, on_send_to_gblocks, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(bg="#f9f9f9")
        self.on_send_to_gblocks = on_send_to_gblocks

        tk.Label(self, text="📥 NCBI 시퀀스 수집기", font=("Helvetica", 13, "bold"), bg="#f9f9f9").pack(pady=6)

        row = tk.Frame(self, bg="#f9f9f9")
        row.pack(fill="x", padx=8, pady=2)
        tk.Label(row, text="Entrez.email", width=14, anchor="w", bg="#f9f9f9").pack(side=tk.LEFT)
        self.email = tk.StringVar(value="you@example.com")
        tk.Entry(row, textvariable=self.email).pack(side=tk.LEFT, fill="x", expand=True)

        row_db = tk.Frame(self, bg="#f9f9f9")
        row_db.pack(fill="x", padx=8, pady=2)
        tk.Label(row_db, text="데이터베이스", width=14, anchor="w", bg="#f9f9f9").pack(side=tk.LEFT)
        self.db = tk.StringVar(value="protein")
        tk.Radiobutton(row_db, text="단백질 (protein)", variable=self.db, value="protein", bg="#f9f9f9").pack(side=tk.LEFT, padx=4)
        tk.Radiobutton(row_db, text="뉴클레오타이드 (nucleotide)", variable=self.db, value="nucleotide", bg="#f9f9f9").pack(side=tk.LEFT, padx=4)

        row2 = tk.Frame(self, bg="#f9f9f9")
        row2.pack(fill="x", padx=8, pady=2)
        tk.Label(row2, text="출력 폴더", width=14, anchor="w", bg="#f9f9f9").pack(side=tk.LEFT)
        self.out_dir = tk.StringVar(value="protein_results")
        tk.Entry(row2, textvariable=self.out_dir).pack(side=tk.LEFT, fill="x", expand=True)
        tk.Button(row2, text="폴더 선택", command=self.browse_out).pack(side=tk.LEFT, padx=4)

        row3 = tk.Frame(self, bg="#f9f9f9")
        row3.pack(fill="x", padx=8, pady=2)
        tk.Label(row3, text="유전자(콤마)", width=14, anchor="w", bg="#f9f9f9").pack(side=tk.LEFT)
        self.genes = tk.StringVar(value="LEP,UCP1,COX1")
        tk.Entry(row3, textvariable=self.genes).pack(side=tk.LEFT, fill="x", expand=True)

        tk.Label(self, text="종 목록 (한 줄당 하나)", anchor="w", bg="#f9f9f9").pack(fill="x", padx=8)
        self.spec_box = tk.Text(self, height=12, width=40)
        self.spec_box.pack(fill="both", expand=False, padx=8, pady=4)

        up = tk.Frame(self, bg="#f9f9f9")
        up.pack(fill="x", padx=8, pady=2)
        tk.Button(up, text="TXT 불러오기 (species[\\tgene])", command=self.load_txt).pack(side=tk.LEFT)
        self.txt_path = tk.StringVar(value="")
        tk.Label(up, textvariable=self.txt_path, bg="#f9f9f9", fg="#555").pack(side=tk.LEFT, padx=6)

        btns = tk.Frame(self, bg="#f9f9f9")
        btns.pack(fill="x", padx=8, pady=6)
        tk.Button(btns, text="NCBI 수집 → 멀티FASTA 생성", command=self.run_fetch).pack(side=tk.LEFT)
        tk.Button(btns, text="멀티FASTA를 Gblocks로 보내기", command=self.send_to_gblocks).pack(side=tk.LEFT, padx=6)

        self.status = tk.Label(self, text="준비됨", bg="#f9f9f9", fg="#333")
        self.status.pack(padx=8, pady=4)
        self.log = tk.Text(self, height=8, width=40, state="disabled")
        self.log.pack(fill="both", expand=True, padx=8, pady=4)

        self.combined_path = None

    def browse_out(self):
        d = filedialog.askdirectory(title="출력 폴더 선택")
        if d:
            self.out_dir.set(d)

    def load_txt(self):
        p = filedialog.askopenfilename(title="TXT 선택", filetypes=[("Text", "*.txt"), ("All", "*.*")])
        if p:
            self.txt_path.set(p)
            self.status.config(text=f"TXT 로드됨: {os.path.basename(p)}")

    def _log(self, msg: str):
        self.log.configure(state="normal")
        self.log.insert("end", msg.rstrip() + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def run_fetch(self):
        email = self.email.get().strip()
        outd = self.out_dir.get().strip() or "protein_results"
        db = self.db.get()
        genes = [g.strip() for g in self.genes.get().split(",") if g.strip()]
        txt = self.txt_path.get().strip()
        species = [s.strip() for s in self.spec_box.get("1.0", "end").splitlines() if s.strip()]

        if not email or "@" not in email:
            messagebox.showwarning("경고", "유효한 Entrez.email을 입력하세요.")
            return
        if not genes and not txt:
            messagebox.showwarning("경고", "유전자 목록을 입력하거나 TXT를 불러오세요.")
            return
        self.status.config(text="⏳ 수집 중…(창을 닫지 마세요)")
        self._log(f"[시작] db={db} | out_dir={outd} | genes={','.join(genes)} | txt={os.path.basename(txt) if txt else '-'}")

        def worker():
            try:
                if txt:
                    self.combined_path = fetch_from_txt(txt, out_dir=outd, email=email, db=db)
                else:
                    self.combined_path = fetch_many(species, genes, out_dir=outd, email=email, db=db)
                self.status.after(0, lambda: self.status.config(text=f"✅ 완료: {self.combined_path}"))
                self._log(f"[완료] {self.combined_path}")
            except Exception as e:
                self.status.after(0, lambda: self.status.config(text="❌ 오류"))
                self._log(f"[오류] {e}")
                traceback = True

        threading.Thread(target=worker, daemon=True).start()

    def send_to_gblocks(self):
        if not self.combined_path or not os.path.isfile(self.combined_path):
            messagebox.showwarning("경고", "먼저 NCBI 수집을 실행해 멀티-FASTA를 생성하세요.")
            return
        try:
            self.on_send_to_gblocks(self.combined_path)
            self.status.config(text=f"➡️ 보냈습니다: {self.combined_path}")
            self._log(f"[전달] Gblocks에 입력 파일 설정 완료")
        except Exception as e:
            messagebox.showerror("에러", str(e))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="Entrez.email (필수)")
    ap.add_argument("--txt", help="입력 TXT (species[\\tgene])")
    ap.add_argument("--out_dir", default="protein_results")
    ap.add_argument("--genes", help="콤마로 구분된 유전자 목록 (TXT가 species-only일 때 사용)")
    ap.add_argument(
        "--db",
        default="protein",
        choices=["protein", "nucleotide"],
        help="데이터베이스 타입: protein 또는 nucleotide (기본값: protein)",
    )
    args = ap.parse_args()

    if args.txt:
        fetch_from_txt(args.txt, args.out_dir, args.email, db=args.db)
    else:
        print("TXT가 없으면 --genes 와 종 목록(코드 수정)이 필요합니다.")
        sys.exit(0)
