"""
protein_fetch_safe.py
---------------------
종 + 유전자명으로 NCBI에서 단백질 FASTA를 자동 수집하여 하나의 멀티-FASTA(.txt)로 저장.

기능 요약
- NCBI 검색(db="protein") → 첫 번째 결과를 efetch로 받아 FASTA 저장
- 오류/간헐적 실패 대비 재시도
- 요청 간격(기본 0.5s) 제한
- 입력: (1) 파이썬 인자 species, genes / (2) 텍스트 파일(한 줄당 'species<TAB>gene' 또는 'species'만)
- 출력: 출력 폴더에 개별 FASTA + 병합된 multi_fasta.txt 생성
- Biopython 필요: pip install biopython

사용 예시
----------
from protein_fetch_safe import fetch_from_txt, fetch_many

# 1) TXT로 입력
fetch_from_txt("pairs.txt", out_dir="protein_results", email="you@example.com")

# 2) 코드로 직접 입력
species = ["Ursus maritimus", "Panthera tigris"]
genes   = ["LEP", "UCP1", "COX1"]
fetch_many(species, genes, out_dir="protein_results", email="you@example.com")
"""

from __future__ import annotations
from typing import List, Tuple, Optional
import os, time, sys, traceback
from Bio import Entrez, SeqIO
from collections import defaultdict

DEFAULT_DELAY = 0.5  # 요청 간격(초)

def _ensure_email(email: Optional[str]):
    if not email or "@" not in email:
        raise ValueError("Entrez.email이 설정되지 않았습니다. 유효한 이메일을 입력하세요.")
    Entrez.email = email

def entrez_search_id(query: str, db: str="protein", retmax: int=10, retries: int=3, delay: float=2.0) -> List[str]:
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

def entrez_fetch_protein_fasta(ncbi_id: str, retries: int=3, delay: float=2.0):
    for attempt in range(retries):
        try:
            with Entrez.efetch(db="protein", id=ncbi_id, rettype="fasta", retmode="text") as handle:
                seq_record = SeqIO.read(handle, "fasta")
            return seq_record
        except Exception as e:
            print(f"⚠️ [다운로드 실패 {attempt+1}/{retries}] ID={ncbi_id} → {e}")
            time.sleep(delay)
    print(f"❌ 다운로드 완전 실패: {ncbi_id}")
    return None

def _sanitize_filename(s: str) -> str:
    return "".join(ch for ch in s.replace(" ", "_") if ch.isalnum() or ch in "._-")

def fetch_one(species: str, gene: str, out_dir: str, delay: float=DEFAULT_DELAY) -> Tuple[bool, Optional[str]]:
    os.makedirs(out_dir, exist_ok=True)
    query = f"{gene}[gene] AND {species}[organism]"
    ids = entrez_search_id(query)
    if not ids:
        print(f"🚫 {species} ({gene}) 검색 결과 없음")
        return False, None
    rec = entrez_fetch_protein_fasta(ids[0])
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
    seen = set(); out = []
    for x in seq:
        if x not in seen:
            seen.add(x); out.append(x)
    return out

def fetch_many(species_list: List[str], genes: List[str], out_dir: str, email: str,
               combine_filename: str="multi_fasta.txt", delay: float=DEFAULT_DELAY) -> str:
    """
    species_list의 각 종에 대해 genes 목록의 첫 번째 히트 FASTA를 받아 개별 파일로 저장하고,
    성공한 것들을 하나의 멀티-FASTA 텍스트(combine_filename)로 병합하여 반환.
    """
    _ensure_email(email)
    os.makedirs(out_dir, exist_ok=True)
    success_paths: List[str] = []
    species_list = [s.strip() for s in species_list if s and s.strip()]
    genes = _unique([g.strip() for g in genes if g and g.strip()])
    for sp in species_list:
        for g in genes:
            ok, p = fetch_one(sp, g, out_dir=out_dir, delay=delay)
            if ok and p:
                success_paths.append(p)

    grouped = defaultdict(list)
    for fp in success_paths:
        # 파일명에서 유전자 추출 (…_CO1.fasta)
        gene = os.path.splitext(os.path.basename(fp))[0].split("_")[-1]
        grouped[gene].append(fp)

    for gene, paths in grouped.items():
        if len(paths) < 2:
            print(f"⚠️ {gene}: 서열이 {len(paths)}개라 Gblocks용 병합을 건너뜁니다.")
            continue
        combined = os.path.join(out_dir, f"multi_fasta_{gene}.txt")
        with open(combined, "w", encoding="utf-8") as out:
            for fp in paths:
                base = os.path.basename(fp)
                prefix = os.path.splitext(base)[0]  # Ursus_maritimus_CO1
                with open(fp, "r", encoding="utf-8") as f:
                    lines = f.readlines()
                # 헤더에 '종_유전자 ' 접두어 추가
                for i, line in enumerate(lines):
                    if line.startswith(">"):
                        lines[i] = f">{prefix} " + line[1:]
                        break
                out.write("".join(lines).rstrip() + "\n")
        print(f"📦 병합 저장: {combined} (entries: {len(paths)})")
    return combined

def fetch_from_txt(txt_path: str, out_dir: str, email: str,
                   combine_filename: str="multi_fasta.txt", delay: float=DEFAULT_DELAY) -> str:
    """
    txt 포맷:
      - 'species<TAB>gene' 형식 또는
      - 'species' 단독 라인들 + 함수 인자로 genes 지정 시 사용

    반환: 병합된 multi_fasta.txt 경로
    """
    _ensure_email(email)
    os.makedirs(out_dir, exist_ok=True)

    pairs: List[Tuple[str,str]] = []
    only_species: List[str] = []

    with open(txt_path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            if "\t" in s:
                sp, ge = s.split("\t", 1)
                sp = sp.strip(); ge = ge.strip()
                if sp and ge:
                    pairs.append((sp, ge))
            else:
                only_species.append(s)

    success_paths: List[str] = []
    for sp, ge in pairs:
        ok, p = fetch_one(sp, ge, out_dir=out_dir, delay=delay)
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

if __name__ == "__main__":
    # 간단한 CLI
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="Entrez.email (필수)")
    ap.add_argument("--txt", help="입력 TXT (species[\\tgene])")
    ap.add_argument("--out_dir", default="protein_results")
    ap.add_argument("--genes", help="콤마로 구분된 유전자 목록 (TXT가 species-only일 때 사용)")
    args = ap.parse_args()

    if args.txt:
        combined = fetch_from_txt(args.txt, args.out_dir, args.email)
    else:
        print("TXT가 없으면 --genes 와 종 목록(코드 수정)이 필요합니다.")
        sys.exit(0)
