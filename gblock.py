import argparse, os, shutil, subprocess, re
from typing import List, Tuple, Dict, Optional
from tempfile import TemporaryDirectory

def parse_fasta(path: str) -> List[Tuple[str, str]]:
    xs, h, buf = [], None, []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s: continue
            if s.startswith(">"):
                if h is not None: xs.append((h, "".join(buf)))
                h, buf = s[1:].strip() or "unnamed", []
            else:
                buf.append(s)
    if h is not None: xs.append((h, "".join(buf)))
    if not xs: raise ValueError("No FASTA entries found.")
    return xs

def write_fasta(entries: List[Tuple[str, str]], path: str, wrap: int = 70) -> None:
    with open(path, "w", encoding="utf-8") as out:
        for h, seq in entries:
            out.write(f">{h}\n")
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    out.write(seq[i:i+wrap] + "\n")
            else:
                out.write(seq + "\n")

def all_equal_length(entries: List[Tuple[str, str]]) -> bool:
    return len({len(s) for _, s in entries}) == 1

def which_or_path(x: Optional[str]) -> Optional[str]:
    if not x: return None
    if os.path.isabs(x) and os.path.isfile(x): return x
    return shutil.which(x)

def run_mafft(in_fa: str, out_fa: str, exe: Optional[str] = None):
    exe = which_or_path(exe) or shutil.which("mafft")
    if not exe: raise RuntimeError("MAFFT not found.")
    cmd = [exe, "--auto", "--quiet", in_fa]
    with open(out_fa, "w", encoding="utf-8") as fout:
        p = subprocess.run(cmd, text=True, stdout=fout, stderr=subprocess.PIPE)
    if p.returncode != 0: raise RuntimeError(f"MAFFT failed: {p.stderr.strip()}")

def run_muscle(in_fa: str, out_fa: str, exe: Optional[str] = None):
    exe = which_or_path(exe) or shutil.which("muscle")
    if not exe: raise RuntimeError("MUSCLE not found.")
    # MUSCLE v5
    cmd = [exe, "-align", in_fa, "-output", out_fa]
    p = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0: raise RuntimeError(f"MUSCLE failed: {p.stderr.strip()}")

def nw_align(a: str, b: str, match: int = 1, mismatch: int = -1, gap: int = -2):
    n, m = len(a), len(b)
    score = [[0]*(m+1) for _ in range(n+1)]
    trace = [[0]*(m+1) for _ in range(n+1)]
    for i in range(1, n+1):
        score[i][0] = i*gap; trace[i][0] = 2
    for j in range(1, m+1):
        score[0][j] = j*gap; trace[0][j] = 3
    for i in range(1, n+1):
        ai = a[i-1]
        for j in range(1, m+1):
            bj = b[j-1]
            sdiag = score[i-1][j-1] + (match if ai==bj else mismatch)
            sup   = score[i-1][j]   + gap
            sleft = score[i][j-1]   + gap
            if sdiag >= sup and sdiag >= sleft:
                score[i][j]=sdiag; trace[i][j]=1
            elif sup >= sleft:
                score[i][j]=sup;   trace[i][j]=2
            else:
                score[i][j]=sleft; trace[i][j]=3
    i,j=n,m; A,B=[],[]
    while i>0 or j>0:
        t=trace[i][j]
        if t==1: A.append(a[i-1]); B.append(b[j-1]); i-=1; j-=1
        elif t==2: A.append(a[i-1]); B.append('-');   i-=1
        else:      A.append('-');     B.append(b[j-1]); j-=1
    return "".join(reversed(A)), "".join(reversed(B))

def gap_pattern_from_aligned_base(base_ref: str, aln_base: str):
    L=len(base_ref); gaps_before=[0]*(L+1); i_base=0; run=0
    for ch in aln_base:
        if ch=='-': run+=1
        else:
            if i_base<L: gaps_before[i_base]+=run
            run=0; i_base+=1
    gaps_before[L]+=run
    return gaps_before

def merge_gap_patterns_max(p1, p2):
    L=max(len(p1),len(p2)); out=[]
    for i in range(L):
        a=p1[i] if i<len(p1) else 0
        b=p2[i] if i<len(p2) else 0
        out.append(max(a,b))
    return out

def expand_to_pattern(base_ref: str, aln_seq: str, target_gaps):
    L=len(base_ref)
    aln_base=[]; j_base=0
    for ch in aln_seq:
        if ch=='-': aln_base.append('-')
        else: aln_base.append(base_ref[j_base]); j_base+=1
    out=[]; i_aln=0
    for k in range(L):
        out.extend('-'*target_gaps[k])
        while i_aln<len(aln_seq):
            ch_a=aln_seq[i_aln]; ch_b=aln_base[i_aln]; i_aln+=1
            out.append(ch_a)
            if ch_b!='-': break
    out.extend('-'*target_gaps[L])
    return "".join(out)

def center_star_align(entries: List[Tuple[str,str]], match=1, mismatch=-1, gap=-2)->List[Tuple[str,str]]:
    center_idx = max(range(len(entries)), key=lambda i: len(entries[i][1]))
    base_ref = entries[center_idx][1]; L=len(base_ref)
    master_gaps=[0]*(L+1)
    aligned_headers=[]; aligned_seqs=[]
    for h,seq in entries:
        aln_base, aln_seq = nw_align(base_ref, seq, match, mismatch, gap)
        gaps_new = gap_pattern_from_aligned_base(base_ref, aln_base)
        merged = merge_gap_patterns_max(master_gaps, gaps_new)
        aligned_seqs = [expand_to_pattern(base_ref, s_prev, merged) for s_prev in aligned_seqs]
        s_new = expand_to_pattern(base_ref, aln_seq, merged)
        master_gaps = merged
        aligned_headers.append(h); aligned_seqs.append(s_new)
    order={h:i for i,(h,_) in enumerate(entries)}
    out=list(zip(aligned_headers, aligned_seqs))
    out.sort(key=lambda x: order[x[0]])
    return out

def auto_align_if_needed(entries: List[Tuple[str,str]], aligner="auto",
                         aligner_exe=None, fallback_align=True)->List[Tuple[str,str]]:
    if len(entries)<=1 or all_equal_length(entries): return entries
    chosen = aligner.lower()
    if chosen=="auto":
        chosen = "mafft" if shutil.which("mafft") else ("muscle" if shutil.which("muscle") else "none")
    if chosen=="mafft":
        with TemporaryDirectory() as w:
            inp=os.path.join(w,"in.fasta"); outp=os.path.join(w,"out.fasta")
            write_fasta(entries, inp, wrap=0); run_mafft(inp,outp,exe=aligner_exe)
            return parse_fasta(outp)
    if chosen=="muscle":
        with TemporaryDirectory() as w:
            inp=os.path.join(w,"in.fasta"); outp=os.path.join(w,"out.fasta")
            write_fasta(entries, inp, wrap=0); run_muscle(inp,outp,exe=aligner_exe)
            return parse_fasta(outp)
    if fallback_align:
        return center_star_align(entries)
    raise RuntimeError("Sequences are unaligned and no aligner found. Install MAFFT or MUSCLE, or enable fallback_align.")

def column_metrics(aligned_seqs: List[str]) -> List[Dict[str,float]]:
    n=len(aligned_seqs); L=len(aligned_seqs[0]); mets=[]
    for i in range(L):
        col=[s[i] for s in aligned_seqs]
        gaps=sum(1 for c in col if c=='-')
        non=[c for c in col if c!='-']
        gap_ratio=gaps/n
        if not non: cons_ratio=0.0
        else:
            cnt={}
            for c in non: cnt[c]=cnt.get(c,0)+1
            cons_ratio=max(cnt.values())/len(non)
        mets.append({"gap_ratio":gap_ratio,"cons_ratio":cons_ratio})
    return mets

def _allow_short_nonconserved_runs(mask: List[bool], max_run: int)->List[bool]:
    """
    Gblocks의 'Maximum # of Contiguous Nonconserved' 근사:
      - False가 연속 max_run 이하로 나오고, 양옆이 True인 경우 그 구간을 True로 승격(블록 내부 흡수)
      - 블록 경계 확장/연결 효과
    """
    if max_run<=0: return mask[:]
    m=mask[:]; L=len(m); i=0
    while i<L:
        if not m[i]:
            j=i
            while j<L and not m[j]: j+=1
            run=j-i  # [i, j)
            left_true = (i-1>=0 and m[i-1])
            right_true= (j<L and m[j])
            if run<=max_run and left_true and right_true:
                for k in range(i,j): m[k]=True
            i=j
        else:
            i+=1
    return m

def find_blocks(mask: List[bool], min_len:int)->List[Tuple[int,int]]:
    blocks=[]; s=None
    for i,keep in enumerate(mask+[False]):
        if keep and s is None: s=i
        elif (not keep) and s is not None:
            if i-s>=min_len: blocks.append((s,i))
            s=None
    return blocks

def soft_trim_block(b, mets, flank_cons, max_gap):
    s,e=b
    while s<e and not (mets[s]["gap_ratio"]<=max_gap and mets[s]["cons_ratio"]>=flank_cons): s+=1
    while e>s and not (mets[e-1]["gap_ratio"]<=max_gap and mets[e-1]["cons_ratio"]>=flank_cons): e-=1
    return (s,e) if e>s else None

def gblock_filter(entries: List[Tuple[str,str]],
                  min_block_len=10, max_gap=0.5, min_cons=0.7, flank_cons:Optional[float]=None,
                  drop_all_gap_columns=False,
                  max_noncons_run:int=0
                  )->List[Tuple[str,str]]:
    if not all_equal_length(entries):
        raise ValueError("All sequences must be equal length (pre-aligned MSA required).")
    headers=[h for h,_ in entries]; seqs=[s for _,s in entries]
    mets=column_metrics(seqs)

    mask=[(m["gap_ratio"]<=max_gap and m["cons_ratio"]>=min_cons) for m in mets]
    if max_noncons_run>0:
        mask = _allow_short_nonconserved_runs(mask, max_noncons_run)

    blocks=find_blocks(mask, min_block_len)

    if flank_cons is not None:
        blocks=[tb for b in blocks if (tb:=soft_trim_block(b,mets,flank_cons,max_gap)) and (tb[1]-tb[0])>=min_block_len]

    if not blocks: return [(h,"") for h in headers]
    keep=[]
    for s,e in blocks: keep.extend(range(s,e))
    if drop_all_gap_columns:
        keep=[i for i in keep if any(seq[i] != '-' for seq in seqs)]
    out=[]
    for h,seq in entries:
        out.append((h, "".join(seq[i] for i in keep)))
    return out

def gblock_pipeline(entries: List[Tuple[str,str]],
                    auto_align=True, aligner="auto", aligner_exe=None, fallback_align=True,
                    min_block_len=10, max_gap=0.5, min_cons=0.7, flank_cons:Optional[float]=None,
                    drop_all_gap_columns=True,
                    max_noncons_run:int=0
                    )->List[Tuple[str,str]]:
    aln = auto_align_if_needed(entries, aligner=aligner, aligner_exe=aligner_exe, fallback_align=fallback_align) if auto_align else entries
    if not all_equal_length(aln): raise ValueError("Alignment failed; sequences still unequal length.")
    return gblock_filter(aln, min_block_len, max_gap, min_cons, flank_cons, drop_all_gap_columns, max_noncons_run)

def gb_allowed_gap_to_max_gap(mode:str)->float:
    mode=mode.lower()
    if mode=="none": return 0.0
    if mode=="half": return 0.5
    if mode=="all":  return 1.0
    try:
        v=float(mode); return max(0.0, min(1.0, v))
    except:
        return 0.5

def convert_gblocks_params_to_internal(N:int,
                                       gb_min_cons_count:int=9,
                                       gb_flank_cons_count:int=14,
                                       gb_max_noncons_run:int=8,
                                       gb_min_block_len:int=10,
                                       gb_allowed_gap:str="None"):
    N=max(1,N)
    mc=max(1, min(gb_min_cons_count, N))
    mf=max(1, min(gb_flank_cons_count, N))
    min_cons = mc / N
    flank_cons = mf / N
    max_gap = gb_allowed_gap_to_max_gap(gb_allowed_gap)
    min_block_len = max(1, gb_min_block_len)
    max_noncons_run = max(0, gb_max_noncons_run)
    return min_block_len, max_gap, min_cons, flank_cons, max_noncons_run

def main():
    p=argparse.ArgumentParser(description="Gblocks-like trimmer with auto/fallback alignment.")
    p.add_argument("--in", dest="in_fasta", required=True, help="Input FASTA/TXT (unaligned or aligned).")
    p.add_argument("--out", dest="out_fasta", required=True, help="Output trimmed FASTA.")

    p.add_argument("--auto_align", type=lambda x:str(x).lower() not in {"0","false","no"}, default=True)
    p.add_argument("--aligner", choices=["auto","mafft","muscle","none"], default="auto")
    p.add_argument("--aligner_exe", default=None)
    p.add_argument("--fallback_align", type=lambda x:str(x).lower() not in {"0","false","no"}, default=True,
                   help="Use built-in center-star aligner if no MAFFT/MUSCLE found (default: True).")

    p.add_argument("--min_block_len", type=int, default=10)
    p.add_argument("--max_gap", type=float, default=0.5)
    p.add_argument("--min_cons", type=float, default=0.7)
    p.add_argument("--flank_cons", type=float, default=None)
    p.add_argument("--drop_all_gap_columns", action="store_true")
    p.add_argument("--max_noncons_run", type=int, default=0,
                   help="Maximum # of contiguous nonconserved columns to absorb within blocks (approx. Gblocks).")

    p.add_argument("--use_gblocks_params", type=lambda x:str(x).lower() not in {"0","false","no"}, default=False,
                   help="Use Gblocks-style integer parameters below (auto-converted using N).")
    p.add_argument("--gb_min_cons_count",   type=int, default=9,  help="Minimum # of sequences for a conserved position.")
    p.add_argument("--gb_flank_cons_count", type=int, default=14, help="Minimum # for a flanking position.")
    p.add_argument("--gb_max_noncons_run",  type=int, default=8,  help="Maximum # of contiguous nonconserved positions.")
    p.add_argument("--gb_min_block_len",    type=int, default=10, help="Minimum length of a block.")
    p.add_argument("--gb_allowed_gap",      type=str, default="None", choices=["None","Half","All"],
                   help="Allowed Gap Positions: None/Half/All. (maps to max_gap 0/0.5/1.0)")
    p.add_argument("--gb_use_similarity",   type=lambda x:str(x).lower() not in {"0","false","no"}, default=True,
                   help="Display only; engine uses identity (no matrices).")

    args=p.parse_args()

    entries=parse_fasta(args.in_fasta)

    if args.use_gblocks_params:
        N=len(entries)
        mb, mg, mc, fc, mnr = convert_gblocks_params_to_internal(
            N,
            args.gb_min_cons_count,
            args.gb_flank_cons_count,
            args.gb_max_noncons_run,
            args.gb_min_block_len,
            args.gb_allowed_gap
        )
        min_block_len = mb
        max_gap       = mg
        min_cons      = mc
        flank_cons    = fc
        max_noncons_run = mnr
    else:
        min_block_len = args.min_block_len
        max_gap       = args.max_gap
        min_cons      = args.min_cons
        flank_cons    = args.flank_cons
        max_noncons_run = args.max_noncons_run

    trimmed=gblock_pipeline(
        entries,
        auto_align=args.auto_align,
        aligner=args.aligner,
        aligner_exe=args.aligner_exe,
        fallback_align=args.fallback_align,
        min_block_len=min_block_len,
        max_gap=max_gap,
        min_cons=min_cons,
        flank_cons=flank_cons,
        drop_all_gap_columns=args.drop_all_gap_columns,
        max_noncons_run=max_noncons_run,
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.out_fasta)) or ".", exist_ok=True)
    write_fasta(trimmed, args.out_fasta, wrap=70)

    in_len = len(entries[0][1])
    out_len = len(trimmed[0][1]) if trimmed else 0
    print("=== Gblock pipeline complete ===")
    print(f"- Input sequences: {len(entries)}")
    print(f"- First seq length (pre):  {in_len}")
    print(f"- First seq length (post): {out_len}")
    print(f"- Params: use_gblocks={args.use_gblocks_params} | min_block_len={min_block_len} "
          f"| max_gap={max_gap} | min_cons={min_cons} | flank_cons={flank_cons} | max_noncons_run={max_noncons_run}")
    if args.use_gblocks_params:
        print(f"- Gblocks-style: mc={args.gb_min_cons_count}, mf={args.gb_flank_cons_count}, "
              f"mnr={args.gb_max_noncons_run}, mb={args.gb_min_block_len}, gap={args.gb_allowed_gap}, "
              f"similarity={args.gb_use_similarity}")

if __name__=="__main__":
    main()
