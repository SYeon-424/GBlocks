import tkinter as tk
from tkinter import filedialog, messagebox
import threading, os

from protein_fetch_safe import fetch_many, fetch_from_txt

class FetchPane(tk.Frame):
    def __init__(self, master, on_send_to_gblocks, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(bg="#f9f9f9")
        self.on_send_to_gblocks = on_send_to_gblocks  # callback(path: str)

        # 입력 위젯
        tk.Label(self, text="📥 NCBI 단백질 수집기", font=("Helvetica", 13, "bold"), bg="#f9f9f9").pack(pady=6)

        row = tk.Frame(self, bg="#f9f9f9"); row.pack(fill="x", padx=8, pady=2)
        tk.Label(row, text="Entrez.email", width=14, anchor="w", bg="#f9f9f9").pack(side=tk.LEFT)
        self.email = tk.StringVar(value="you@example.com")
        tk.Entry(row, textvariable=self.email).pack(side=tk.LEFT, fill="x", expand=True)

        row2 = tk.Frame(self, bg="#f9f9f9"); row2.pack(fill="x", padx=8, pady=2)
        tk.Label(row2, text="출력 폴더", width=14, anchor="w", bg="#f9f9f9").pack(side=tk.LEFT)
        self.out_dir = tk.StringVar(value="protein_results")
        tk.Entry(row2, textvariable=self.out_dir).pack(side=tk.LEFT, fill="x", expand=True)
        tk.Button(row2, text="폴더 선택", command=self.browse_out).pack(side=tk.LEFT, padx=4)

        row3 = tk.Frame(self, bg="#f9f9f9"); row3.pack(fill="x", padx=8, pady=2)
        tk.Label(row3, text="유전자(콤마)", width=14, anchor="w", bg="#f9f9f9").pack(side=tk.LEFT)
        self.genes = tk.StringVar(value="LEP,UCP1,COX1")
        tk.Entry(row3, textvariable=self.genes).pack(side=tk.LEFT, fill="x", expand=True)

        tk.Label(self, text="종 목록 (한 줄당 하나)", anchor="w", bg="#f9f9f9").pack(fill="x", padx=8)
        self.spec_box = tk.Text(self, height=12, width=40)
        self.spec_box.pack(fill="both", expand=False, padx=8, pady=4)

        # TXT 업로드
        up = tk.Frame(self, bg="#f9f9f9"); up.pack(fill="x", padx=8, pady=2)
        tk.Button(up, text="TXT 불러오기 (species[\\tgene])", command=self.load_txt).pack(side=tk.LEFT)
        self.txt_path = tk.StringVar(value="")
        tk.Label(up, textvariable=self.txt_path, bg="#f9f9f9", fg="#555").pack(side=tk.LEFT, padx=6)

        # 동작 버튼
        btns = tk.Frame(self, bg="#f9f9f9"); btns.pack(fill="x", padx=8, pady=6)
        tk.Button(btns, text="NCBI 수집 → 멀티FASTA 생성", command=self.run_fetch).pack(side=tk.LEFT)
        tk.Button(btns, text="멀티FASTA를 Gblocks로 보내기", command=self.send_to_gblocks).pack(side=tk.LEFT, padx=6)

        # 상태/로그
        self.status = tk.Label(self, text="준비됨", bg="#f9f9f9", fg="#333")
        self.status.pack(padx=8, pady=4)
        self.log = tk.Text(self, height=8, width=40, state="disabled")
        self.log.pack(fill="both", expand=True, padx=8, pady=4)

        self.combined_path = None

    def browse_out(self):
        d = filedialog.askdirectory(title="출력 폴더 선택")
        if d: self.out_dir.set(d)

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
        outd  = self.out_dir.get().strip() or "protein_results"
        genes = [g.strip() for g in self.genes.get().split(",") if g.strip()]
        txt   = self.txt_path.get().strip()
        species = [s.strip() for s in self.spec_box.get("1.0", "end").splitlines() if s.strip()]

        if not email or "@" not in email:
            messagebox.showwarning("경고", "유효한 Entrez.email을 입력하세요.")
            return
        if not genes and not txt:
            messagebox.showwarning("경고", "유전자 목록을 입력하거나 TXT를 불러오세요.")
            return
        self.status.config(text="⏳ 수집 중…(창을 닫지 마세요)")
        self._log(f"[시작] out_dir={outd} | genes={','.join(genes)} | txt={os.path.basename(txt) if txt else '-'}")

        def worker():
            try:
                if txt:
                    self.combined_path = fetch_from_txt(txt, out_dir=outd, email=email)
                else:
                    self.combined_path = fetch_many(species, genes, out_dir=outd, email=email)
                self.status.after(0, lambda: self.status.config(text=f"✅ 완료: {self.combined_path}"))
                self._log(f"[완료] {self.combined_path}")
            except Exception as e:
                self.status.after(0, lambda: self.status.config(text="❌ 오류"))
                self._log(f"[오류] {e}"); traceback = True

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
