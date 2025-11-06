import tkinter as tk
from tkinter import ttk
import os
from app import GblockApp
from fetch_pane import FetchPane

def main():
    root = tk.Tk()
    root.withdraw()

    windows = {"fetch": None}

    gb_win = tk.Toplevel(root)
    gb_win.title("Gblocks")
    gb_win.geometry("980x820")

    def open_fetch_window():
        if windows["fetch"] is not None and windows["fetch"].winfo_exists():
            try:
                windows["fetch"].deiconify()
                windows["fetch"].lift()
                windows["fetch"].focus_force()
            except Exception:
                pass
            return

        fwin = tk.Toplevel(root)
        fwin.title("NCBI Fetch (Protein)")
        fwin.geometry("900x700")

        def on_send_to_gblocks(path: str):
            gblocks_ui.file_path.set(path)
            gblocks_ui.status.config(text=f"입력 파일 설정됨: {path}")
            try:
                gb_win.deiconify()
                gb_win.lift()
                gb_win.focus_force()
            except Exception:
                pass

        pane_root = tk.Frame(fwin, bg="#f9f9f9")
        pane_root.pack(fill="both", expand=True)
        fetch_ui = FetchPane(pane_root, on_send_to_gblocks=on_send_to_gblocks)
        fetch_ui.pack(fill="both", expand=True)

        fwin.protocol("WM_DELETE_WINDOW", fwin.destroy)
        windows["fetch"] = fwin

    gblocks_ui = GblockApp(gb_win, on_open_fetch=open_fetch_window)

    gb_win.protocol("WM_DELETE_WINDOW", root.destroy)

    root.mainloop()

if __name__ == "__main__":
    main()
