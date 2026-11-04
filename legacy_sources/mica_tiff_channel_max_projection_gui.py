# -*- coding: utf-8 -*-
"""
Mica interleaved multi-channel Z-stack TIFF -> per-channel Maximum Intensity Projection

想定ページ順:
  2ch: Z1-Ch1, Z1-Ch2, Z2-Ch1, Z2-Ch2, ...
  3ch: Z1-Ch1, Z1-Ch2, Z1-Ch3, Z2-Ch1, ...
  4ch: Z1-Ch1, Z1-Ch2, Z1-Ch3, Z1-Ch4, Z2-Ch1, ...

各チャンネルについて Z方向の Maximum Intensity Projection を作成し、
元ファイル名 + 指定postfix + .tif で保存します。

例:
  input: sample.tif
  Ch1 postfix: _DAPI_MAX
  Ch2 postfix: _ZO1_MAX

  output:
    sample_DAPI_MAX.tif
    sample_ZO1_MAX.tif

必要ライブラリ:
  pip install numpy tifffile
"""

import os
import threading
import traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import tifffile


def read_interleaved_tiff(path):
    with tifffile.TiffFile(path) as tif:
        if len(tif.pages) < 1:
            raise ValueError("TIFFにページがありません。")

        first = tif.pages[0].asarray()
        if first.ndim != 2:
            raise ValueError(
                f"各ページが2D画像ではありません (shape={first.shape})。"
                "このスクリプトはMicaのmultipage grayscale TIFFを想定しています。"
            )

        h, w = first.shape
        dtype = first.dtype
        n_pages = len(tif.pages)

        stack = np.empty((n_pages, h, w), dtype=dtype)
        stack[0] = first

        for i in range(1, n_pages):
            arr = tif.pages[i].asarray()
            if arr.shape != (h, w):
                raise ValueError(
                    f"ページサイズが不一致です: page 1={(h,w)}, page {i+1}={arr.shape}"
                )
            stack[i] = arr

        return stack


def make_max_projections(path, n_channels, postfixes, output_dir=None):
    path = os.path.abspath(path)
    stack = read_interleaved_tiff(path)
    n_pages = stack.shape[0]

    if n_pages % n_channels != 0:
        raise ValueError(
            f"総ページ数 {n_pages} がチャンネル数 {n_channels} で割り切れません。"
            "\nページ順またはチャンネル数を確認してください。"
        )

    n_z = n_pages // n_channels
    zcyx = stack.reshape(n_z, n_channels, stack.shape[1], stack.shape[2])

    stem = os.path.splitext(os.path.basename(path))[0]
    if output_dir is None:
        output_dir = os.path.dirname(path)
    os.makedirs(output_dir, exist_ok=True)

    results = []
    for c in range(n_channels):
        ch_stack = zcyx[:, c, :, :]
        max_img = np.max(ch_stack, axis=0)

        postfix = postfixes[c].strip() or f"_Ch{c+1}_MAX"
        out_name = f"{stem}{postfix}.tif"
        out_path = os.path.join(output_dir, out_name)

        tifffile.imwrite(out_path, max_img, photometric="minisblack", metadata=None)

        results.append({
            "channel": c + 1,
            "n_z": n_z,
            "shape": max_img.shape,
            "dtype": str(max_img.dtype),
            "output": out_path,
        })
    return results


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Mica TIFF Channel MAX Projection")
        self.geometry("860x650")
        self.minsize(760, 560)

        self.files = []
        self.output_dir = tk.StringVar(value="")
        self.n_channels = tk.IntVar(value=2)
        self.postfix_vars = [
            tk.StringVar(value="_Ch1_MAX"),
            tk.StringVar(value="_Ch2_MAX"),
            tk.StringVar(value="_Ch3_MAX"),
            tk.StringVar(value="_Ch4_MAX"),
        ]
        self._build_ui()
        self._update_postfix_state()

    def _build_ui(self):
        pad = 10

        top = ttk.Frame(self, padding=pad)
        top.pack(fill="x")
        ttk.Button(top, text="TIFFファイルを選択", command=self.select_files).pack(side="left", padx=(0,8))
        ttk.Button(top, text="フォルダ内TIFFを追加", command=self.select_folder).pack(side="left", padx=(0,8))
        ttk.Button(top, text="リストをクリア", command=self.clear_files).pack(side="left")

        list_frame = ttk.LabelFrame(self, text="処理対象ファイル", padding=pad)
        list_frame.pack(fill="both", expand=True, padx=pad, pady=(0,pad))
        self.listbox = tk.Listbox(list_frame, height=10)
        self.listbox.pack(side="left", fill="both", expand=True)
        sb = ttk.Scrollbar(list_frame, orient="vertical", command=self.listbox.yview)
        sb.pack(side="right", fill="y")
        self.listbox.configure(yscrollcommand=sb.set)

        settings = ttk.LabelFrame(self, text="設定", padding=pad)
        settings.pack(fill="x", padx=pad, pady=(0,pad))

        row1 = ttk.Frame(settings)
        row1.pack(fill="x", pady=4)
        ttk.Label(row1, text="チャンネル数:").pack(side="left")
        ch_spin = ttk.Spinbox(
            row1, from_=2, to=4, width=5,
            textvariable=self.n_channels,
            command=self._update_postfix_state
        )
        ch_spin.pack(side="left", padx=(8,20))
        ch_spin.bind("<KeyRelease>", lambda e: self._update_postfix_state())

        ttk.Label(row1, text="ページ順: Z1-Ch1, Z1-Ch2, ... → Z2-Ch1, Z2-Ch2, ...").pack(side="left")

        pf_frame = ttk.Frame(settings)
        pf_frame.pack(fill="x", pady=(8,4))
        self.postfix_entries = []
        for i in range(4):
            f = ttk.Frame(pf_frame)
            f.pack(side="left", padx=(0,12))
            ttk.Label(f, text=f"Ch{i+1} postfix").pack(anchor="w")
            ent = ttk.Entry(f, textvariable=self.postfix_vars[i], width=16)
            ent.pack()
            self.postfix_entries.append(ent)

        out_frame = ttk.Frame(settings)
        out_frame.pack(fill="x", pady=(8,4))
        ttk.Label(out_frame, text="出力先:").pack(side="left")
        ttk.Entry(out_frame, textvariable=self.output_dir).pack(side="left", fill="x", expand=True, padx=8)
        ttk.Button(out_frame, text="選択", command=self.select_output_dir).pack(side="left")
        ttk.Button(out_frame, text="元ファイルと同じ場所", command=lambda: self.output_dir.set("")).pack(side="left", padx=(8,0))

        ttk.Label(settings, text="出力先を空欄にすると、各入力TIFFと同じフォルダに保存します。").pack(anchor="w", pady=(4,0))

        run_frame = ttk.Frame(self, padding=(pad,0,pad,pad))
        run_frame.pack(fill="x")
        self.run_btn = ttk.Button(run_frame, text="MAX画像を作成", command=self.start_processing)
        self.run_btn.pack(side="left")
        self.progress = ttk.Progressbar(run_frame, mode="determinate")
        self.progress.pack(side="left", fill="x", expand=True, padx=(12,0))

        log_frame = ttk.LabelFrame(self, text="ログ", padding=pad)
        log_frame.pack(fill="both", expand=True, padx=pad, pady=(0,pad))
        self.log = tk.Text(log_frame, height=10, wrap="word")
        self.log.pack(fill="both", expand=True)

    def _update_postfix_state(self):
        try:
            n = int(self.n_channels.get())
        except Exception:
            n = 2
        n = max(2, min(4, n))
        for i, ent in enumerate(self.postfix_entries):
            ent.configure(state=("normal" if i < n else "disabled"))

    def select_files(self):
        paths = filedialog.askopenfilenames(
            title="Mica TIFFを選択",
            filetypes=[("TIFF files", "*.tif *.tiff"), ("All files", "*.*")],
        )
        self._add_files(paths)

    def select_folder(self):
        folder = filedialog.askdirectory(title="TIFFフォルダを選択")
        if not folder:
            return
        paths = [
            os.path.join(folder, name)
            for name in os.listdir(folder)
            if name.lower().endswith((".tif", ".tiff"))
            and os.path.isfile(os.path.join(folder, name))
        ]
        self._add_files(sorted(paths))

    def _add_files(self, paths):
        existing = set(self.files)
        for p in paths:
            p = os.path.abspath(p)
            if p not in existing:
                self.files.append(p)
                self.listbox.insert("end", p)
                existing.add(p)

    def clear_files(self):
        self.files.clear()
        self.listbox.delete(0, "end")
        self.log.delete("1.0", "end")
        self.progress["value"] = 0

    def select_output_dir(self):
        folder = filedialog.askdirectory(title="出力フォルダを選択")
        if folder:
            self.output_dir.set(folder)

    def write_log(self, msg):
        self.after(0, self._write_log_ui, msg)

    def _write_log_ui(self, msg):
        self.log.insert("end", msg + "\n")
        self.log.see("end")

    def start_processing(self):
        if not self.files:
            messagebox.showwarning("入力なし", "TIFFファイルを選択してください。")
            return
        try:
            n_channels = int(self.n_channels.get())
        except Exception:
            messagebox.showerror("設定エラー", "チャンネル数は2〜4を指定してください。")
            return
        if n_channels not in (2, 3, 4):
            messagebox.showerror("設定エラー", "チャンネル数は2〜4を指定してください。")
            return

        postfixes = [self.postfix_vars[i].get() for i in range(n_channels)]
        output_dir = self.output_dir.get().strip() or None

        self.run_btn.configure(state="disabled")
        self.progress.configure(maximum=len(self.files), value=0)
        self.log.delete("1.0", "end")

        threading.Thread(
            target=self._process_worker,
            args=(n_channels, postfixes, output_dir),
            daemon=True,
        ).start()

    def _process_worker(self, n_channels, postfixes, output_dir):
        ok = 0
        ng = 0
        for i, path in enumerate(self.files, start=1):
            try:
                self.write_log(f"[{i}/{len(self.files)}] {os.path.basename(path)}")
                results = make_max_projections(path, n_channels, postfixes, output_dir)
                for r in results:
                    self.write_log(f"  Ch{r['channel']}: Z={r['n_z']} → {os.path.basename(r['output'])}")
                ok += 1
            except Exception as e:
                ng += 1
                self.write_log(f"  ERROR: {e}")
                self.write_log(traceback.format_exc())
            self.after(0, lambda v=i: self.progress.configure(value=v))

        self.write_log("")
        self.write_log(f"完了: 成功 {ok} ファイル / エラー {ng} ファイル")
        self.after(0, lambda: self.run_btn.configure(state="normal"))

        if ng == 0:
            self.after(0, lambda: messagebox.showinfo("完了", f"{ok}ファイルの処理が完了しました。"))
        else:
            self.after(0, lambda: messagebox.showwarning("完了", f"成功 {ok}、エラー {ng}。\n詳細はログを確認してください。"))


if __name__ == "__main__":
    App().mainloop()
