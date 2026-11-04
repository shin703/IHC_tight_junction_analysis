# -*- coding: utf-8 -*-
"""
Caco-2 DAPI StarDist nuclei review GUI - viewport/cached annotation version
========================================================================

親フォルダを選択し、直下のTIFFをsample名でペアリングして一覧表示します。

例:
    con1_DAPI.tiff : DAPI
    con1_ZO1.tiff  : ZO-1
    con2_DAPI.tiff : DAPI
    con2_ZO1.tiff  : ZO-1

出力は親フォルダ/ZO1_junction_analysis/<sample名>/ に保存します。
GUI上でGroupを指定でき、全解析時にsample summaryとgroup summaryを出力します。

DAPI核は StarDist 2D_versatile_fluo で初期検出し、
貼付コードと同じ方式で手動修正します。

高速化の要点
------------
- 全画像を毎回resizeしない
- canvas上で見えているviewportだけを切り出して描画
- brightness変換済みbase画像をcache
- 核ごとのbbox / contour / centroidをcache
- 編集された核だけcontour cacheを更新
- ZO-1はボタンを押した時だけlazy-loadし、viewportだけoverlay
- 代表画像用ZO-1 QC GUIで前処理を1段階ずつ確認・調整可能
- 全未確認DAPIを一括StarDistし、未確認draftとして保存可能
- StarDist/ZO-1解析設定はユーザー共通設定として全フォルダ間で保持
- 一括StarDistでは接触核を8近傍で自動Merge
- junction intensityは各laneのBG-subtracted ZO-1上位N pixel平均を全lane/pairで平均
- 描画中は画像を毎mouse event再生成せず、mouse release後に更新

操作
----
左ドラッグ:
    未選択 -> 新しい核
    1核選択 -> その核へ追加

Ctrl + 左ドラッグ:
    消しゴム

ダブルクリック:
    核選択 / 解除

右ドラッグ:
    pan

Mouse wheel:
    zoom

Delete:
    選択核削除

Ctrl+Z:
    Undo

[ / ]:
    brush size

S:
    保存

N:
    保存して次へ

出力（各子フォルダ）
--------------------
ZO1_junction_analysis/
    dapi_nuclei_annotation.npz
    dapi_nuclei_labels.tif
    dapi_nuclei_binary.tif
    dapi_nuclei_review_overlay.png

必要:
pip install numpy==1.26.4 opencv-python pillow tifffile scipy
pip install tensorflow stardist csbdeep
# LZW TIFFをtifffileで読む場合:
pip install --no-deps imagecodecs==2024.12.30
"""

import csv
import json
import math
import os
import queue
import re
import sys
import threading
import traceback
from pathlib import Path


# ============================================================
# Safe stdout/stderr for windowed PyInstaller builds
# ============================================================
# PyInstaller with console=False sets sys.stdout / sys.stderr to None on
# Windows. StarDist internally uses tqdm for tiled prediction, and tqdm
# expects a writable stream. Without this guard, StarDist can fail with:
#     AttributeError: 'NoneType' object has no attribute 'write'
# Keep the devnull handles alive for the lifetime of the application.
_DEVNULL_STREAMS = []


def ensure_writable_standard_streams():
    global _DEVNULL_STREAMS

    for name in ("stdout", "stderr"):
        stream = getattr(sys, name, None)

        if stream is not None and hasattr(stream, "write"):
            continue

        try:
            replacement = open(
                os.devnull,
                "w",
                encoding="utf-8",
                buffering=1,
            )
        except Exception:
            # Last-resort in-memory writer. This path should almost never be used.
            import io
            replacement = io.StringIO()

        setattr(sys, name, replacement)
        _DEVNULL_STREAMS.append(replacement)


ensure_writable_standard_streams()

import cv2
import numpy as np
import tifffile

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from PIL import Image, ImageTk


# ============================================================
# Japanese / English UI translation
# ============================================================

_APP_LANGUAGE = "ja"
_LANGUAGE_PREF_PATH = Path.home() / "mica_zo1_gui_language.json"

# Japanese/English pairs used across the integrated GUI, editor windows,
# message boxes and file dialogs. Existing scientific terms such as DAPI,
# ZO-1, StarDist, MAX and ROI are intentionally preserved.
_I18N_PAIRS = [
    ("親フォルダを選択", "Select Parent Folder"),
    ("親フォルダ未選択", "No parent folder selected"),
    ("入力ファイル postfix / Channel postfix", "Input file postfix / Channel postfix"),
    ("再スキャン / Rescan", "Rescan"),
    ("MAX設定から取得", "Copy from MAX settings"),
    ("選択サンプルを開く", "Open Selected Sample"),
    ("最初の未保存サンプルを開く", "Open First Unreviewed Sample"),
    ("一覧を再読み込み", "Reload List"),
    ("全未確認を一括StarDist", "Run StarDist on All Unreviewed"),
    ("選択サンプルにGroupを設定", "Set Group for Selected Samples"),
    ("選択サンプルのGroupをクリア", "Clear Group for Selected Samples"),
    ("複数選択可（Ctrl/Shift）。Groupはsummary CSVに保存されます。", "Multiple selection is supported (Ctrl/Shift). Group labels are saved in the summary CSV."),
    ("台形分割数", "Corridor lane count"),
    ("代表画像でZO-1処理を確認・調整", "Review / Adjust ZO-1 Processing on Representative Image"),
    ("選択サンプルを解析", "Analyze Selected Samples"),
    ("全保存済みを解析＋Group summary", "Analyze All Reviewed + Group Summary"),
    ("解析待機", "Ready for analysis"),
    ("核アノテーション", "Nuclei annotation"),
    ("核数", "Nuclei"),
    ("処理対象ファイル", "Input files"),
    ("MAX Projection 設定", "MAX Projection Settings"),
    ("チャンネル数:", "Number of channels:"),
    ("出力先:", "Output folder:"),
    ("選択", "Browse"),
    ("各元ファイルと同じ場所", "Same folder as each source file"),
    ("出力先が空欄の場合は各入力TIFFと同じフォルダに保存します。", "If the output folder is blank, each result is saved next to its source TIFF."),
    ("MAX画像を作成", "Create MAX Images"),
    ("出力フォルダを解析タブで開く", "Open Output Folder in Analysis Tab"),
    ("ログ", "Log"),
    ("TIFFファイルを選択", "Select TIFF Files"),
    ("フォルダ内TIFFを追加", "Add TIFFs from Folder"),
    ("リストをクリア", "Clear List"),
    ("想定ページ順: Z1-Ch1, Z1-Ch2, ... → Z2-Ch1, Z2-Ch2, ...\n2ch時の既定出力は <sample>_DAPI.tif / <sample>_ZO1.tif なので、そのまま解析タブで認識できます。", "Expected page order: Z1-Ch1, Z1-Ch2, ... → Z2-Ch1, Z2-Ch2, ...\nFor 2 channels, the default outputs are <sample>_DAPI.tif / <sample>_ZO1.tif, which are recognized directly by the Analysis tab."),
    ("ページ順: Z1-Ch1, Z1-Ch2, ... → Z2-Ch1, Z2-Ch2, ...", "Page order: Z1-Ch1, Z1-Ch2, ... → Z2-Ch1, Z2-Ch2, ..."),
    ("Mica TIFFを選択", "Select Mica TIFF"),
    ("TIFFフォルダを選択", "Select TIFF Folder"),
    ("出力フォルダを選択", "Select Output Folder"),
    ("入力なし", "No input"),
    ("TIFFファイルを選択してください。", "Please select TIFF files."),
    ("設定エラー", "Settings Error"),
    ("チャンネル数は2〜4を指定してください。", "Set the number of channels to 2–4."),
    ("完了", "Completed"),
    ("解析フォルダ", "Analysis Folder"),
    ("解析へ渡せる単一の出力フォルダがありません。", "There is no single output folder that can be passed to the Analysis tab."),
    ("処理エラー", "Processing Error"),
    ("未計算", "Not calculated"),
    ("処理中...", "Processing..."),
    ("Mouse: 画像外", "Mouse: outside image"),
    ("Mouse: 画像上にカーソルを置いてください", "Mouse: move the cursor over the image"),
    ("現在のZO-1前処理設定を本解析パラメータへ反映しました。", "The current ZO-1 preprocessing settings were applied to the main analysis parameters."),
    ("DAPI/ZO-1が揃った全サンプルが確認済み・保存済みです。", "All samples with paired DAPI/ZO-1 images have already been reviewed and saved."),
    ("最後のサンプルまで到達しました。", "Reached the last sample."),
    ("有効な核が2個未満です。", "Fewer than two valid nuclei were found."),
    ("核ラベル数がuint16上限に達しました。", "The number of nucleus labels reached the uint16 limit."),
    ("核数がuint16上限を超えました。", "The number of nuclei exceeds the uint16 limit."),
    ("核ラベル画像サイズが一致しません。", "The nucleus-label image size does not match."),
    ("核を1つだけ選択してください。", "Please select exactly one nucleus."),
    ("Fill holeする核を選択してください。", "Select a nucleus for Fill Hole."),
    ("選択核に埋める内部holeはありませんでした。", "No internal hole was found in the selected nucleus."),
    ("直接接している別ラベル核はありませんでした。", "No directly touching nucleus with a different label was found."),
    ("すべての核検出を削除しますか？", "Delete all nucleus detections?"),
    ("StarDist処理中は保存できません。", "Cannot save while StarDist is running."),
    ("DAPI核検出・修正結果を保存して閉じますか？", "Save the DAPI nucleus detection/review results and close?"),
    ("Brightness maxはminより大きくしてください。", "Brightness max must be greater than min."),
    ("1 pixel (µm)は0より大きくしてください。", "1 pixel (µm) must be greater than 0."),
    ("Probability thresholdは0～1です。", "Probability threshold must be between 0 and 1."),
    ("NMS thresholdは0～1です。", "NMS threshold must be between 0 and 1."),
    ("核面積は0以上にしてください。", "Nucleus area must be 0 or greater."),
    ("Max nucleus areaはMin nucleus areaより大きくしてください。", "Max nucleus area must be greater than Min nucleus area."),
    ("Normalize p-min < p-max、かつ0～100にしてください。", "Set Normalize p-min < p-max, both within 0–100."),
    ("Rolling ball radiusは1 px以上にしてください。", "Rolling-ball radius must be at least 1 px."),
    ("White top-hat radiusは1 px以上にしてください。", "White top-hat radius must be at least 1 px."),
    ("Top-hat thresholdは0以上にしてください。", "Top-hat threshold must be 0 or greater."),
    ("ZO-1 top-hat thresholdは0以上にしてください。", "ZO-1 top-hat threshold must be 0 or greater."),
    ("Min object sizeは1 px以上にしてください。", "Min object size must be at least 1 px."),
    ("DAPI postfix と ZO-1 postfix の両方を入力してください。", "Enter both the DAPI postfix and the ZO-1 postfix."),
    ("DAPI postfix と ZO-1 postfix は異なる文字列にしてください。", "The DAPI and ZO-1 postfixes must be different."),
    ("先に親フォルダを選択してください。", "Select a parent folder first."),
    ("未確認のDAPI画像はありません。", "There are no unreviewed DAPI images."),
    ("サンプルを選択してください。", "Select a sample."),
    ("代表画像として使うサンプルを選択してください。", "Select a sample to use as the representative image."),
    ("選択サンプルにZO-1画像がありません。", "The selected sample does not have a ZO-1 image."),
    ("台形分割数は4以上にしてください。", "Corridor lane count must be at least 4."),
    ("Neighbor distance factorは1.0以上にしてください。", "Neighbor distance factor must be at least 1.0."),
    ("ZO-1 rolling ball radiusは1 px以上にしてください。", "ZO-1 rolling-ball radius must be at least 1 px."),
    ("Lane peak top N pixelsは1以上にしてください。", "Lane peak top N pixels must be at least 1."),
    ("Intact cutoffは0～100%にしてください。", "Intact cutoff must be between 0 and 100%."),
    ("解析するサンプルを選択してください。", "Select samples to analyze."),
    ("DAPIまたはZO-1がありません。", "DAPI or ZO-1 is missing."),
    ("先にDAPI核をStarDistで確認・修正して保存してください。", "Review/correct the DAPI nuclei with StarDist and save them first."),
    ("保存済みDAPI核アノテーションのあるサンプルがありません。", "No samples have saved DAPI nucleus annotations."),
    ("Groupを設定するサンプルを選択してください。", "Select samples for which to set a Group."),
    ("Group名を入力してください。", "Enter a Group name."),
    ("Groupをクリアするサンプルを選択してください。", "Select samples whose Group should be cleared."),
    ("TIFFにページがありません。", "The TIFF contains no pages."),
    ("保存済みDAPI核アノテーションがありません:\n", "No saved DAPI nucleus annotation was found:\n"),
    ("保存済み核labelとDAPI画像サイズが一致しません。\nlabels=", "Saved nucleus-label size does not match the DAPI image.\nlabels="),
    ("DAPIとZO-1画像サイズが一致しません。\nDAPI=", "DAPI and ZO-1 image sizes do not match.\nDAPI="),
    ("保存済み核ラベルとDAPIサイズが一致しません。\n保存: ", "Saved nucleus-label size does not match the DAPI image.\nSaved: "),
    ("DAPIとZO-1画像サイズが一致しません。\nDAPI: ", "DAPI and ZO-1 image sizes do not match.\nDAPI: "),
    ("ZO-1 overlayとDAPI画像サイズが一致しません。\nDAPI: ", "ZO-1 overlay and DAPI image sizes do not match.\nDAPI: "),
    ("Brightness min/maxを整数で入力してください。", "Enter integer values for Brightness min/max."),
    ("現在の検出結果と手動修正を新しいStarDist結果で置き換えますか？", "Replace the current detections and manual edits with new StarDist results?"),
    ("パラメータを確認してください: ", "Check the parameters: "),
    ("解析設定を確認してください: ", "Check the analysis settings: "),
    ("代表画像", "Representative image"),
    ("再計算 / Apply parameters", "Recalculate / Apply Parameters"),
    ("← 前の処理", "← Previous Step"),
    ("次の処理 →", "Next Step →"),
    ("この設定を本解析へ反映", "Apply These Settings to Main Analysis"),
    ("反映後、親GUIの「選択サンプルを解析」または「全保存済みサンプルを解析」を実行してください。", "After applying, run 'Analyze Selected Samples' or 'Analyze All Reviewed' in the parent GUI."),
    ("右ドラッグ: pan\nMouse wheel: zoom\n各処理は同じviewportで比較できます。", "Right drag: pan\nMouse wheel: zoom\nAll processing steps can be compared in the same viewport."),
    ("Thresholdはwhite top-hat後のraw輝度スケールに適用します。", "The threshold is applied to the raw-intensity scale after white top-hat filtering."),
    ("保存 [S]", "Save [S]"),
    ("保存して閉じる", "Save and Close"),
    ("保存して次へ [N]", "Save and Next [N]"),
    ("左ドラッグ: 新しい核を描画\n1核選択後の左ドラッグ: その核へ追加\nCtrl+左ドラッグ: 消しゴム\nダブルクリック: 核選択/解除\n右ドラッグ: 移動\nホイール: 拡大縮小\n[ / ]: ブラシサイズ\nDelete: 選択核削除\nCtrl+Z: Undo", "Left drag: draw a new nucleus\nLeft drag after selecting one nucleus: add to that nucleus\nCtrl + left drag: eraser\nDouble click: select/deselect nucleus\nRight drag: pan\nMouse wheel: zoom\n[ / ]: brush size\nDelete: delete selected nucleus\nCtrl+Z: Undo"),
    ("Sample", "Sample"),
    ("Group", "Group"),
    ("StarDist detection", "StarDist Detection"),
    ("Run StarDist / replace detections", "Run StarDist / Replace Detections"),
    ("Display", "Display"),
    ("Show nucleus contours", "Show Nucleus Contours"),
    ("Show nucleus IDs", "Show Nucleus IDs"),
    ("Brightness min", "Brightness Min"),
    ("Brightness max", "Brightness Max"),
    ("Apply", "Apply"),
    ("Auto", "Auto"),
    ("Fit image", "Fit Image"),
    ("Manual correction", "Manual Correction"),
    ("Delete selected", "Delete Selected"),
    ("Merge selected", "Merge Selected"),
    ("Fill hole selected", "Fill Hole in Selected"),
    ("Merge touching nuclei (8-neighbor)", "Merge Touching Nuclei (8-neighbor)"),
    ("Split disconnected selected", "Split Disconnected Selected"),
    ("Undo", "Undo"),
    ("Clear all nuclei", "Clear All Nuclei"),
    ("Save", "Save"),
    ("Status / mouse", "Status / Mouse"),
    ("StarDist processing...", "StarDist Processing..."),
    ("StarDist error", "StarDist Error"),
    ("Manually edited (not saved)", "Manually Edited (Not Saved)"),
    ("ZO-1 preprocessing parameters", "ZO-1 Preprocessing Parameters"),
    ("Processing step", "Processing Step"),
    ("QC / statistics", "QC / Statistics"),
    ("Use settings", "Use Settings"),
    ("View", "View"),
    ("Tight-junction analysis", "Tight-junction Analysis"),
    ("Max pair distance factor", "Max Pair Distance Factor"),
    ("ZO-1 rolling ball radius (px)", "ZO-1 Rolling-ball Radius (px)"),
    ("White top-hat radius (px)", "White Top-hat Radius (px)"),
    ("Top-hat threshold", "Top-hat Threshold"),
    ("Min object size (px)", "Min Object Size (px)"),
    ("Intact cutoff (%)", "Intact Cutoff (%)"),
    ("Lane peak top N pixels", "Lane Peak Top N Pixels"),
    ("確認済み・保存済み", "Reviewed / Saved"),
    ("StarDist済み・未確認", "StarDist Draft / Unreviewed"),
    ("未処理", "Not Processed"),
    ("Batch StarDist draft loaded / 未確認", "Batch StarDist Draft Loaded / Unreviewed"),
    ("なし", "None"),
    ("DAPI StarDist核検出 / 確認", "DAPI StarDist detection / review"),
    ("ZO-1前処理QC", "ZO-1 preprocessing QC"),
    ("1. Raw ZO-1", "1. Raw ZO-1"),
    ("2. Rolling-ball background", "2. Rolling-ball background"),
    ("3. Background-subtracted", "3. Background-subtracted"),
    ("4. White top-hat", "4. White top-hat"),
    ("5. Threshold binary", "5. Threshold binary"),
    ("6. Skeleton", "6. Skeleton"),
    ("7. Raw + skeleton", "7. Raw + skeleton"),
    ("8. Top-hat + skeleton", "8. Top-hat + skeleton"),
    ("9. DAPI + nuclei + skeleton", "9. DAPI + nuclei + skeleton"),
]

# Additional Japanese labels for originally English-only UI text so Japanese mode
# can also provide a more complete localized interface.
_I18N_EN_TO_JA_EXTRA = {
    "Sample": "サンプル",
    "Group": "グループ",
    "StarDist Detection": "StarDist核検出",
    "Run StarDist / Replace Detections": "StarDist実行 / 検出を置換",
    "Display": "表示",
    "Show Nucleus Contours": "核輪郭を表示",
    "Show Nucleus IDs": "核IDを表示",
    "Brightness Min": "表示輝度 Min",
    "Brightness Max": "表示輝度 Max",
    "Apply": "適用",
    "Auto": "自動",
    "Fit Image": "画像全体を表示",
    "Manual Correction": "手動修正",
    "Delete Selected": "選択核を削除",
    "Merge Selected": "選択核をMerge",
    "Fill Hole in Selected": "選択核のHoleを埋める",
    "Merge Touching Nuclei (8-neighbor)": "接触核をMerge（8近傍）",
    "Split Disconnected Selected": "非連結の選択核を分割",
    "Undo": "元に戻す",
    "Clear All Nuclei": "全核を削除",
    "Save": "保存",
    "Status / Mouse": "状態 / マウス",
    "StarDist Processing...": "StarDist処理中...",
    "StarDist Error": "StarDistエラー",
    "Manually Edited (Not Saved)": "手動修正済み（未保存）",
    "Representative image": "代表画像",
    "ZO-1 Preprocessing Parameters": "ZO-1前処理パラメータ",
    "Processing Step": "処理ステップ",
    "QC / Statistics": "QC / 統計",
    "Use Settings": "設定を使用",
    "View": "表示",
    "Tight-junction Analysis": "Tight-junction解析",
    "Max Pair Distance Factor": "最大ペア距離係数",
    "ZO-1 Rolling-ball Radius (px)": "ZO-1 rolling ball半径 (px)",
    "White Top-hat Radius (px)": "White top-hat半径 (px)",
    "Top-hat Threshold": "Top-hat閾値",
    "Min Object Size (px)": "最小object size (px)",
    "Intact Cutoff (%)": "Intact cutoff (%)",
    "Lane Peak Top N Pixels": "Lane peak上位N pixel",
    "Step 1: Mica multipage TIFF → Channel MAX projection": "Step 1: Mica multipage TIFF → チャンネル別MAX投影",
    "Group:": "グループ:",
    "DAPI StarDist detection / review": "DAPI StarDist核検出 / 確認",
    "ZO-1 preprocessing QC": "ZO-1前処理QC",
    "1. Raw ZO-1": "1. Raw ZO-1",
    "2. Rolling-ball background": "2. Rolling-ball背景",
    "3. Background-subtracted": "3. Background subtraction後",
    "4. White top-hat": "4. White top-hat",
    "5. Threshold binary": "5. Threshold二値化",
    "6. Skeleton": "6. Skeleton",
    "7. Raw + skeleton": "7. Raw + skeleton",
    "8. Top-hat + skeleton": "8. Top-hat + skeleton",
    "9. DAPI + nuclei + skeleton": "9. DAPI + 核 + skeleton",
}

_JA_TO_EN = {ja: en for ja, en in _I18N_PAIRS}
_EN_TO_JA = {en: ja for ja, en in _I18N_PAIRS}
_EN_TO_JA.update(_I18N_EN_TO_JA_EXTRA)
for _en, _ja in _I18N_EN_TO_JA_EXTRA.items():
    _JA_TO_EN[_ja] = _en

# Dynamic fragments used inside f-strings/messages.
_DYNAMIC_JA_EN = [
    ("完了: 成功 ", "Completed: success "),
    (" ファイル / エラー ", " files / errors "),
    (" ファイル", " files"),
    ("ファイルの処理が完了しました。", " files were processed successfully."),
    ("成功 ", "Success "),
    ("、エラー ", ", errors "),
    ("。\n詳細はログを確認してください。", ".\nSee the log for details."),
    ("解析可能 ", "Ready "),
    ("確認済み ", "Reviewed "),
    ("StarDist未確認 ", "StarDist drafts "),
    ("不足: ", "Missing: "),
    ("自動ペアリング: ", "Automatic pairing: "),
    ("DAPI StarDist核検出 / 確認 - ", "DAPI StarDist detection / review - "),
    ("ZO-1前処理QC - ", "ZO-1 preprocessing QC - "),
    ("核数: ", "Nuclei: "),
    ("選択: ", "Selected: "),
    ("枚を連続StarDistします。", " images will be processed sequentially with StarDist."),
    ("枚", " images"),
    ("未確認DAPI ", "Unreviewed DAPI: "),
    ("現在の検出結果", "current detections"),
    ("処理中", "Processing"),
    ("保存済み核を読み込み中", "Loading saved nuclei"),
    ("ZO-1を読み込み中", "Loading ZO-1"),
    ("ZO-1二値化 / skeleton", "ZO-1 binarization / skeleton"),
    ("核ペアを作成中", "Creating nucleus pairs"),
    ("核ペア解析 ", "Nucleus-pair analysis "),
    ("overlay作成中", "Creating overlay"),
    ("解析開始: ", "Analysis started: "),
    ("解析終了: ", "Analysis completed: "),
    ("エラー: ", "Errors: "),
    ("サンプルのtight-junction解析が終了しました。", "Tight-junction analysis of the samples has completed."),
    ("を親フォルダに出力しました。", " was written to the parent folder."),
    ("ページサイズが不一致です: ", "Page size mismatch: "),
    ("総ページ数 ", "Total page count "),
    (" がチャンネル数 ", " is not divisible by the number of channels "),
    (" で割り切れません。", "."),
    ("ページ順またはチャンネル数を確認してください。", "Check the page order or number of channels."),
    ("各ページが2D画像ではありません", "Each page must be a 2D image"),
    ("この処理はMicaのmultipage grayscale TIFFを想定しています。", "This process expects a Mica multipage grayscale TIFF."),
    ("StarDistを読み込めません。", "Could not import StarDist."),
    ("必要です。", " is required."),
    ("実行しますか？", "Continue?"),
    ("開始しますか？", "Start?"),
]
_DYNAMIC_EN_JA = [
    (en, ja)
    for ja, en in _DYNAMIC_JA_EN
    if en not in {".", " files", " images"}
    and len(en.strip()) >= 4
]


def set_app_language(lang):
    global _APP_LANGUAGE
    _APP_LANGUAGE = "en" if str(lang).lower().startswith("en") else "ja"


def get_app_language():
    return _APP_LANGUAGE


def tr(text, lang=None):
    """Translate a UI string to the selected language without touching data values."""
    if text is None:
        return text
    if not isinstance(text, str):
        return text

    target = (lang or _APP_LANGUAGE)
    target = "en" if str(target).lower().startswith("en") else "ja"

    if target == "en":
        if text in _JA_TO_EN:
            return _JA_TO_EN[text]
        out = text
        for ja, en in sorted(_DYNAMIC_JA_EN, key=lambda x: len(x[0]), reverse=True):
            out = out.replace(ja, en)
        return out

    # Japanese mode: reverse both exact and dynamic translations. This allows
    # runtime switching back from English after widget strings were rewritten.
    if text in _EN_TO_JA:
        return _EN_TO_JA[text]
    out = text
    for en, ja in sorted(_DYNAMIC_EN_JA, key=lambda x: len(x[0]), reverse=True):
        out = out.replace(en, ja)
    return out


def _translate_treeview(tree, lang):
    try:
        cols = tree.cget("columns")
        for col in cols:
            try:
                h = tree.heading(col, "text")
                tree.heading(col, text=tr(h, lang))
            except Exception:
                pass

        # Translate status-like textual values, but only when a translation
        # actually changes the string. Sample names/file names remain untouched.
        for iid in tree.get_children(""):
            try:
                vals = list(tree.item(iid, "values"))
                changed = False
                for i, v in enumerate(vals):
                    if isinstance(v, str):
                        nv = tr(v, lang)
                        if nv != v:
                            vals[i] = nv
                            changed = True
                if changed:
                    tree.item(iid, values=vals)
            except Exception:
                pass
    except Exception:
        pass


def translate_widget_tree(widget, lang=None):
    """Translate existing Tk/ttk widgets, including pop-up editor windows."""
    target = lang or _APP_LANGUAGE

    try:
        if isinstance(widget, (tk.Tk, tk.Toplevel)):
            title = widget.title()
            if title:
                widget.title(tr(title, target))
    except Exception:
        pass

    # Widgets with direct text options.
    try:
        keys = widget.keys()
    except Exception:
        keys = ()

    if "text" in keys:
        try:
            # Avoid rewriting Entry/Text user content: only widgets that actually
            # expose a label-style text option enter this branch.
            current = widget.cget("text")
            new = tr(current, target)
            if new != current:
                widget.configure(text=new)
        except Exception:
            pass

    # Labels driven by StringVar (status labels etc.). Do not touch Entry values.
    if isinstance(widget, (tk.Label, ttk.Label)):
        try:
            tv = widget.cget("textvariable")
            if tv:
                current = widget.getvar(tv)
                if isinstance(current, str):
                    new = tr(current, target)
                    if new != current:
                        widget.setvar(tv, new)
        except Exception:
            pass

    if isinstance(widget, ttk.Treeview):
        _translate_treeview(widget, target)

    if isinstance(widget, ttk.Notebook):
        try:
            for tab_id in widget.tabs():
                current = widget.tab(tab_id, "text")
                new = tr(current, target)
                if new != current:
                    widget.tab(tab_id, text=new)
        except Exception:
            pass

    if isinstance(widget, ttk.Combobox):
        try:
            vals = list(widget.cget("values"))
            new_vals = [tr(str(v), target) for v in vals]
            if new_vals != vals:
                widget.configure(values=new_vals)
            current = widget.get()
            new_current = tr(current, target)
            if new_current != current:
                widget.set(new_current)
        except Exception:
            pass

    if isinstance(widget, tk.Listbox):
        try:
            for i in range(widget.size()):
                current = widget.get(i)
                new = tr(current, target)
                if new != current:
                    widget.delete(i)
                    widget.insert(i, new)
        except Exception:
            pass

    try:
        for child in widget.winfo_children():
            translate_widget_tree(child, target)
    except Exception:
        pass


# Translate message boxes/file-dialog titles automatically, including strings
# generated by the existing analysis code after the language has been switched.
_ORIG_MESSAGEBOX = {}
for _name in ("showinfo", "showwarning", "showerror", "askyesno", "askokcancel", "askretrycancel", "askquestion"):
    if hasattr(messagebox, _name):
        _ORIG_MESSAGEBOX[_name] = getattr(messagebox, _name)


def _make_messagebox_wrapper(name):
    original = _ORIG_MESSAGEBOX[name]
    def wrapped(title, message=None, *args, **kwargs):
        title = tr(title)
        if message is not None:
            message = tr(message)
        return original(title, message, *args, **kwargs)
    return wrapped


for _name in _ORIG_MESSAGEBOX:
    setattr(messagebox, _name, _make_messagebox_wrapper(_name))

_ORIG_FILEDIALOG = {}
for _name in ("askopenfilename", "askopenfilenames", "askdirectory", "asksaveasfilename"):
    if hasattr(filedialog, _name):
        _ORIG_FILEDIALOG[_name] = getattr(filedialog, _name)


def _make_filedialog_wrapper(name):
    original = _ORIG_FILEDIALOG[name]
    def wrapped(*args, **kwargs):
        if "title" in kwargs:
            kwargs["title"] = tr(kwargs["title"])
        return original(*args, **kwargs)
    return wrapped


for _name in _ORIG_FILEDIALOG:
    setattr(filedialog, _name, _make_filedialog_wrapper(_name))


try:
    from scipy.ndimage import find_objects, binary_fill_holes
    from scipy.spatial import Delaunay, distance, cKDTree
except Exception:
    find_objects = None
    binary_fill_holes = None
    Delaunay = None
    distance = None
    cKDTree = None

try:
    from skimage.morphology import remove_small_objects, skeletonize
except Exception:
    remove_small_objects = None
    skeletonize = None


# ============================================================
# Settings
# ============================================================

DAPI_FILENAME = "1_CH1.tif"
ZO1_FILENAME = "1_CH3.tif"

OUTPUT_DIRNAME = "ZO1_junction_analysis"
ANNOTATION_FILENAME = "dapi_nuclei_annotation.npz"

# Batch StarDist result before human review.
# Keeping this separate prevents automatic batch detection from being mistaken
# for a manually reviewed annotation.
DRAFT_ANNOTATION_FILENAME = "dapi_nuclei_stardist_draft.npz"

LABEL_FILENAME = "dapi_nuclei_labels.tif"
BINARY_FILENAME = "dapi_nuclei_binary.tif"
OVERLAY_FILENAME = "dapi_nuclei_review_overlay.png"

# Legacy per-parent settings filename (read only once for migration).
SETTINGS_FILENAME = "zo1_nuclei_review_settings.json"

# Current settings are shared across ALL parent folders.
# Stored in the user's home directory so switching datasets does not change them.
GLOBAL_SETTINGS_FILENAME = "caco2_zo1_gui_global_settings.json"

ERROR_LOG_FILENAME = "ZO1_open_error.log"

DEFAULT_PIXEL_SIZE_UM = 0.37744
DEFAULT_PROB_THRESH = 0.50
DEFAULT_NMS_THRESH = 0.40
DEFAULT_MIN_AREA_UM2 = 15.0
DEFAULT_MAX_AREA_UM2 = 1000.0
DEFAULT_NORMALIZE_PMIN = 1.0
DEFAULT_NORMALIZE_PMAX = 99.8

# ============================================================
# Tight-junction analysis defaults
# ============================================================

DEFAULT_JUNCTION_LANES = 20

# Candidate neighbor filtering only. This is NOT the junction metric.
# Candidate pairs are ALL centroid pairs within:
# median nearest-neighbor distance * factor.
# A slightly wider default allows junctions to be measured across larger cells/gaps.
DEFAULT_PAIR_DISTANCE_FACTOR = 3.0

# Strict blocker rule:
# if even one pixel belonging to a THIRD nucleus falls inside the center-to-center
# pair corridor, that pair is excluded.
DEFAULT_MAX_OTHER_NUCLEUS_FRACTION = 0.0

# ZO-1 preprocessing:
# raw image -> rolling-ball-style background subtraction -> fixed RAW threshold
#
# Radius is expressed in ORIGINAL image pixels.
DEFAULT_ZO1_ROLLING_BALL_RADIUS_PX = 200

# After rolling-ball subtraction, enhance thin bright junctional ridges
# using a white top-hat filter. Radius is in ORIGINAL image pixels.
DEFAULT_ZO1_TOPHAT_RADIUS_PX = 7

# Fixed threshold applied to the WHITE TOP-HAT image on the original
# raw-intensity scale (not on a 0-255 normalized image).
DEFAULT_ZO1_TOPHAT_THRESHOLD = 15.0

DEFAULT_ZO1_MIN_OBJECT_SIZE_PX = 8

# After small-object removal, a lane is positive if at least this many
# ZO-1 skeleton pixels are present inside the lane.
DEFAULT_LANE_MIN_SKELETON_PIXELS = 1

# Junction intensity metric:
# for each lane, average the brightest N pixels in the ROLLING-BALL
# background-subtracted ZO-1 image. Then average across all lanes.
# This does not require the lane to pass the skeleton threshold.
DEFAULT_LANE_PEAK_TOP_N_PIXELS = 5

# Summary-only classification threshold.
DEFAULT_INTACT_COVERAGE_THRESHOLD_PCT = 70.0

PAIR_CSV_FILENAME = "junction_pair_analysis.csv"
JUNCTION_SUMMARY_FILENAME = "junction_analysis_summary.csv"
JUNCTION_OVERLAY_FILENAME = "junction_analysis_overlay.png"
ZO1_BG_SUB_FILENAME = "zo1_background_subtracted.tif"
ZO1_TOPHAT_FILENAME = "zo1_white_tophat.tif"
ZO1_BINARY_FILENAME = "zo1_binary_mask.tif"
ZO1_SKELETON_FILENAME = "zo1_skeleton.tif"

# QC overlays requested for checking the ZO-1 skeleton.
ZO1_RAW_SKELETON_OVERLAY_FILENAME = "zo1_raw_plus_skeleton_overlay.png"
ZO1_TOPHAT_SKELETON_OVERLAY_FILENAME = "zo1_tophat_plus_skeleton_overlay.png"
ZO1_NUCLEI_SKELETON_OVERLAY_FILENAME = "zo1_skeleton_plus_nuclei_overlay.png"

ALL_JUNCTION_SUMMARY_FILENAME = "ZO1_junction_summary.csv"
GROUP_JUNCTION_SUMMARY_FILENAME = "ZO1_junction_summary_by_group.csv"
SAMPLE_GROUPS_FILENAME = "zo1_sample_groups.json"


# ============================================================
# Image I/O
# ============================================================

def _first_2d_plane(arr):
    arr = np.asarray(arr)
    arr = np.squeeze(arr)

    if arr.ndim == 2:
        return arr

    # RGB/RGBA
    if arr.ndim == 3 and arr.shape[-1] in (3, 4):
        rgb = arr[..., :3].astype(np.float32)
        return (
            0.299 * rgb[..., 0]
            + 0.587 * rgb[..., 1]
            + 0.114 * rgb[..., 2]
        )

    # Z/T stack -> first plane after squeezing
    while arr.ndim > 2:
        arr = arr[0]

    if arr.ndim != 2:
        raise ValueError(f"2D画像へ変換できません: shape={arr.shape}")

    return arr


def read_tiff_gray(path):
    """
    高速なOpenCV decodeを最初に試す。
    Unicode path対策として np.fromfile + cv2.imdecode を使用。
    失敗時は tifffile -> Pillow にfallback。
    """
    path = Path(path)
    errors = []

    # 1) OpenCV
    try:
        buf = np.fromfile(str(path), dtype=np.uint8)
        img = cv2.imdecode(buf, cv2.IMREAD_UNCHANGED)
        if img is not None:
            if img.ndim == 3:
                if img.shape[-1] >= 3:
                    img = cv2.cvtColor(img[..., :3], cv2.COLOR_BGR2GRAY)
                else:
                    img = _first_2d_plane(img)
            img = _first_2d_plane(img)
            return np.asarray(img)
    except Exception as exc:
        errors.append(f"OpenCV: {exc}")

    # 2) tifffile
    try:
        with tifffile.TiffFile(str(path)) as tif:
            if len(tif.series) > 0:
                arr = tif.series[0].asarray()
            else:
                arr = tif.asarray()
        return np.asarray(_first_2d_plane(arr))
    except Exception as exc:
        errors.append(f"tifffile: {exc}")

    # 3) Pillow
    try:
        with Image.open(str(path)) as im:
            im.seek(0)
            arr = np.asarray(im)
        return np.asarray(_first_2d_plane(arr))
    except Exception as exc:
        errors.append(f"Pillow: {exc}")

    raise RuntimeError(
        f"TIFFを読み込めません:\n{path}\n\n" + "\n".join(errors)
    )


def normalize_8bit(img):
    x = np.asarray(img, dtype=np.float32)
    good = x[np.isfinite(x)]
    if good.size == 0:
        return np.zeros(x.shape, dtype=np.uint8)

    p1, p99 = np.percentile(good, [1.0, 99.5])
    if p99 <= p1:
        p1, p99 = float(good.min()), float(good.max())

    if p99 <= p1:
        return np.zeros(x.shape, dtype=np.uint8)

    out = np.clip((x - p1) / (p99 - p1), 0.0, 1.0)
    return np.round(out * 255.0).astype(np.uint8)


def imwrite_unicode(path, image, params=None):
    """
    Unicode-safe image writer for Windows.

    cv2.imwrite() can silently return False for paths containing Japanese
    characters on some OpenCV/Windows combinations. Encode in memory first,
    then write bytes with numpy.tofile(), which handles Unicode paths.

    Raises RuntimeError if encoding or writing fails.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    ext = path.suffix.lower()
    if ext == ".jpeg":
        ext = ".jpg"

    if ext not in {".png", ".jpg", ".bmp", ".tif", ".tiff", ".webp"}:
        raise ValueError(f"Unsupported image extension: {path.suffix}")

    if params is None:
        params = []

    ok, encoded = cv2.imencode(
        ext,
        np.asarray(image),
        params
    )

    if not ok:
        raise RuntimeError(
            f"OpenCV image encoding failed: {path}"
        )

    try:
        encoded.tofile(str(path))
    except Exception as exc:
        raise RuntimeError(
            f"Image write failed: {path}\n{exc}"
        ) from exc

    if not path.exists() or path.stat().st_size <= 0:
        raise RuntimeError(
            f"Image file was not created correctly: {path}"
        )

    return True


# ============================================================
# File helpers
# ============================================================

def find_named_tiff(folder, wanted):
    folder = Path(folder)

    p = folder / wanted
    if p.exists():
        return p

    wanted_stem = Path(wanted).stem.lower()

    for q in folder.iterdir():
        if (
            q.is_file()
            and q.suffix.lower() in {".tif", ".tiff"}
            and q.stem.lower() == wanted_stem
        ):
            return q

    return None


def _analysis_paths_from_dir(d):
    d = Path(d)
    return {
        "dir": d,
        "npz": d / ANNOTATION_FILENAME,
        "draft_npz": d / DRAFT_ANNOTATION_FILENAME,
        "labels": d / LABEL_FILENAME,
        "binary": d / BINARY_FILENAME,
        "overlay": d / OVERLAY_FILENAME,
    }


def output_paths(folder):
    """
    Legacy child-folder layout:
        sample_folder/ZO1_junction_analysis/...
    """
    return _analysis_paths_from_dir(
        Path(folder) / OUTPUT_DIRNAME
    )


def _safe_name_for_path(name):
    s = str(name).strip()
    s = re.sub(r'[<>:"/\\|?*]+', "_", s)
    s = s.rstrip(" .")
    return s or "sample"


def sample_display_name(sample):
    """
    Flat-file mode uses sample['name'] (e.g. con1).
    Legacy mode falls back to the child-folder name.
    """
    if isinstance(sample, dict):
        name = str(sample.get("name", "")).strip()
        if name:
            return name

        folder = sample.get("folder")
        if folder is not None:
            return Path(folder).name

    return str(sample)


def sample_output_paths(sample):
    """
    Flat input:
        ROOT/con1_DAPI.tiff
        ROOT/con1_ZO1.tiff

    Output:
        ROOT/ZO1_junction_analysis/con1/...

    Legacy child-folder input is still supported.
    """
    if (
        isinstance(sample, dict)
        and sample.get("flat_input", False)
    ):
        root = Path(sample["folder"])
        name = _safe_name_for_path(
            sample_display_name(sample)
        )
        return _analysis_paths_from_dir(
            root / OUTPUT_DIRNAME / name
        )

    if isinstance(sample, dict):
        return output_paths(sample["folder"])

    return output_paths(sample)


def natural_sort_key(text):
    """
    con1, con2, con10 の順に並べるための自然順ソート。
    """
    return [
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", str(text))
    ]


def split_sample_channel_by_postfix(stem, dapi_postfix, zo1_postfix):
    """Return (sample_name, channel) for a TIFF stem using user postfixes.

    Matching is case-insensitive and the postfix is treated literally.
    Example:
        stem='con1_DAPI_MAX', dapi_postfix='_DAPI_MAX'
        -> ('con1', 'DAPI')
    """
    stem = str(stem)
    dapi_pf = str(dapi_postfix or "").strip()
    zo1_pf = str(zo1_postfix or "").strip()

    if not dapi_pf or not zo1_pf:
        return None

    if dapi_pf.casefold() == zo1_pf.casefold():
        return None

    folded = stem.casefold()

    for channel, postfix in (("DAPI", dapi_pf), ("ZO1", zo1_pf)):
        pf_fold = postfix.casefold()
        if folded.endswith(pf_fold):
            sample = stem[: len(stem) - len(postfix)].strip()
            if sample:
                return sample, channel

    return None


def load_json(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            x = json.load(f)
        return x if isinstance(x, dict) else {}
    except Exception:
        return {}


def save_json(path, x):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(x, f, indent=2, ensure_ascii=False)


def save_stardist_draft(
    path,
    labels,
    sample_name,
    dapi_name,
    zo1_name,
    params,
    touching_merge_groups=0,
    touching_labels_removed=0,
):
    """
    Save an automatic StarDist result that has NOT yet been human reviewed.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    labels = np.asarray(
        labels,
        dtype=np.uint16
    )

    np.savez_compressed(
        path,
        nuclei_labels=labels,
        nuclei_binary=(labels > 0).astype(np.uint8) * 255,
        sample_folder=str(sample_name),
        dapi_filename=str(dapi_name),
        zo1_filename=str(zo1_name),
        pixel_size_um=float(params["pixel_size_um"]),
        prob_thresh=float(params["prob_thresh"]),
        nms_thresh=float(params["nms_thresh"]),
        min_area_um2=float(params["min_area_um2"]),
        max_area_um2=float(params["max_area_um2"]),
        normalize_pmin=float(params["normalize_pmin"]),
        normalize_pmax=float(params["normalize_pmax"]),
        auto_merge_touching_8neighbor=True,
        touching_merge_groups=int(touching_merge_groups),
        touching_labels_removed=int(touching_labels_removed),
        review_status="stardist_draft_unreviewed",
    )


# ============================================================
# Tight-junction analysis
# ============================================================

def _write_csv(path, rows, fieldnames=None):
    path = Path(path)

    if fieldnames is None:
        if rows:
            fieldnames = list(rows[0].keys())
        else:
            fieldnames = []

    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_group_summary(path, summaries):
    """
    Sample-level summaryを群ごとに集約する。

    統計単位は sample（画像）であり、pairを独立nとして扱わない。
    主なendpointについて groupごとの mean / SD を出力する。
    """
    metrics = [
        "mean_junction_coverage_pct",
        "median_junction_coverage_pct",
        "mean_lane_peak_intensity",
        "median_pair_lane_peak_intensity",
        "intact_pair_fraction",
        "n_valid_pairs",
        "n_reviewed_nuclei",
    ]

    grouped = {}

    for row in summaries:
        group = str(
            row.get("group", "")
        ).strip() or "Unassigned"

        grouped.setdefault(
            group,
            []
        ).append(row)

    out_rows = []

    for group in sorted(
        grouped,
        key=natural_sort_key
    ):
        items = grouped[group]

        out = {
            "group": group,
            "n_samples": len(items),
        }

        for metric in metrics:
            vals = []

            for item in items:
                try:
                    v = float(
                        item.get(metric, np.nan)
                    )
                except Exception:
                    v = np.nan

                if np.isfinite(v):
                    vals.append(v)

            arr = np.asarray(
                vals,
                dtype=np.float64
            )

            out[f"{metric}_mean"] = (
                float(np.mean(arr))
                if arr.size
                else np.nan
            )

            out[f"{metric}_sd"] = (
                float(np.std(arr, ddof=1))
                if arr.size > 1
                else np.nan
            )

        out_rows.append(out)

    _write_csv(
        path,
        out_rows
    )


def _nucleus_regions_from_labels(labels):
    """
    Build compact per-nucleus geometry using bbox-local pixels only.

    Returns:
        regions[label] = {
            bbox: (x1,y1,x2,y2),
            centroid_xy: [x,y],
            coords_xy: Nx2 float array
        }
    """
    labels = np.asarray(labels)
    regions = {}

    if find_objects is not None:
        objects = find_objects(labels)
        iterable = enumerate(objects, start=1)
    else:
        labs = [int(x) for x in np.unique(labels) if int(x) > 0]
        iterable = []
        for lab in labs:
            ys, xs = np.where(labels == lab)
            if len(xs) == 0:
                continue
            slc = (
                slice(int(ys.min()), int(ys.max()) + 1),
                slice(int(xs.min()), int(xs.max()) + 1)
            )
            iterable.append((lab, slc))

    for lab, slc in iterable:
        if slc is None:
            continue

        ys, xs = slc
        x1, x2 = int(xs.start), int(xs.stop)
        y1, y2 = int(ys.start), int(ys.stop)

        local = labels[y1:y2, x1:x2] == int(lab)
        yy, xx = np.where(local)

        if len(xx) == 0:
            continue

        gx = xx.astype(np.float64) + x1
        gy = yy.astype(np.float64) + y1

        coords = np.column_stack([gx, gy])

        regions[int(lab)] = {
            "bbox": (x1, y1, x2, y2),
            "centroid_xy": np.array(
                [float(np.mean(gx)), float(np.mean(gy))],
                dtype=np.float64
            ),
            "coords_xy": coords,
            "area_px": int(len(xx)),
        }

    return regions


def _delaunay_edges(centroids_xy):
    n = len(centroids_xy)

    if n < 2:
        return set()

    if n == 2:
        return {(0, 1)}

    if Delaunay is None:
        raise RuntimeError(
            "Tight-junction解析には scipy が必要です。\n"
            "pip install scipy"
        )

    try:
        tri = Delaunay(centroids_xy)
    except Exception:
        # Degenerate / near-collinear fallback:
        # connect each point to its two closest neighbors.
        if distance is None:
            raise

        D = distance.cdist(centroids_xy, centroids_xy)
        np.fill_diagonal(D, np.inf)

        edges = set()

        for i in range(n):
            order = np.argsort(D[i])[: min(2, n - 1)]
            for j in order:
                a, b = sorted((int(i), int(j)))
                edges.add((a, b))

        return edges

    edges = set()

    for simplex in tri.simplices:
        for i, j in ((0, 1), (1, 2), (2, 0)):
            a, b = int(simplex[i]), int(simplex[j])
            if a > b:
                a, b = b, a
            edges.add((a, b))

    return edges


def _adaptive_pair_distance(centroids_xy, factor):
    n = len(centroids_xy)

    if n < 2:
        return np.inf

    if distance is None:
        return np.inf

    D = distance.cdist(centroids_xy, centroids_xy)
    np.fill_diagonal(D, np.inf)

    nn = np.min(D, axis=1)
    nn = nn[np.isfinite(nn)]

    if nn.size == 0:
        return np.inf

    return float(np.median(nn) * float(factor))


def _candidate_pairs_within_distance(
    centroids_xy,
    max_distance
):
    """
    Return ALL nucleus pairs whose centroid distance <= max_distance.

    This intentionally replaces Delaunay-only candidate generation so that
    somewhat farther but unobstructed cell pairs can still be analyzed.
    cKDTree avoids constructing a full NxN distance matrix for candidate search.
    """
    pts = np.asarray(
        centroids_xy,
        dtype=np.float64
    )

    n = len(pts)

    if n < 2:
        return set()

    if not np.isfinite(max_distance):
        return {
            (i, j)
            for i in range(n)
            for j in range(i + 1, n)
        }

    if cKDTree is not None:
        tree = cKDTree(pts)

        pairs = tree.query_pairs(
            r=float(max_distance),
            output_type="set"
        )

        return {
            (int(i), int(j))
            if int(i) < int(j)
            else (int(j), int(i))
            for i, j in pairs
        }

    # Fallback when cKDTree is unavailable.
    if distance is None:
        raise RuntimeError(
            "Tight-junction解析には scipy が必要です。\\n"
            "pip install scipy"
        )

    D = distance.cdist(
        pts,
        pts
    )

    rows, cols = np.where(
        np.triu(
            D <= float(max_distance),
            k=1
        )
    )

    return {
        (int(i), int(j))
        for i, j in zip(rows, cols)
    }


def _build_pair_corridor(region_a, region_b, image_shape):
    """
    Approximate the two outer common tangents of two nuclear ROIs.

    Pair coordinates:
        u = nucleus A -> nucleus B
        v = perpendicular to u

    The longitudinal limits are the two cross-sections passing through the
    nuclear centroids (A center: u=0, B center: u=d).

    The min/max support of each nucleus along v are linearly connected to form
    the two outer tangent-like boundaries. This gives a quadrilateral/trapezoid
    that extends from nucleus-center to nucleus-center.
    """
    cA = region_a["centroid_xy"]
    cB = region_b["centroid_xy"]

    vec = cB - cA
    d = float(np.linalg.norm(vec))

    if d < 1e-6:
        return None

    u = vec / d
    v = np.array([-u[1], u[0]], dtype=np.float64)

    ptsA = region_a["coords_xy"]
    ptsB = region_b["coords_xy"]

    relA = ptsA - cA
    relB = ptsB - cB

    uA = relA @ u
    vA = relA @ v

    uB_local = relB @ u
    vB = relB @ v

    # B positions expressed from A origin.
    uB_global = d + uB_local

    # Longitudinal analysis range:
    # use the line through nucleus-A center (u=0) to the line through
    # nucleus-B center (u=d), rather than the nearest facing nuclear surfaces.
    #
    # This intentionally allows a ZO-1 junction that projects into a large
    # nuclear ROI to remain inside the analysis corridor.
    x_start = 0.0
    x_end = float(d)

    if x_end <= x_start + 1.0:
        return None

    a_low = float(np.min(vA))
    a_high = float(np.max(vA))
    b_low = float(np.min(vB))
    b_high = float(np.max(vB))

    def low_at(x):
        return a_low + (b_low - a_low) * (x / d)

    def high_at(x):
        return a_high + (b_high - a_high) * (x / d)

    corners_pair = np.array(
        [
            [x_start, low_at(x_start)],
            [x_end, low_at(x_end)],
            [x_end, high_at(x_end)],
            [x_start, high_at(x_start)],
        ],
        dtype=np.float64
    )

    corners_xy = (
        cA
        + corners_pair[:, 0:1] * u
        + corners_pair[:, 1:2] * v
    )

    h, w = image_shape

    minx = max(
        0,
        int(np.floor(np.min(corners_xy[:, 0]))) - 2
    )
    maxx = min(
        w - 1,
        int(np.ceil(np.max(corners_xy[:, 0]))) + 2
    )
    miny = max(
        0,
        int(np.floor(np.min(corners_xy[:, 1]))) - 2
    )
    maxy = min(
        h - 1,
        int(np.ceil(np.max(corners_xy[:, 1]))) + 2
    )

    if maxx <= minx or maxy <= miny:
        return None

    yy, xx = np.mgrid[
        miny:maxy + 1,
        minx:maxx + 1
    ]

    px = xx.astype(np.float64) - cA[0]
    py = yy.astype(np.float64) - cA[1]

    xp = px * u[0] + py * u[1]
    yp = px * v[0] + py * v[1]

    low = a_low + (b_low - a_low) * (xp / d)
    high = a_high + (b_high - a_high) * (xp / d)

    width = high - low

    corridor = (
        (xp >= x_start)
        & (xp <= x_end)
        & (width > 1e-6)
        & (yp >= low)
        & (yp <= high)
    )

    t = np.full(
        xp.shape,
        np.nan,
        dtype=np.float32
    )

    t[corridor] = (
        (yp[corridor] - low[corridor])
        / width[corridor]
    ).astype(np.float32)

    return {
        "cA": cA,
        "cB": cB,
        "distance_px": d,
        "x_start": x_start,
        "x_end": x_end,
        "corners_xy": corners_xy,
        "bbox": (minx, miny, maxx, maxy),
        "mask_local": corridor,
        "lane_fraction_local": t,
    }


def _rolling_ball_background_fast(
    image,
    radius_px=DEFAULT_ZO1_ROLLING_BALL_RADIUS_PX,
):
    """
    Fast rolling-ball-style background estimation.

    A circular/elliptical morphological opening is used as a fast practical
    rolling-ball equivalent. For large radii the image is internally downsampled
    so the effective structuring-element radius stays <= 50 px, then the
    estimated background is resized back to the original resolution.

    Importantly, radius_px is always interpreted in ORIGINAL-image pixels.
    The returned background remains on the original raw intensity scale.
    """
    raw = np.asarray(image, dtype=np.float32)

    radius = max(1, int(round(float(radius_px))))
    h, w = raw.shape

    # Keep the internal opening kernel reasonably small for speed.
    # Example: radius 200 px -> downsample x4 -> radius 50 px internally.
    ds = max(1, int(np.ceil(radius / 50.0)))

    if ds > 1:
        small_w = max(8, int(round(w / ds)))
        small_h = max(8, int(round(h / ds)))

        work = cv2.resize(
            raw,
            (small_w, small_h),
            interpolation=cv2.INTER_AREA
        )

        r_work = max(
            1,
            int(round(radius / ds))
        )
    else:
        work = raw
        r_work = radius

    k = 2 * r_work + 1

    kernel = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE,
        (k, k)
    )

    background_small = cv2.morphologyEx(
        work,
        cv2.MORPH_OPEN,
        kernel
    )

    if ds > 1:
        background = cv2.resize(
            background_small,
            (w, h),
            interpolation=cv2.INTER_LINEAR
        )
    else:
        background = background_small

    return background.astype(
        np.float32,
        copy=False
    )


def _white_tophat_raw(
    image,
    radius_px=DEFAULT_ZO1_TOPHAT_RADIUS_PX,
):
    """
    White top-hat on a RAW-scale image:
        top_hat = image - opening(image)

    A circular/elliptical structuring element is used.
    Thin bright ZO-1 ridges are preserved, while broader cytoplasmic brightness
    is strongly suppressed.
    """
    x = np.asarray(image, dtype=np.float32)

    radius = max(
        1,
        int(round(float(radius_px)))
    )

    k = 2 * radius + 1

    kernel = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE,
        (k, k)
    )

    opened = cv2.morphologyEx(
        x,
        cv2.MORPH_OPEN,
        kernel
    )

    top_hat = x - opened

    np.maximum(
        top_hat,
        0,
        out=top_hat
    )

    return top_hat.astype(
        np.float32,
        copy=False
    )


def _preprocess_zo1_for_junction(
    zo1_raw,
    rolling_ball_radius_px=DEFAULT_ZO1_ROLLING_BALL_RADIUS_PX,
    tophat_radius_px=DEFAULT_ZO1_TOPHAT_RADIUS_PX,
    tophat_threshold=DEFAULT_ZO1_TOPHAT_THRESHOLD,
    min_object_size=DEFAULT_ZO1_MIN_OBJECT_SIZE_PX,
):
    """
    ZO-1 preprocessing for junction detection.

    Analysis scale:
        original RAW TIFF intensity

    Steps:
        1. rolling-ball-style background subtraction (default radius 200 px)
        2. white top-hat on the background-subtracted image (default radius 7 px)
        3. fixed threshold on the top-hat RAW scale (default 15)
        4. remove tiny objects
        5. skeletonize

    Why:
        Broad weak cytoplasmic ZO-1 signal can exceed a simple threshold and
        create a medial-axis skeleton through the cell center. White top-hat
        suppresses broad bright regions and preferentially keeps thin bright
        junctional ridges.

    8-bit normalized images are produced only for QC display.
    """
    if remove_small_objects is None or skeletonize is None:
        raise RuntimeError(
            "Tight-junction解析には scikit-image が必要です.\n"
            "pip install scikit-image"
        )

    raw = np.asarray(
        zo1_raw,
        dtype=np.float32
    )

    background = _rolling_ball_background_fast(
        raw,
        radius_px=rolling_ball_radius_px
    )

    bg_sub = raw - background

    np.maximum(
        bg_sub,
        0,
        out=bg_sub
    )

    top_hat = _white_tophat_raw(
        bg_sub,
        radius_px=tophat_radius_px
    )

    threshold = float(
        tophat_threshold
    )

    if threshold < 0:
        raise ValueError(
            "ZO-1 top-hat thresholdは0以上にしてください。"
        )

    binary = top_hat >= threshold

    binary = remove_small_objects(
        binary,
        min_size=max(
            1,
            int(min_object_size)
        )
    )

    skel = skeletonize(binary)

    # Visualization only.
    raw8 = normalize_8bit(raw)
    bg_sub8 = normalize_8bit(bg_sub)
    top_hat8 = normalize_8bit(top_hat)

    return (
        raw8,
        bg_sub.astype(np.float32, copy=False),
        bg_sub8,
        top_hat.astype(np.float32, copy=False),
        top_hat8,
        binary.astype(bool),
        skel.astype(bool),
        threshold
    )


def _coverage_color_bgr(coverage):
    """
    red -> yellow -> green for 0 -> 50 -> 100%.
    """
    x = float(np.clip(coverage, 0, 100)) / 100.0

    if x < 0.5:
        r = 255
        g = int(round(510 * x))
    else:
        g = 255
        r = int(round(510 * (1.0 - x)))

    return (0, g, r)


def analyze_tight_junction_sample(
    sample_folder,
    params,
    progress_callback=None,
):
    """
    Analyze one saved/reviewed DAPI annotation + ZO-1 image.

    sample_folder can be either:
      - legacy Path / child folder
      - flat-file sample dict created by MainApp.scan_samples()

    Primary output:
        junction_coverage_pct per adjacent nucleus pair
    """
    sample_group = ""

    if isinstance(sample_folder, dict):
        sample = sample_folder
        folder = Path(sample["folder"])
        display_name = sample_display_name(sample)
        sample_group = str(
            sample.get("group", "")
        ).strip()

        paths = sample_output_paths(sample)

        dapi_path = (
            Path(sample["dapi"])
            if sample.get("dapi") is not None
            else None
        )
        zo1_path = (
            Path(sample["zo1"])
            if sample.get("zo1") is not None
            else None
        )

    else:
        folder = Path(sample_folder)
        display_name = folder.name
        paths = output_paths(folder)

        dapi_path = find_named_tiff(
            folder,
            DAPI_FILENAME
        )
        zo1_path = find_named_tiff(
            folder,
            ZO1_FILENAME
        )

    if dapi_path is None or zo1_path is None:
        raise FileNotFoundError(
            f"DAPI/ZO-1 pairがありません: {display_name}"
        )

    if not paths["npz"].exists():
        raise FileNotFoundError(
            f"保存済みDAPI核アノテーションがありません:\n"
            f"{paths['npz']}"
        )

    if progress_callback:
        progress_callback(
            f"{display_name}: 保存済み核を読み込み中"
        )

    data = np.load(
        paths["npz"],
        allow_pickle=True
    )

    if "nuclei_labels" not in data:
        raise RuntimeError(
            f"nuclei_labelsがありません: {paths['npz'].name}"
        )

    labels = data["nuclei_labels"].astype(np.uint16)

    if progress_callback:
        progress_callback(
            f"{display_name}: ZO-1を読み込み中"
        )

    dapi = read_tiff_gray(dapi_path)
    zo1 = read_tiff_gray(zo1_path)

    if labels.shape != dapi.shape:
        raise ValueError(
            "保存済み核labelとDAPI画像サイズが一致しません。\n"
            f"labels={labels.shape}, DAPI={dapi.shape}"
        )

    if zo1.shape != dapi.shape:
        raise ValueError(
            "DAPIとZO-1画像サイズが一致しません。\n"
            f"DAPI={dapi.shape}, ZO-1={zo1.shape}"
        )

    paths["dir"].mkdir(
        parents=True,
        exist_ok=True
    )

    n_lanes = max(
        4,
        int(params.get(
            "junction_lanes",
            DEFAULT_JUNCTION_LANES
        ))
    )

    pair_factor = max(
        1.0,
        float(params.get(
            "pair_distance_factor",
            DEFAULT_PAIR_DISTANCE_FACTOR
        ))
    )

    max_other = float(
        np.clip(
            params.get(
                "max_other_nucleus_fraction",
                DEFAULT_MAX_OTHER_NUCLEUS_FRACTION
            ),
            0.0,
            1.0
        )
    )

    lane_min_skel = max(
        1,
        int(params.get(
            "lane_min_skeleton_pixels",
            DEFAULT_LANE_MIN_SKELETON_PIXELS
        ))
    )

    lane_peak_top_n = max(
        1,
        int(params.get(
            "lane_peak_top_n_pixels",
            DEFAULT_LANE_PEAK_TOP_N_PIXELS
        ))
    )

    intact_cutoff = float(
        np.clip(
            params.get(
                "intact_coverage_threshold_pct",
                DEFAULT_INTACT_COVERAGE_THRESHOLD_PCT
            ),
            0,
            100
        )
    )

    if progress_callback:
        progress_callback(
            f"{display_name}: ZO-1二値化 / skeleton"
        )

    (
        zo1_8,
        zo1_bg_sub,
        zo1_bg_sub8,
        zo1_tophat,
        zo1_tophat8,
        zo1_binary,
        zo1_skel,
        zo1_thr,
    ) = _preprocess_zo1_for_junction(
        zo1,
        rolling_ball_radius_px=int(
            params.get(
                "zo1_rolling_ball_radius_px",
                DEFAULT_ZO1_ROLLING_BALL_RADIUS_PX
            )
        ),
        tophat_radius_px=int(
            params.get(
                "zo1_tophat_radius_px",
                DEFAULT_ZO1_TOPHAT_RADIUS_PX
            )
        ),
        tophat_threshold=float(
            params.get(
                "zo1_tophat_threshold",
                DEFAULT_ZO1_TOPHAT_THRESHOLD
            )
        ),
        min_object_size=int(
            params.get(
                "zo1_min_object_size_px",
                DEFAULT_ZO1_MIN_OBJECT_SIZE_PX
            )
        ),
    )

    tifffile.imwrite(
        str(paths["dir"] / ZO1_BG_SUB_FILENAME),
        zo1_bg_sub.astype(np.float32)
    )

    tifffile.imwrite(
        str(paths["dir"] / ZO1_TOPHAT_FILENAME),
        zo1_tophat.astype(np.float32)
    )

    tifffile.imwrite(
        str(paths["dir"] / ZO1_BINARY_FILENAME),
        zo1_binary.astype(np.uint8) * 255
    )

    tifffile.imwrite(
        str(paths["dir"] / ZO1_SKELETON_FILENAME),
        zo1_skel.astype(np.uint8) * 255
    )

    if progress_callback:
        progress_callback(
            f"{display_name}: 核ペアを作成中"
        )

    regions = _nucleus_regions_from_labels(
        labels
    )

    label_ids = sorted(regions.keys())

    if len(label_ids) < 2:
        raise RuntimeError(
            "有効な核が2個未満です。"
        )

    centroids = np.array(
        [
            regions[lab]["centroid_xy"]
            for lab in label_ids
        ],
        dtype=np.float64
    )

    max_pair_dist = _adaptive_pair_distance(
        centroids,
        pair_factor
    )

    # All pairs within the adaptive distance are candidates.
    # This is deliberately broader than Delaunay adjacency.
    edges = _candidate_pairs_within_distance(
        centroids,
        max_pair_dist
    )

    pair_rows = []
    draw_items = []

    reject_counts = {
        "too_far": 0,
        "invalid_corridor": 0,
        "third_nucleus": 0,
        "small_corridor": 0,
    }

    edge_list = sorted(edges)

    for edge_idx, (ia, ib) in enumerate(
        edge_list,
        start=1
    ):
        if progress_callback and (
            edge_idx == 1
            or edge_idx % 50 == 0
            or edge_idx == len(edge_list)
        ):
            progress_callback(
                f"{display_name}: 核ペア解析 "
                f"{edge_idx}/{len(edge_list)}"
            )

        lab_a = label_ids[ia]
        lab_b = label_ids[ib]

        reg_a = regions[lab_a]
        reg_b = regions[lab_b]

        d = float(
            np.linalg.norm(
                reg_b["centroid_xy"]
                - reg_a["centroid_xy"]
            )
        )

        if np.isfinite(max_pair_dist) and d > max_pair_dist:
            reject_counts["too_far"] += 1
            continue

        corridor = _build_pair_corridor(
            reg_a,
            reg_b,
            labels.shape
        )

        if corridor is None:
            reject_counts["invalid_corridor"] += 1
            continue

        minx, miny, maxx, maxy = corridor["bbox"]
        cmask = corridor["mask_local"]

        if int(np.sum(cmask)) < n_lanes:
            reject_counts["small_corridor"] += 1
            continue

        lab_crop = labels[
            miny:maxy + 1,
            minx:maxx + 1
        ]

        other = (
            cmask
            & (lab_crop != 0)
            & (lab_crop != lab_a)
            & (lab_crop != lab_b)
        )

        other_pixels = int(
            np.sum(other)
        )

        other_fraction = float(
            other_pixels
            / max(1, np.sum(cmask))
        )

        # Strict line-of-sight rule:
        # if ANY third nucleus enters the center-to-center corridor,
        # do not measure this pair.
        if other_pixels > 0:
            reject_counts["third_nucleus"] += 1
            continue

        t = corridor["lane_fraction_local"]

        zsk = zo1_skel[
            miny:maxy + 1,
            minx:maxx + 1
        ]

        zbin = zo1_binary[
            miny:maxy + 1,
            minx:maxx + 1
        ]

        z8crop = zo1_8[
            miny:maxy + 1,
            minx:maxx + 1
        ]

        # Intensity quantification uses rolling-ball background-subtracted
        # ZO-1, NOT white top-hat and NOT only skeleton-positive pixels.
        zbgcrop = zo1_bg_sub[
            miny:maxy + 1,
            minx:maxx + 1
        ]

        lane_positive = []
        lane_skel_counts = []
        lane_peak_intensities = []

        for lane in range(n_lanes):
            lo = lane / n_lanes
            hi = (lane + 1) / n_lanes

            if lane == n_lanes - 1:
                lm = (
                    cmask
                    & (t >= lo)
                    & (t <= hi)
                )
            else:
                lm = (
                    cmask
                    & (t >= lo)
                    & (t < hi)
                )

            n_skel = int(
                np.sum(zsk & lm)
            )

            lane_skel_counts.append(
                n_skel
            )

            lane_positive.append(
                n_skel >= lane_min_skel
            )

            # Robust lane peak intensity:
            # mean of brightest N pixels in this geometrical lane.
            # Every valid lane contributes even when skeleton/binary is absent.
            lane_vals = zbgcrop[lm]

            if lane_vals.size > 0:
                n_top = min(
                    lane_peak_top_n,
                    int(lane_vals.size)
                )

                if n_top == int(lane_vals.size):
                    lane_peak = float(
                        np.mean(lane_vals)
                    )
                else:
                    top_vals = np.partition(
                        lane_vals,
                        lane_vals.size - n_top
                    )[-n_top:]

                    lane_peak = float(
                        np.mean(top_vals)
                    )
            else:
                lane_peak = np.nan

            lane_peak_intensities.append(
                lane_peak
            )

        positive_lanes = int(
            np.sum(lane_positive)
        )

        coverage = (
            100.0
            * positive_lanes
            / n_lanes
        )

        lane_peak_arr = np.asarray(
            lane_peak_intensities,
            dtype=np.float64
        )

        valid_lane_peaks = lane_peak_arr[
            np.isfinite(lane_peak_arr)
        ]

        pair_mean_lane_peak = (
            float(np.mean(valid_lane_peaks))
            if valid_lane_peaks.size
            else np.nan
        )

        pair_median_lane_peak = (
            float(np.median(valid_lane_peaks))
            if valid_lane_peaks.size
            else np.nan
        )

        zvals = z8crop[cmask]

        zo1_mean = (
            float(np.mean(zvals))
            if zvals.size
            else np.nan
        )

        positive_area_fraction = float(
            np.sum(zbin & cmask)
            / max(1, np.sum(cmask))
        )

        row = {
            "sample_folder": display_name,
            "sample_name": display_name,
            "group": sample_group,
            "pair_id": len(pair_rows) + 1,
            "nucleus_A": int(lab_a),
            "nucleus_B": int(lab_b),
            "centroid_A_x": float(
                reg_a["centroid_xy"][0]
            ),
            "centroid_A_y": float(
                reg_a["centroid_xy"][1]
            ),
            "centroid_B_x": float(
                reg_b["centroid_xy"][0]
            ),
            "centroid_B_y": float(
                reg_b["centroid_xy"][1]
            ),
            "nuclear_distance_px": d,
            "corridor_length_px": float(
                corridor["x_end"]
                - corridor["x_start"]
            ),
            "corridor_area_px": int(
                np.sum(cmask)
            ),
            "other_nucleus_fraction": other_fraction,
            "n_lanes": n_lanes,
            "positive_lanes": positive_lanes,
            "junction_coverage_pct": coverage,
            "intact_ge_cutoff": bool(
                coverage >= intact_cutoff
            ),

            # Junctional intensity independent of detection threshold:
            # rolling-ball BG-subtracted ZO-1, top-N pixel mean per lane,
            # then averaged across all geometrical lanes.
            "lane_peak_top_n_pixels": int(
                lane_peak_top_n
            ),
            "mean_lane_peak_intensity": pair_mean_lane_peak,
            "median_lane_peak_intensity": pair_median_lane_peak,
            "lane_peak_intensities": ";".join(
                (
                    f"{x:.6g}"
                    if np.isfinite(x)
                    else "nan"
                )
                for x in lane_peak_intensities
            ),

            "zo1_corridor_mean_8bit": zo1_mean,
            "zo1_positive_area_fraction": positive_area_fraction,
            "lane_positive_pattern": "".join(
                "1" if x else "0"
                for x in lane_positive
            ),
            "lane_skeleton_pixels": ";".join(
                str(x)
                for x in lane_skel_counts
            ),
        }

        pair_rows.append(row)

        draw_items.append(
            (
                row,
                corridor
            )
        )

    _write_csv(
        paths["dir"] / PAIR_CSV_FILENAME,
        pair_rows,
        fieldnames=[
            "sample_folder",
            "sample_name",
            "group",
            "pair_id",
            "nucleus_A",
            "nucleus_B",
            "centroid_A_x",
            "centroid_A_y",
            "centroid_B_x",
            "centroid_B_y",
            "nuclear_distance_px",
            "corridor_length_px",
            "corridor_area_px",
            "other_nucleus_fraction",
            "n_lanes",
            "positive_lanes",
            "junction_coverage_pct",
            "intact_ge_cutoff",
            "lane_peak_top_n_pixels",
            "mean_lane_peak_intensity",
            "median_lane_peak_intensity",
            "lane_peak_intensities",
            "zo1_corridor_mean_8bit",
            "zo1_positive_area_fraction",
            "lane_positive_pattern",
            "lane_skeleton_pixels",
        ]
    )

    coverages = np.array(
        [
            float(r["junction_coverage_pct"])
            for r in pair_rows
        ],
        dtype=np.float64
    )

    pair_lane_peak_intensities = np.array(
        [
            float(r["mean_lane_peak_intensity"])
            for r in pair_rows
            if np.isfinite(
                float(r["mean_lane_peak_intensity"])
            )
        ],
        dtype=np.float64
    )

    summary = {
        "sample_folder": display_name,
        "sample_name": display_name,
        "group": sample_group,
        "dapi_file": dapi_path.name,
        "zo1_file": zo1_path.name,
        "n_reviewed_nuclei": int(len(label_ids)),
        "n_distance_candidate_pairs": int(len(edges)),
        "n_valid_pairs": int(len(pair_rows)),
        "pair_distance_factor": pair_factor,
        "pair_candidate_rule": "all_centroid_pairs_within_adaptive_distance",
        "third_nucleus_rule": "exclude_if_any_third_nucleus_pixel_in_corridor",
        "max_pair_distance_px": (
            float(max_pair_dist)
            if np.isfinite(max_pair_dist)
            else ""
        ),
        "junction_lanes": n_lanes,
        "lane_peak_top_n_pixels": int(
            lane_peak_top_n
        ),
        "mean_lane_peak_intensity": (
            float(
                np.mean(
                    pair_lane_peak_intensities
                )
            )
            if pair_lane_peak_intensities.size
            else np.nan
        ),
        "median_pair_lane_peak_intensity": (
            float(
                np.median(
                    pair_lane_peak_intensities
                )
            )
            if pair_lane_peak_intensities.size
            else np.nan
        ),
        "sd_pair_lane_peak_intensity": (
            float(
                np.std(
                    pair_lane_peak_intensities,
                    ddof=1
                )
            )
            if pair_lane_peak_intensities.size > 1
            else np.nan
        ),
        "zo1_rolling_ball_radius_px": int(
            params.get(
                "zo1_rolling_ball_radius_px",
                DEFAULT_ZO1_ROLLING_BALL_RADIUS_PX
            )
        ),
        "zo1_tophat_radius_px": int(
            params.get(
                "zo1_tophat_radius_px",
                DEFAULT_ZO1_TOPHAT_RADIUS_PX
            )
        ),
        "zo1_tophat_threshold": float(zo1_thr),
        "corridor_longitudinal_definition": "nucleus-center_to_nucleus-center",
        "mean_junction_coverage_pct": (
            float(np.mean(coverages))
            if coverages.size
            else np.nan
        ),
        "median_junction_coverage_pct": (
            float(np.median(coverages))
            if coverages.size
            else np.nan
        ),
        "sd_junction_coverage_pct": (
            float(np.std(coverages, ddof=1))
            if coverages.size > 1
            else np.nan
        ),
        "intact_coverage_cutoff_pct": intact_cutoff,
        "intact_pair_fraction": (
            float(
                np.mean(
                    coverages >= intact_cutoff
                )
            )
            if coverages.size
            else np.nan
        ),
        "rejected_too_far": int(
            reject_counts["too_far"]
        ),
        "rejected_invalid_corridor": int(
            reject_counts["invalid_corridor"]
        ),
        "rejected_third_nucleus": int(
            reject_counts["third_nucleus"]
        ),
        "rejected_small_corridor": int(
            reject_counts["small_corridor"]
        ),
    }

    _write_csv(
        paths["dir"] / JUNCTION_SUMMARY_FILENAME,
        [summary]
    )

    if progress_callback:
        progress_callback(
            f"{display_name}: overlay作成中"
        )

    # --------------------------------------------------------
    # Requested QC overlay 1:
    # RAW ZO-1 grayscale + ZO-1 skeleton
    #
    # raw ZO-1 = grayscale
    # skeleton = red
    # --------------------------------------------------------
    raw_skel_overlay = cv2.cvtColor(
        zo1_8,
        cv2.COLOR_GRAY2BGR
    )

    # thicken only for visibility in the QC PNG;
    # the ANALYSIS skeleton itself remains one-pixel wide.
    skel_vis = cv2.dilate(
        zo1_skel.astype(np.uint8),
        np.ones((3, 3), dtype=np.uint8),
        iterations=1
    ).astype(bool)

    raw_skel_overlay[skel_vis] = (
        0,
        0,
        255
    )

    imwrite_unicode(
        str(
            paths["dir"]
            / ZO1_RAW_SKELETON_OVERLAY_FILENAME
        ),
        raw_skel_overlay
    )

    # White top-hat QC:
    # top-hat = grayscale, analysis skeleton = red
    tophat_skel_overlay = cv2.cvtColor(
        zo1_tophat8,
        cv2.COLOR_GRAY2BGR
    )

    tophat_skel_overlay[skel_vis] = (
        0,
        0,
        255
    )

    imwrite_unicode(
        str(
            paths["dir"]
            / ZO1_TOPHAT_SKELETON_OVERLAY_FILENAME
        ),
        tophat_skel_overlay
    )

    # --------------------------------------------------------
    # Requested QC overlay 2:
    # DAPI + reviewed nucleus position + ZO-1 skeleton
    #
    # DAPI       = blue-dominant grayscale
    # nuclei     = green contour + yellow centroid
    # ZO-1 skel  = red
    # --------------------------------------------------------
    d8_qc = normalize_8bit(dapi)

    nuclei_skel_overlay = np.zeros(
        (d8_qc.shape[0], d8_qc.shape[1], 3),
        dtype=np.uint8
    )

    # BGR image: DAPI emphasized in blue.
    nuclei_skel_overlay[..., 0] = d8_qc
    nuclei_skel_overlay[..., 1] = (
        0.25 * d8_qc
    ).astype(np.uint8)
    nuclei_skel_overlay[..., 2] = (
        0.10 * d8_qc
    ).astype(np.uint8)

    nuclei_skel_overlay[skel_vis] = (
        0,
        0,
        255
    )

    for lab, region in regions.items():
        x1, y1, x2, y2 = region["bbox"]

        local = (
            labels[y1:y2, x1:x2] == lab
        ).astype(np.uint8)

        contours, _ = cv2.findContours(
            local,
            cv2.RETR_EXTERNAL,
            cv2.CHAIN_APPROX_SIMPLE
        )

        shifted = []

        for cnt in contours:
            c = cnt.copy()
            c[:, 0, 0] += x1
            c[:, 0, 1] += y1
            shifted.append(c)

        cv2.drawContours(
            nuclei_skel_overlay,
            shifted,
            -1,
            (0, 255, 0),
            1
        )

        cx, cy = region["centroid_xy"]

        cv2.circle(
            nuclei_skel_overlay,
            (
                int(round(cx)),
                int(round(cy))
            ),
            2,
            (0, 255, 255),
            -1,
            cv2.LINE_AA
        )

    imwrite_unicode(
        str(
            paths["dir"]
            / ZO1_NUCLEI_SKELETON_OVERLAY_FILENAME
        ),
        nuclei_skel_overlay
    )

    # Debug overlay:
    # DAPI=blue-ish grayscale, ZO-1=red, reviewed nucleus contours=green,
    # pair trapezoids colored by coverage.
    d8 = normalize_8bit(dapi)

    overlay = np.zeros(
        (d8.shape[0], d8.shape[1], 3),
        dtype=np.uint8
    )

    overlay[..., 0] = d8  # B
    overlay[..., 1] = (
        0.35 * d8
    ).astype(np.uint8)
    overlay[..., 2] = np.maximum(
        (
            0.20 * d8
        ).astype(np.uint8),
        zo1_8
    )

    # nuclei contours
    for lab, region in regions.items():
        x1, y1, x2, y2 = region["bbox"]

        local = (
            labels[y1:y2, x1:x2] == lab
        ).astype(np.uint8)

        contours, _ = cv2.findContours(
            local,
            cv2.RETR_EXTERNAL,
            cv2.CHAIN_APPROX_SIMPLE
        )

        shifted = []

        for cnt in contours:
            c = cnt.copy()
            c[:, 0, 0] += x1
            c[:, 0, 1] += y1
            shifted.append(c)

        cv2.drawContours(
            overlay,
            shifted,
            -1,
            (0, 255, 0),
            1
        )

    # pair trapezoids and internal lane guides
    for row, corridor in draw_items:
        coverage = float(
            row["junction_coverage_pct"]
        )

        color = _coverage_color_bgr(
            coverage
        )

        corners = np.round(
            corridor["corners_xy"]
        ).astype(np.int32).reshape(-1, 1, 2)

        cv2.polylines(
            overlay,
            [corners],
            True,
            color,
            1,
            cv2.LINE_AA
        )

        # Draw 4 internal lane guide lines only (20 lanes remains in analysis)
        c = corridor["corners_xy"]

        # c0->c3 left side, c1->c2 right side
        for frac in (0.2, 0.4, 0.6, 0.8):
            p_left = (
                c[0] * (1 - frac)
                + c[3] * frac
            )
            p_right = (
                c[1] * (1 - frac)
                + c[2] * frac
            )

            cv2.line(
                overlay,
                tuple(np.round(p_left).astype(int)),
                tuple(np.round(p_right).astype(int)),
                color,
                1,
                cv2.LINE_AA
            )

        mid = np.mean(
            corridor["corners_xy"],
            axis=0
        )

        cv2.putText(
            overlay,
            f"{coverage:.0f}",
            tuple(np.round(mid).astype(int)),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.35,
            color,
            1,
            cv2.LINE_AA
        )

    imwrite_unicode(
        str(
            paths["dir"]
            / JUNCTION_OVERLAY_FILENAME
        ),
        overlay
    )

    return summary



def merge_touching_labels_8neighbor(labels):
    """
    Automatically merge distinct labels that directly touch in an 8-neighborhood.

    Used by BATCH StarDist before draft saving.

    Returns:
        merged_labels (uint16, contiguous relabeling)
        n_merged_groups
        n_labels_removed
    """
    L = np.asarray(labels, dtype=np.uint16)

    if L.size == 0:
        return L.copy(), 0, 0

    max_lab = int(L.max())

    if max_lab <= 1:
        return L.copy(), 0, 0

    present = np.unique(L)
    present = present[present > 0].astype(np.int32)

    if present.size <= 1:
        return L.copy(), 0, 0

    parent = np.arange(
        max_lab + 1,
        dtype=np.int32
    )

    def find(a):
        a = int(a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = int(parent[a])
        return a

    def union(a, b):
        ra = find(a)
        rb = find(b)

        if ra == rb:
            return

        if ra < rb:
            parent[rb] = ra
        else:
            parent[ra] = rb

    base = np.uint64(max_lab + 1)

    def process_views(a, b):
        m = (
            (a > 0)
            & (b > 0)
            & (a != b)
        )

        if not np.any(m):
            return

        av = a[m].astype(
            np.uint64,
            copy=False
        )
        bv = b[m].astype(
            np.uint64,
            copy=False
        )

        lo = np.minimum(av, bv)
        hi = np.maximum(av, bv)

        codes = np.unique(
            lo * base + hi
        )

        for code in codes:
            aa = int(code // base)
            bb = int(code % base)

            if (
                aa > 0
                and bb > 0
                and aa != bb
            ):
                union(aa, bb)

    # Horizontal / vertical / both diagonals = 8-neighborhood.
    process_views(
        L[:, :-1],
        L[:, 1:]
    )
    process_views(
        L[:-1, :],
        L[1:, :]
    )
    process_views(
        L[:-1, :-1],
        L[1:, 1:]
    )
    process_views(
        L[1:, :-1],
        L[:-1, 1:]
    )

    roots = np.array(
        [find(lab) for lab in present],
        dtype=np.int32
    )

    unique_roots, root_counts = np.unique(
        roots,
        return_counts=True
    )

    n_groups = int(
        np.sum(root_counts >= 2)
    )

    n_removed = int(
        np.sum(
            root_counts[root_counts >= 2] - 1
        )
    )

    if n_removed <= 0:
        # StarDist filtering already returns contiguous labels.
        return L.copy(), 0, 0

    root_map = np.arange(
        max_lab + 1,
        dtype=np.int32
    )

    for lab in present:
        root_map[int(lab)] = find(lab)

    merged_root = root_map[L]

    # Relabel surviving groups contiguously 1..N.
    surviving = np.unique(
        merged_root
    )
    surviving = surviving[
        surviving > 0
    ]

    contiguous = np.zeros(
        max_lab + 1,
        dtype=np.uint16
    )

    for new_lab, root in enumerate(
        surviving,
        start=1
    ):
        contiguous[int(root)] = new_lab

    out = contiguous[
        merged_root
    ].astype(np.uint16)

    return out, n_groups, n_removed


# ============================================================
# StarDist
# ============================================================

class StarDistNucleiSegmenter:
    _model = None

    @classmethod
    def get_model(cls):
        if cls._model is not None:
            return cls._model

        try:
            from stardist.models import StarDist2D
        except Exception as exc:
            raise RuntimeError(
                "StarDistを読み込めません。\n"
                "pip install tensorflow stardist csbdeep\n\n"
                f"{exc}"
            ) from exc

        try:
            cls._model = StarDist2D.from_pretrained("2D_versatile_fluo")
            return cls._model
        except OSError as exc:
            # Windows symlink workaround from the reference code.
            if getattr(exc, "winerror", None) != 1314:
                raise

            root = (
                Path.home()
                / ".keras"
                / "models"
                / "StarDist2D"
                / "2D_versatile_fluo"
            )
            extracted = root / "2D_versatile_fluo_extracted"
            target = root / "2D_versatile_fluo"

            if target.is_symlink():
                target.unlink()

            if extracted.exists() and not target.exists():
                extracted.rename(target)

            cls._model = StarDist2D.from_pretrained("2D_versatile_fluo")
            return cls._model

    @staticmethod
    def percentile_normalize(image, pmin=1.0, pmax=99.8):
        x = np.asarray(image, dtype=np.float32)
        lo, hi = np.percentile(x, [float(pmin), float(pmax)])

        if hi <= lo:
            lo, hi = float(x.min()), float(x.max())

        if hi <= lo:
            return np.zeros_like(x, dtype=np.float32)

        return np.clip((x - lo) / (hi - lo), 0, 1).astype(np.float32)

    @staticmethod
    def _filter_relabel(labels, min_area_px, max_area_px):
        labels = np.asarray(labels, dtype=np.int32)
        counts = np.bincount(labels.ravel())

        remap = np.zeros(len(counts), dtype=np.uint16)
        next_lab = 1

        for lab in range(1, len(counts)):
            area = int(counts[lab])

            if area < min_area_px:
                continue

            if max_area_px > 0 and area > max_area_px:
                continue

            if next_lab >= np.iinfo(np.uint16).max:
                raise RuntimeError("核数がuint16上限を超えました。")

            remap[lab] = next_lab
            next_lab += 1

        return remap[labels].astype(np.uint16)

    @classmethod
    def segment(
        cls,
        dapi_image,
        pixel_size_um,
        prob_thresh=0.5,
        nms_thresh=0.4,
        min_area_um2=15.0,
        max_area_um2=1000.0,
        normalize_pmin=1.0,
        normalize_pmax=99.8,
    ):
        normalized = cls.percentile_normalize(
            dapi_image,
            normalize_pmin,
            normalize_pmax,
        )

        h, w = normalized.shape

        n_tiles = (
            max(1, int(np.ceil(h / 1024))),
            max(1, int(np.ceil(w / 1024))),
        )

        # In a PyInstaller windowed EXE, stdout/stderr can be None.
        # StarDist/tqdm requires a writable stream during tiled prediction.
        ensure_writable_standard_streams()

        labels, _ = cls.get_model().predict_instances(
            normalized,
            prob_thresh=float(prob_thresh),
            nms_thresh=float(nms_thresh),
            n_tiles=n_tiles,
        )

        px_area = max(1e-12, float(pixel_size_um) ** 2)

        min_area_px = max(
            1,
            int(round(float(min_area_um2) / px_area))
        )

        max_area_px = (
            max(0, int(round(float(max_area_um2) / px_area)))
            if float(max_area_um2) > 0
            else 0
        )

        return cls._filter_relabel(
            labels,
            min_area_px,
            max_area_px,
        )


# ============================================================
# Nuclei manager
# ============================================================

class NucleiLabelManager:
    """
    DisplayManager互換のinstance label manager。
    元コードの操作系を維持しつつ、brush自体は局所領域だけ作る。
    """

    def __init__(self, shape, ann_path=None):
        self.shape = tuple(shape)
        self.label_mask = np.zeros(self.shape, dtype=np.uint16)

        self.next_label_id = 1
        self.selected_labels = set()

        # DisplayManager compatibility
        self.roi_done = set()
        self.centerlines = {}

        self.current_pen_label = None
        self.eraser_touched_labels = set()

        self.undo_stack = []
        self.max_undo = 20

        self.bbox_cache = {}
        self.dirty_labels = set()

        self.metadata = {}

        if ann_path is not None and Path(ann_path).exists():
            self.load(ann_path)

    # ---------------- load/save ----------------

    def load(self, path):
        d = np.load(path, allow_pickle=True)

        if "nuclei_labels" not in d:
            raise RuntimeError(f"nuclei_labelsがありません: {path}")

        labels = d["nuclei_labels"].astype(np.uint16)

        if labels.shape != self.shape:
            raise ValueError(
                "保存済み核ラベルとDAPIサイズが一致しません。\n"
                f"保存: {labels.shape}\n画像: {self.shape}"
            )

        self.label_mask = labels
        self.next_label_id = int(labels.max()) + 1
        self.selected_labels.clear()

        for key in (
            "sample_folder",
            "dapi_filename",
            "zo1_filename",
            "pixel_size_um",
            "prob_thresh",
            "nms_thresh",
            "min_area_um2",
            "max_area_um2",
            "normalize_pmin",
            "normalize_pmax",
        ):
            if key in d:
                v = d[key]
                try:
                    v = v.item()
                except Exception:
                    pass
                self.metadata[key] = v

        self.rebuild_all_bboxes()
        self.mark_all_dirty()

    def save(self, path, sample_name, dapi_name, zo1_name, params):
        np.savez_compressed(
            path,
            nuclei_labels=self.label_mask.astype(np.uint16),
            nuclei_binary=(self.label_mask > 0).astype(np.uint8) * 255,
            sample_folder=str(sample_name),
            dapi_filename=str(dapi_name),
            zo1_filename=str(zo1_name),
            pixel_size_um=float(params["pixel_size_um"]),
            prob_thresh=float(params["prob_thresh"]),
            nms_thresh=float(params["nms_thresh"]),
            min_area_um2=float(params["min_area_um2"]),
            max_area_um2=float(params["max_area_um2"]),
            normalize_pmin=float(params["normalize_pmin"]),
            normalize_pmax=float(params["normalize_pmax"]),
        )

        self.metadata = dict(params)

    # ---------------- labels/meta ----------------

    def count(self):
        return len(self.bbox_cache)

    def new_label(self):
        lab = int(self.next_label_id)

        if lab >= np.iinfo(np.uint16).max:
            raise RuntimeError("核ラベル数がuint16上限に達しました。")

        self.next_label_id += 1
        return lab

    def mark_dirty(self, lab):
        lab = int(lab)
        if lab > 0:
            self.dirty_labels.add(lab)

    def mark_all_dirty(self):
        self.dirty_labels.update(self.bbox_cache.keys())

    def consume_dirty_labels(self):
        out = set(self.dirty_labels)
        self.dirty_labels.clear()
        return out

    def rebuild_all_bboxes(self):
        self.bbox_cache.clear()

        if int(self.label_mask.max()) == 0:
            return

        # scipy.ndimage.find_objects is much faster than np.where per label.
        if find_objects is not None:
            objects = find_objects(self.label_mask)

            for i, slc in enumerate(objects, start=1):
                if slc is None:
                    continue

                ys, xs = slc
                self.bbox_cache[i] = (
                    int(xs.start),
                    int(ys.start),
                    int(xs.stop),
                    int(ys.stop),
                )
            return

        # Fallback
        for lab in np.unique(self.label_mask):
            lab = int(lab)
            if lab > 0:
                self.update_bbox(lab)

    def update_bbox(self, lab):
        lab = int(lab)

        ys, xs = np.where(self.label_mask == lab)

        if len(xs) == 0:
            self.bbox_cache.pop(lab, None)
            return

        self.bbox_cache[lab] = (
            int(xs.min()),
            int(ys.min()),
            int(xs.max()) + 1,
            int(ys.max()) + 1,
        )

    def update_dirty_bboxes(self):
        dirty = list(self.dirty_labels)
        for lab in dirty:
            self.update_bbox(lab)

    def set_label_mask(self, labels):
        labels = np.asarray(labels, dtype=np.uint16)

        if labels.shape != self.shape:
            raise ValueError("核ラベル画像サイズが一致しません。")

        self.label_mask = labels.copy()
        self.next_label_id = int(labels.max()) + 1
        self.selected_labels.clear()

        self.rebuild_all_bboxes()
        self.mark_all_dirty()

    # ---------------- undo ----------------

    def push_undo(self):
        self.undo_stack.append({
            "label_mask": self.label_mask.copy(),
            "selected_labels": set(self.selected_labels),
            "next_label_id": int(self.next_label_id),
        })

        if len(self.undo_stack) > self.max_undo:
            self.undo_stack.pop(0)

    def undo(self):
        if not self.undo_stack:
            return False

        state = self.undo_stack.pop()

        self.label_mask = state["label_mask"]
        self.selected_labels = state["selected_labels"]
        self.next_label_id = state["next_label_id"]

        self.rebuild_all_bboxes()
        self.mark_all_dirty()

        return True

    # ---------------- brush ----------------

    def draw_brush(self, ix, iy, screen_radius, zoom, add=True):
        h, w = self.shape

        radius = max(
            1,
            int(round(float(screen_radius) / max(float(zoom), 1e-6)))
        )

        x0 = max(0, int(ix) - radius - 1)
        x1 = min(w, int(ix) + radius + 2)
        y0 = max(0, int(iy) - radius - 1)
        y1 = min(h, int(iy) + radius + 2)

        if x1 <= x0 or y1 <= y0:
            return

        brush = np.zeros((y1 - y0, x1 - x0), dtype=np.uint8)

        cv2.circle(
            brush,
            (int(ix) - x0, int(iy) - y0),
            radius,
            1,
            -1,
        )

        local = self.label_mask[y0:y1, x0:x1]

        if add:
            if self.current_pen_label is None:
                if len(self.selected_labels) == 1:
                    self.current_pen_label = int(
                        next(iter(self.selected_labels))
                    )
                else:
                    self.current_pen_label = self.new_label()

            lab = int(self.current_pen_label)

            writable = (
                (brush > 0)
                & ((local == 0) | (local == lab))
            )

            local[writable] = lab
            self.mark_dirty(lab)

        else:
            touched = np.unique(local[brush > 0])
            touched = touched[touched > 0]

            local[brush > 0] = 0

            for lab in touched:
                lab = int(lab)
                self.eraser_touched_labels.add(lab)
                self.mark_dirty(lab)

    def finish_pen(self):
        if self.current_pen_label is not None:
            self.mark_dirty(self.current_pen_label)

        self.current_pen_label = None
        self.update_dirty_bboxes()

    def _split_label_components(self, lab):
        lab = int(lab)

        binary = (self.label_mask == lab).astype(np.uint8)

        n, comp = cv2.connectedComponents(binary, 8)

        if n <= 2:
            if int(binary.sum()) == 0:
                self.selected_labels.discard(lab)
                self.bbox_cache.pop(lab, None)
            else:
                self.update_bbox(lab)
            return

        self.label_mask[self.label_mask == lab] = 0

        first = True

        for cid in range(1, n):
            new_lab = lab if first else self.new_label()
            first = False

            self.label_mask[comp == cid] = new_lab
            self.mark_dirty(new_lab)
            self.update_bbox(new_lab)

    def finish_eraser(self):
        for lab in list(self.eraser_touched_labels):
            self._split_label_components(lab)

        self.eraser_touched_labels.clear()
        self.update_dirty_bboxes()

    # ---------------- commands ----------------

    def delete_selected(self):
        if not self.selected_labels:
            return False

        self.push_undo()

        for lab in list(self.selected_labels):
            lab = int(lab)
            self.label_mask[self.label_mask == lab] = 0
            self.bbox_cache.pop(lab, None)
            self.mark_dirty(lab)

        self.selected_labels.clear()
        return True

    def merge_selected(self):
        if len(self.selected_labels) < 2:
            return False

        self.push_undo()

        labs = sorted(int(x) for x in self.selected_labels)
        target = labs[0]

        for lab in labs[1:]:
            self.label_mask[self.label_mask == lab] = target
            self.bbox_cache.pop(lab, None)
            self.mark_dirty(lab)

        self.selected_labels = {target}

        self.mark_dirty(target)
        self.update_bbox(target)

        return True

    def fill_holes_selected(self):
        """
        Fill internal holes of selected nucleus labels.

        Other nucleus labels are never overwritten. Only background pixels
        completely enclosed by the selected nucleus are added.
        """
        if not self.selected_labels:
            return False, 0

        if binary_fill_holes is None:
            raise RuntimeError(
                "Fill holeには scipy が必要です。\n"
                "pip install scipy"
            )

        self.push_undo()

        total_added = 0

        for lab in sorted(int(x) for x in self.selected_labels):
            bbox = self.bbox_cache.get(lab)
            if bbox is None:
                continue

            x1, y1, x2, y2 = bbox

            # One-pixel pad gives binary_fill_holes a clear outside background.
            px1 = max(0, x1 - 1)
            py1 = max(0, y1 - 1)
            px2 = min(self.shape[1], x2 + 1)
            py2 = min(self.shape[0], y2 + 1)

            local_labels = self.label_mask[py1:py2, px1:px2]
            binary = (local_labels == lab)

            filled = binary_fill_holes(binary)
            add = filled & (~binary) & (local_labels == 0)

            n_add = int(np.sum(add))
            if n_add <= 0:
                continue

            local_labels[add] = lab
            total_added += n_add

            self.mark_dirty(lab)
            self.update_bbox(lab)

        if total_added <= 0:
            # Do not keep a useless undo step.
            if self.undo_stack:
                self.undo_stack.pop()
            return False, 0

        return True, total_added

    def find_touching_groups(self):
        """
        Return groups of distinct nucleus labels that directly touch
        in an 8-neighborhood.

        No distance threshold or dilation is used:
        only labels sharing an edge or corner are considered touching.
        """
        L = self.label_mask

        if not self.bbox_cache or int(L.max()) <= 1:
            return []

        max_lab = int(L.max())
        parent = np.arange(max_lab + 1, dtype=np.int32)

        def find(a):
            a = int(a)
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = int(parent[a])
            return a

        def union(a, b):
            ra = find(a)
            rb = find(b)
            if ra == rb:
                return
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

        base = np.uint64(max_lab + 1)

        def process_pair_views(a, b):
            mask = (a > 0) & (b > 0) & (a != b)
            if not np.any(mask):
                return

            av = a[mask].astype(np.uint64, copy=False)
            bv = b[mask].astype(np.uint64, copy=False)

            lo = np.minimum(av, bv)
            hi = np.maximum(av, bv)

            codes = np.unique(lo * base + hi)

            for code in codes:
                aa = int(code // base)
                bb = int(code % base)
                if aa > 0 and bb > 0 and aa != bb:
                    union(aa, bb)

        # 8-neighborhood: horizontal, vertical, both diagonals.
        process_pair_views(L[:, :-1], L[:, 1:])
        process_pair_views(L[:-1, :], L[1:, :])
        process_pair_views(L[:-1, :-1], L[1:, 1:])
        process_pair_views(L[1:, :-1], L[:-1, 1:])

        groups = {}

        for lab in self.bbox_cache.keys():
            root = find(lab)
            groups.setdefault(root, []).append(int(lab))

        return [
            sorted(g)
            for g in groups.values()
            if len(g) >= 2
        ]

    def merge_touching_groups(self, groups):
        """
        Merge touching-label groups in one full-image remap pass.
        Returns (changed, n_groups, n_labels_removed).
        """
        groups = [
            sorted(set(int(x) for x in g if int(x) > 0))
            for g in groups
        ]
        groups = [g for g in groups if len(g) >= 2]

        if not groups:
            return False, 0, 0

        self.push_undo()

        max_lab = int(self.label_mask.max())
        remap = np.arange(max_lab + 1, dtype=np.uint16)

        n_removed = 0

        for group in groups:
            target = min(group)

            for lab in group:
                if lab <= max_lab:
                    remap[lab] = target

            n_removed += len(group) - 1

        self.label_mask = remap[self.label_mask]

        self.selected_labels.clear()
        self.next_label_id = int(self.label_mask.max()) + 1

        self.rebuild_all_bboxes()
        self.dirty_labels.clear()
        self.mark_all_dirty()

        return True, len(groups), n_removed

    def split_selected_disconnected(self):
        if len(self.selected_labels) != 1:
            return False, "核を1つだけ選択してください。"

        lab = int(next(iter(self.selected_labels)))
        binary = (self.label_mask == lab).astype(np.uint8)

        n, _ = cv2.connectedComponents(binary, 8)

        if n <= 2:
            return (
                False,
                "選択核は非連結成分に分かれていません。\n"
                "先に消しゴムで分離してください。"
            )

        self.push_undo()
        self._split_label_components(lab)
        self.selected_labels.clear()

        return True, ""

    def clear_all(self):
        if not self.bbox_cache:
            return False

        self.push_undo()

        self.label_mask.fill(0)
        self.selected_labels.clear()
        self.next_label_id = 1
        self.bbox_cache.clear()
        self.dirty_labels.clear()

        return True

    def get_status_text(self, brush_size):
        selected = ", ".join(
            str(x) for x in sorted(self.selected_labels)
        )

        return (
            f"核数: {self.count()}\n"
            f"選択: {selected if selected else 'なし'}\n"
            f"Brush: {int(brush_size)} screen px"
        )


# ============================================================
# Display manager
# Viewport drawing from the user's reference code
# ============================================================

class DisplayManager:
    def __init__(self, canvas_owner, gray8, roi_manager):
        self.owner = canvas_owner
        self.gray8 = gray8
        self.roi = roi_manager

        self.h, self.w = gray8.shape

        self.tk_img = None

        self.base_cache = None
        self.base_cache_params = None

        self.roi_meta_cache = {}
        self.roi_meta_dirty = True

        # ZO-1 is loaded lazily by the editor only when the user asks for it.
        # It is an 8-bit display image and is cropped to the current viewport.
        self.zo1_gray8 = None

    def set_zo1_overlay_image(self, gray8):
        if gray8 is None:
            self.zo1_gray8 = None
            self.tk_img = None
            return

        gray8 = np.asarray(gray8, dtype=np.uint8)

        if gray8.shape != (self.h, self.w):
            raise ValueError(
                "ZO-1 overlayとDAPI画像サイズが一致しません。\n"
                f"DAPI: {(self.h, self.w)}\nZO-1: {gray8.shape}"
            )

        self.zo1_gray8 = gray8
        self.tk_img = None

    # ---------------- invalidate ----------------

    def invalidate_all(self):
        self.tk_img = None
        self.base_cache = None
        self.base_cache_params = None
        self.roi_meta_dirty = True
        self.roi.mark_all_dirty()

    def invalidate_base(self):
        self.tk_img = None
        self.base_cache = None
        self.base_cache_params = None

    def invalidate_roi_layer(self):
        self.tk_img = None
        self.roi_meta_dirty = True

    def invalidate_resize(self):
        self.tk_img = None

    # ---------------- base ----------------

    def make_base_gray(self, display_min, display_max):
        vmin = max(0, int(display_min))
        vmax = max(0, int(display_max))

        if vmax <= vmin:
            vmax = vmin + 1

        params = (vmin, vmax)

        if (
            self.base_cache is not None
            and self.base_cache_params == params
        ):
            return self.base_cache

        img = self.gray8.astype(np.float32)

        out = np.clip(
            (img - vmin) / (vmax - vmin),
            0,
            1,
        )

        out = (out * 255).astype(np.uint8)

        self.base_cache = out
        self.base_cache_params = params

        return out

    # ---------------- viewport ----------------

    def get_viewport_image_rect(
        self,
        canvas,
        zoom,
        offset_x,
        offset_y
    ):
        cw = max(1, canvas.winfo_width())
        ch = max(1, canvas.winfo_height())

        x1 = int(np.floor((0 - offset_x) / zoom))
        y1 = int(np.floor((0 - offset_y) / zoom))
        x2 = int(np.ceil((cw - offset_x) / zoom))
        y2 = int(np.ceil((ch - offset_y) / zoom))

        x1 = max(0, x1)
        y1 = max(0, y1)
        x2 = min(self.w, x2)
        y2 = min(self.h, y2)

        if x2 <= x1 or y2 <= y1:
            return None

        return {
            "ix1": x1,
            "iy1": y1,
            "ix2": x2,
            "iy2": y2,
            "canvas_x": offset_x + x1 * zoom,
            "canvas_y": offset_y + y1 * zoom,
        }

    # ---------------- ROI cache ----------------

    def update_roi_meta(self, lab):
        lab = int(lab)

        bbox = self.roi.bbox_cache.get(lab)

        if bbox is None:
            self.roi_meta_cache.pop(lab, None)
            return

        x1, y1, x2, y2 = bbox

        if x2 <= x1 or y2 <= y1:
            self.roi_meta_cache.pop(lab, None)
            return

        crop = (
            self.roi.label_mask[y1:y2, x1:x2] == lab
        ).astype(np.uint8)

        if int(crop.sum()) == 0:
            self.roi_meta_cache.pop(lab, None)
            return

        contours, _ = cv2.findContours(
            crop,
            cv2.RETR_LIST,
            cv2.CHAIN_APPROX_SIMPLE
        )

        global_contours = []

        for cnt in contours:
            c = cnt.copy()
            c[:, 0, 0] += x1
            c[:, 0, 1] += y1
            global_contours.append(c)

        m = cv2.moments(crop, binaryImage=True)

        if m["m00"] > 0:
            centroid = (
                float(x1 + m["m10"] / m["m00"]),
                float(y1 + m["m01"] / m["m00"]),
            )
        else:
            centroid = (
                float((x1 + x2) / 2),
                float((y1 + y2) / 2),
            )

        self.roi_meta_cache[lab] = {
            "bbox": bbox,
            "contours": global_contours,
            "centroid": centroid,
        }

    def update_dirty_roi_meta(self):
        dirty = self.roi.consume_dirty_labels()

        for lab in dirty:
            self.update_roi_meta(lab)

        # Remove deleted labels without np.unique(full image).
        current = set(self.roi.bbox_cache.keys())

        for lab in list(self.roi_meta_cache.keys()):
            if lab not in current:
                self.roi_meta_cache.pop(lab, None)

        self.roi_meta_dirty = False

    @staticmethod
    def bbox_intersects_view(
        bbox,
        ix1,
        iy1,
        ix2,
        iy2
    ):
        x1, y1, x2, y2 = bbox

        return not (
            x2 < ix1
            or x1 > ix2
            or y2 < iy1
            or y1 > iy2
        )

    def draw_roi_layer_on_view(
        self,
        view_rgb,
        ix1,
        iy1,
        ix2,
        iy2,
        show_contours=True,
        show_numbers=False,
    ):
        self.update_dirty_roi_meta()

        if not show_contours:
            return

        for lab, meta in self.roi_meta_cache.items():
            if not self.bbox_intersects_view(
                meta["bbox"],
                ix1,
                iy1,
                ix2,
                iy2
            ):
                continue

            color = (
                (255, 255, 0)
                if lab in self.roi.selected_labels
                else (0, 255, 0)
            )

            shifted = []

            for cnt in meta["contours"]:
                c = cnt.copy()
                c[:, 0, 0] -= ix1
                c[:, 0, 1] -= iy1
                shifted.append(c)

            cv2.drawContours(
                view_rgb,
                shifted,
                -1,
                color,
                1,
            )

            if show_numbers:
                cx, cy = meta["centroid"]

                cv2.putText(
                    view_rgb,
                    str(lab),
                    (
                        int(round(cx - ix1)),
                        int(round(cy - iy1))
                    ),
                    cv2.FONT_HERSHEY_SIMPLEX,
                    0.40,
                    color,
                    1,
                    cv2.LINE_AA,
                )

    # ---------------- redraw ----------------

    def redraw(
        self,
        canvas,
        zoom,
        offset_x,
        offset_y,
        display_min,
        display_max,
        show_contours=True,
        show_numbers=False,
        show_zo1_overlay=False,
        zo1_overlay_alpha=0.70,
    ):
        if (
            canvas.winfo_width() <= 1
            or canvas.winfo_height() <= 1
        ):
            return

        rect = self.get_viewport_image_rect(
            canvas,
            zoom,
            offset_x,
            offset_y,
        )

        canvas.delete("all")

        if rect is None:
            return

        ix1 = rect["ix1"]
        iy1 = rect["iy1"]
        ix2 = rect["ix2"]
        iy2 = rect["iy2"]

        base = self.make_base_gray(
            display_min,
            display_max,
        )

        view_gray = base[iy1:iy2, ix1:ix2]

        view_rgb = cv2.cvtColor(
            view_gray,
            cv2.COLOR_GRAY2RGB
        )

        # Optional ZO-1 display overlay.
        # Intensity-dependent alpha keeps weak background from tinting the
        # entire image. Red is used so green/yellow nucleus contours remain clear.
        if (
            show_zo1_overlay
            and self.zo1_gray8 is not None
        ):
            zo1_view = self.zo1_gray8[iy1:iy2, ix1:ix2]

            if zo1_view.shape == view_gray.shape:
                z = zo1_view.astype(np.float32) / 255.0
                alpha = np.clip(
                    z * float(zo1_overlay_alpha),
                    0.0,
                    float(zo1_overlay_alpha)
                )

                basef = view_rgb.astype(np.float32)

                basef[..., 0] = (
                    basef[..., 0] * (1.0 - alpha)
                    + 255.0 * alpha
                )
                basef[..., 1] *= (1.0 - alpha)
                basef[..., 2] *= (1.0 - alpha)

                view_rgb = np.clip(
                    basef,
                    0,
                    255
                ).astype(np.uint8)

        self.draw_roi_layer_on_view(
            view_rgb,
            ix1,
            iy1,
            ix2,
            iy2,
            show_contours=show_contours,
            show_numbers=show_numbers,
        )

        disp_w = max(
            1,
            int(view_rgb.shape[1] * zoom)
        )

        disp_h = max(
            1,
            int(view_rgb.shape[0] * zoom)
        )

        resized = cv2.resize(
            view_rgb,
            (disp_w, disp_h),
            interpolation=cv2.INTER_NEAREST
        )

        self.tk_img = ImageTk.PhotoImage(
            Image.fromarray(resized)
        )

        canvas.create_image(
            rect["canvas_x"],
            rect["canvas_y"],
            anchor=tk.NW,
            image=self.tk_img
        )


# ============================================================
# DAPI editor
# ============================================================

class NucleiAnnotationWindow(tk.Toplevel):
    def __init__(
        self,
        master,
        sample,
        global_settings,
        saved_callback,
        next_callback,
    ):
        super().__init__(master)

        self.sample = sample
        self.folder = Path(sample["folder"])
        self.sample_name = sample_display_name(sample)
        self.dapi_path = Path(sample["dapi"])
        self.zo1_path = Path(sample["zo1"])

        self.paths = sample_output_paths(sample)

        self.saved_callback = saved_callback
        self.next_callback = next_callback

        self.title(
            f"DAPI StarDist detection / review - "
            f"{self.sample_name} - {self.dapi_path.name}"
        )

        self.geometry("1400x900")
        self.minsize(950, 650)
        self.panel_width = 270

        # DAPIだけを読む。ZO-1はここでは読まない。
        self.raw = read_tiff_gray(self.dapi_path)
        self.gray8 = normalize_8bit(self.raw)

        self.h, self.w = self.gray8.shape

        # Reviewed annotation always has priority.
        # If it does not exist, load a batch-StarDist draft when available.
        self.loaded_from_draft = False

        if self.paths["npz"].exists():
            initial_annotation_path = self.paths["npz"]
        elif self.paths["draft_npz"].exists():
            initial_annotation_path = self.paths["draft_npz"]
            self.loaded_from_draft = True
        else:
            initial_annotation_path = None

        self.manager = NucleiLabelManager(
            self.gray8.shape,
            initial_annotation_path
        )

        saved = self.manager.metadata
        gs = global_settings or {}

        self.pixel_size = tk.DoubleVar(
            value=float(
                gs.get(
                    "pixel_size_um",
                    saved.get(
                        "pixel_size_um",
                        DEFAULT_PIXEL_SIZE_UM
                    )
                )
            )
        )

        self.prob_thresh = tk.DoubleVar(
            value=float(
                gs.get(
                    "prob_thresh",
                    saved.get(
                        "prob_thresh",
                        DEFAULT_PROB_THRESH
                    )
                )
            )
        )

        self.nms_thresh = tk.DoubleVar(
            value=float(
                gs.get(
                    "nms_thresh",
                    saved.get(
                        "nms_thresh",
                        DEFAULT_NMS_THRESH
                    )
                )
            )
        )

        self.min_area_um2 = tk.DoubleVar(
            value=float(
                gs.get(
                    "min_area_um2",
                    saved.get(
                        "min_area_um2",
                        DEFAULT_MIN_AREA_UM2
                    )
                )
            )
        )

        self.max_area_um2 = tk.DoubleVar(
            value=float(
                gs.get(
                    "max_area_um2",
                    saved.get(
                        "max_area_um2",
                        DEFAULT_MAX_AREA_UM2
                    )
                )
            )
        )

        self.normalize_pmin = tk.DoubleVar(
            value=float(
                gs.get(
                    "normalize_pmin",
                    saved.get(
                        "normalize_pmin",
                        DEFAULT_NORMALIZE_PMIN
                    )
                )
            )
        )

        self.normalize_pmax = tk.DoubleVar(
            value=float(
                gs.get(
                    "normalize_pmax",
                    saved.get(
                        "normalize_pmax",
                        DEFAULT_NORMALIZE_PMAX
                    )
                )
            )
        )

        self.display_min = tk.IntVar(value=0)
        self.display_max = tk.IntVar(value=255)

        self.show_contours = tk.BooleanVar(value=True)
        self.show_ids = tk.BooleanVar(value=False)

        self.show_zo1_overlay = False
        self.zo1_overlay_alpha = 0.70
        self.zo1_overlay_button_text = tk.StringVar(
            value="ZO-1 overlay: OFF"
        )

        self.auto_display_mode = 0

        # Same annotation navigation style as reference code.
        self.zoom = 1.0
        self.offset_x = 0.0
        self.offset_y = 0.0

        self.last_drag = None
        self.drag_mode = None

        self.brush_screen_radius = 10

        self.left_press_after_id = None
        self.left_press_event = None
        self.left_pressed = False
        self.left_started = False
        self.double_click_delay_ms = 180

        self.processing_queue = queue.Queue()
        self.processing_running = False

        # Persist current StarDist parameters globally, independent of dataset.
        self._detection_settings_after_id = None

        self.unsaved_changes = False

        self.display = DisplayManager(
            self,
            self.gray8,
            self.manager
        )

        self.build_ui()

        for _var in (
            self.pixel_size,
            self.prob_thresh,
            self.nms_thresh,
            self.min_area_um2,
            self.max_area_um2,
            self.normalize_pmin,
            self.normalize_pmax,
        ):
            _var.trace_add(
                "write",
                self._schedule_detection_settings_save
            )

        self.protocol(
            "WM_DELETE_WINDOW",
            self.on_window_close
        )

        self.after(
            60,
            lambda: self.redraw(False)
        )

        if self.paths["npz"].exists():
            self.detect_status.config(
                text=f"Reviewed/Saved: {self.manager.count()} nuclei",
                bg="green"
            )
        elif self.loaded_from_draft:
            self.detect_status.config(
                text=(
                    f"Batch StarDist loaded: "
                    f"{self.manager.count()} nuclei / 未確認"
                ),
                bg="dark orange"
            )
        else:
            # Initial StarDist once only when neither reviewed annotation nor
            # batch draft exists.
            self.after(
                300,
                lambda: self.run_stardist(
                    confirm_replace=False
                )
            )

    def _widget_is_inside_panel(self, widget):
        w = widget

        while w is not None:
            if w == getattr(
                self,
                "panel_container",
                None
            ):
                return True

            try:
                w = w.master
            except Exception:
                break

        return False

    def _on_annotation_window_wheel(self, event):
        """
        Scroll the RIGHT control panel only when the mouse is over it.
        The left image canvas continues to use the wheel for zoom.
        """
        try:
            widget = self.winfo_containing(
                event.x_root,
                event.y_root
            )
        except Exception:
            return

        if not self._widget_is_inside_panel(
            widget
        ):
            return

        if getattr(event, "num", None) == 4:
            steps = -3
        elif getattr(event, "num", None) == 5:
            steps = 3
        else:
            delta = int(
                getattr(
                    event,
                    "delta",
                    0
                )
            )

            if delta == 0:
                return

            steps = (
                -3
                if delta > 0
                else 3
            )

        self.panel_canvas.yview_scroll(
            steps,
            "units"
        )

        return "break"

    def _schedule_detection_settings_save(
        self,
        *_args
    ):
        if self._detection_settings_after_id is not None:
            try:
                self.after_cancel(
                    self._detection_settings_after_id
                )
            except Exception:
                pass

        self._detection_settings_after_id = self.after(
            500,
            self._persist_detection_settings_if_valid
        )

    def _persist_detection_settings_if_valid(self):
        self._detection_settings_after_id = None

        try:
            params = self.get_detection_parameters()
        except Exception:
            return

        callback = getattr(
            self.master,
            "persist_stardist_settings",
            None
        )

        if callback is not None:
            callback(params)

    # ---------------- UI ----------------

    def build_ui(self):
        main = ttk.Frame(self)
        main.pack(fill=tk.BOTH, expand=True)

        canvas_frame = ttk.Frame(main)
        canvas_frame.pack(
            side=tk.LEFT,
            fill=tk.BOTH,
            expand=True
        )

        self.canvas = tk.Canvas(
            canvas_frame,
            bg="black",
            cursor="cross"
        )

        self.canvas.pack(
            fill=tk.BOTH,
            expand=True
        )

        panel_container = ttk.Frame(
            main,
            width=self.panel_width
        )

        self.panel_container = panel_container

        panel_container.pack(
            side=tk.RIGHT,
            fill=tk.Y,
            padx=4,
            pady=4
        )

        panel_container.pack_propagate(False)

        panel_canvas = tk.Canvas(
            panel_container,
            width=self.panel_width,
            highlightthickness=0
        )

        self.panel_canvas = panel_canvas

        panel_canvas.pack(
            side=tk.LEFT,
            fill=tk.BOTH,
            expand=True
        )

        scroll = ttk.Scrollbar(
            panel_container,
            orient=tk.VERTICAL,
            command=panel_canvas.yview
        )

        scroll.pack(
            side=tk.RIGHT,
            fill=tk.Y
        )

        panel_canvas.configure(
            yscrollcommand=scroll.set
        )

        panel = ttk.Frame(panel_canvas)

        pw = panel_canvas.create_window(
            (0, 0),
            window=panel,
            anchor="nw"
        )

        panel.bind(
            "<Configure>",
            lambda e: panel_canvas.configure(
                scrollregion=panel_canvas.bbox("all")
            )
        )

        panel_canvas.bind(
            "<Configure>",
            lambda e: panel_canvas.itemconfigure(
                pw,
                width=e.width
            )
        )

        # Right control panel: wheel scroll.
        # Left image canvas keeps its own wheel binding for zoom.
        self.bind(
            "<MouseWheel>",
            self._on_annotation_window_wheel,
            add="+"
        )
        self.bind(
            "<Button-4>",
            self._on_annotation_window_wheel,
            add="+"
        )
        self.bind(
            "<Button-5>",
            self._on_annotation_window_wheel,
            add="+"
        )

        # Sample
        info = ttk.LabelFrame(panel, text="Sample")
        info.pack(fill=tk.X, pady=5)

        for text in (
            f"Folder: {self.sample_name}",
            f"DAPI: {self.dapi_path.name}",
            f"ZO-1: {self.zo1_path.name}",
        ):
            ttk.Label(
                info,
                text=text,
                wraplength=245,
                justify=tk.LEFT
            ).pack(
                anchor="w",
                padx=5,
                pady=2
            )

        # StarDist
        detect = ttk.LabelFrame(
            panel,
            text="StarDist detection"
        )

        detect.pack(fill=tk.X, pady=5)

        entries = [
            ("1 pixel (µm)", self.pixel_size),
            ("Probability threshold", self.prob_thresh),
            ("NMS threshold", self.nms_thresh),
            ("Min nucleus area (µm²)", self.min_area_um2),
            ("Max nucleus area (µm²)", self.max_area_um2),
            ("Normalize p-min", self.normalize_pmin),
            ("Normalize p-max", self.normalize_pmax),
        ]

        for row, (label, var) in enumerate(entries):
            ttk.Label(
                detect,
                text=label
            ).grid(
                row=row,
                column=0,
                sticky="w",
                padx=5,
                pady=2
            )

            ttk.Entry(
                detect,
                textvariable=var,
                width=10
            ).grid(
                row=row,
                column=1,
                sticky="ew",
                padx=5,
                pady=2
            )

        detect.columnconfigure(1, weight=1)

        self.run_button = ttk.Button(
            detect,
            text="Run StarDist / replace detections",
            command=self.run_stardist
        )

        self.run_button.grid(
            row=len(entries),
            column=0,
            columnspan=2,
            sticky="ew",
            padx=5,
            pady=6
        )

        if self.paths["npz"].exists():
            initial_detect_text = "Reviewed result loaded"
            initial_detect_bg = "green"
        elif self.loaded_from_draft:
            initial_detect_text = "Batch StarDist draft loaded / 未確認"
            initial_detect_bg = "dark orange"
        else:
            initial_detect_text = "Not detected"
            initial_detect_bg = "gray"

        self.detect_status = tk.Label(
            detect,
            text=initial_detect_text,
            fg="white",
            bg=initial_detect_bg,
        )

        self.detect_status.grid(
            row=len(entries) + 1,
            column=0,
            columnspan=2,
            sticky="ew",
            padx=5,
            pady=4
        )

        # Display
        disp = ttk.LabelFrame(panel, text="Display")
        disp.pack(fill=tk.X, pady=5)

        ttk.Checkbutton(
            disp,
            text="Show nucleus contours",
            variable=self.show_contours,
            command=self.refresh_display_layer
        ).pack(
            fill=tk.X,
            padx=5,
            pady=2
        )

        ttk.Checkbutton(
            disp,
            text="Show nucleus IDs",
            variable=self.show_ids,
            command=self.refresh_display_layer
        ).pack(
            fill=tk.X,
            padx=5,
            pady=2
        )

        ttk.Button(
            disp,
            textvariable=self.zo1_overlay_button_text,
            command=self.toggle_zo1_overlay
        ).pack(
            fill=tk.X,
            padx=5,
            pady=(4, 6)
        )

        ttk.Label(
            disp,
            text="Brightness min"
        ).pack(
            anchor="w",
            padx=5
        )

        ttk.Entry(
            disp,
            textvariable=self.display_min
        ).pack(
            fill=tk.X,
            padx=5,
            pady=2
        )

        ttk.Label(
            disp,
            text="Brightness max"
        ).pack(
            anchor="w",
            padx=5
        )

        ttk.Entry(
            disp,
            textvariable=self.display_max
        ).pack(
            fill=tk.X,
            padx=5,
            pady=2
        )

        disp_row = ttk.Frame(disp)

        disp_row.pack(
            fill=tk.X,
            padx=5,
            pady=5
        )

        ttk.Button(
            disp_row,
            text="Apply",
            command=self.apply_display
        ).pack(
            side=tk.LEFT,
            fill=tk.X,
            expand=True
        )

        ttk.Button(
            disp_row,
            text="Auto",
            command=self.auto_display
        ).pack(
            side=tk.LEFT,
            fill=tk.X,
            expand=True,
            padx=(4, 0)
        )

        ttk.Button(
            disp,
            text="Fit image",
            command=self.fit_image
        ).pack(
            fill=tk.X,
            padx=5,
            pady=(0, 5)
        )

        # Edit
        edit = ttk.LabelFrame(
            panel,
            text="Manual correction"
        )

        edit.pack(fill=tk.X, pady=5)

        ttk.Button(
            edit,
            text="Delete selected",
            command=self.delete_selected
        ).pack(fill=tk.X, padx=5, pady=2)

        ttk.Button(
            edit,
            text="Merge selected",
            command=self.merge_selected
        ).pack(fill=tk.X, padx=5, pady=2)

        ttk.Button(
            edit,
            text="Fill hole selected",
            command=self.fill_holes_selected
        ).pack(fill=tk.X, padx=5, pady=2)

        ttk.Button(
            edit,
            text="Merge touching nuclei (8-neighbor)",
            command=self.merge_touching_nuclei
        ).pack(fill=tk.X, padx=5, pady=2)

        ttk.Button(
            edit,
            text="Split disconnected selected",
            command=self.split_selected
        ).pack(fill=tk.X, padx=5, pady=2)

        ttk.Button(
            edit,
            text="Undo",
            command=self.undo
        ).pack(fill=tk.X, padx=5, pady=2)

        ttk.Button(
            edit,
            text="Clear all nuclei",
            command=self.clear_all
        ).pack(
            fill=tk.X,
            padx=5,
            pady=(2, 5)
        )

        # Save
        save = ttk.LabelFrame(panel, text="Save")
        save.pack(fill=tk.X, pady=5)

        ttk.Button(
            save,
            text="保存 [S]",
            command=self.save_annotation
        ).pack(fill=tk.X, padx=5, pady=2)

        ttk.Button(
            save,
            text="保存して閉じる",
            command=self.save_and_close
        ).pack(fill=tk.X, padx=5, pady=2)

        ttk.Button(
            save,
            text="保存して次へ [N]",
            command=self.save_and_next
        ).pack(
            fill=tk.X,
            padx=5,
            pady=(2, 5)
        )

        # Status
        status = ttk.LabelFrame(
            panel,
            text="Status / mouse"
        )

        status.pack(
            fill=tk.BOTH,
            expand=True,
            pady=5
        )

        self.status_label = ttk.Label(
            status,
            text="",
            anchor="nw",
            justify=tk.LEFT,
            wraplength=245,
        )

        self.status_label.pack(
            fill=tk.X,
            padx=5,
            pady=5
        )

        ttk.Label(
            status,
            text=(
                "左ドラッグ: 新しい核を描画\n"
                "1核選択後の左ドラッグ: その核へ追加\n"
                "Ctrl+左ドラッグ: 消しゴム\n"
                "ダブルクリック: 核選択/解除\n"
                "右ドラッグ: 移動\n"
                "ホイール: 拡大縮小\n"
                "[ / ]: ブラシサイズ\n"
                "Delete: 選択核削除\n"
                "Ctrl+Z: Undo"
            ),
            anchor="nw",
            justify=tk.LEFT,
            wraplength=245,
        ).pack(
            fill=tk.X,
            padx=5,
            pady=5
        )

        # bindings
        self.canvas.bind(
            "<Configure>",
            lambda e: self.redraw(False)
        )

        self.canvas.bind("<MouseWheel>", self.on_wheel)
        self.canvas.bind("<Button-4>", self.on_wheel_linux)
        self.canvas.bind("<Button-5>", self.on_wheel_linux)

        self.canvas.bind("<ButtonPress-1>", self.on_left_press)
        self.canvas.bind("<B1-Motion>", self.on_left_motion)
        self.canvas.bind("<ButtonRelease-1>", self.on_left_release)
        self.canvas.bind("<Double-Button-1>", self.on_left_double_click)

        self.canvas.bind("<ButtonPress-3>", self.on_right_down)
        self.canvas.bind("<B3-Motion>", self.on_right_drag)
        self.canvas.bind("<ButtonRelease-3>", self.on_right_up)

        self.bind("<Delete>", lambda e: self.delete_selected())
        self.bind("<Control-z>", lambda e: self.undo())

        self.bind("[", lambda e: self.change_brush_size(-2))
        self.bind("]", lambda e: self.change_brush_size(2))

        self.bind("<KeyPress-s>", self.on_s)
        self.bind("<KeyPress-n>", self.on_n)

        self.update_status()

    # ---------------- display ----------------

    def redraw(self, force_resize=False):
        self.display.redraw(
            canvas=self.canvas,
            zoom=self.zoom,
            offset_x=self.offset_x,
            offset_y=self.offset_y,
            display_min=int(self.display_min.get()),
            display_max=int(self.display_max.get()),
            show_contours=bool(self.show_contours.get()),
            show_numbers=bool(self.show_ids.get()),
            show_zo1_overlay=bool(self.show_zo1_overlay),
            zo1_overlay_alpha=float(self.zo1_overlay_alpha),
        )

    def refresh_display_layer(self):
        self.display.invalidate_roi_layer()
        self.redraw(False)

    def toggle_zo1_overlay(self):
        """
        Toggle ZO-1 overlay.

        ZO-1 TIFF is loaded only on the first ON action to keep ordinary
        DAPI annotation startup and editing light.
        """
        if self.show_zo1_overlay:
            self.show_zo1_overlay = False
            self.zo1_overlay_button_text.set(
                "ZO-1 overlay: OFF"
            )
            self.redraw(False)
            return

        try:
            if self.display.zo1_gray8 is None:
                self.zo1_overlay_button_text.set(
                    "ZO-1 overlay: loading..."
                )
                self.update_idletasks()

                zo1_raw = read_tiff_gray(
                    self.zo1_path
                )

                if zo1_raw.shape != self.raw.shape:
                    raise ValueError(
                        "DAPIとZO-1画像サイズが一致しません。\n"
                        f"DAPI: {self.raw.shape}\n"
                        f"ZO-1: {zo1_raw.shape}"
                    )

                zo1_gray8 = normalize_8bit(
                    zo1_raw
                )

                # Only retain 8-bit display data to reduce RAM use.
                self.display.set_zo1_overlay_image(
                    zo1_gray8
                )

                del zo1_raw

            self.show_zo1_overlay = True
            self.zo1_overlay_button_text.set(
                "ZO-1 overlay: ON"
            )
            self.redraw(False)

        except Exception:
            self.show_zo1_overlay = False
            self.zo1_overlay_button_text.set(
                "ZO-1 overlay: OFF"
            )
            messagebox.showerror(
                "ZO-1 overlay",
                traceback.format_exc(),
                parent=self
            )

    def apply_display(self):
        try:
            vmin = int(self.display_min.get())
            vmax = int(self.display_max.get())
        except Exception:
            messagebox.showerror(
                "Display",
                "Brightness min/maxを整数で入力してください。",
                parent=self
            )
            return

        if vmax <= vmin:
            messagebox.showerror(
                "Display",
                "Brightness maxはminより大きくしてください。",
                parent=self
            )
            return

        self.display.invalidate_base()
        self.redraw(False)

    def auto_display(self):
        presets = [
            (0.1, 99.9),
            (1.0, 99.5),
            (2.0, 99.0),
            (5.0, 98.0),
            (10.0, 99.0),
        ]

        pmin, pmax = presets[self.auto_display_mode]

        self.display_min.set(
            int(np.percentile(self.gray8, pmin))
        )

        self.display_max.set(
            int(np.percentile(self.gray8, pmax))
        )

        if self.display_max.get() <= self.display_min.get():
            self.display_max.set(
                min(
                    255,
                    self.display_min.get() + 1
                )
            )

        self.auto_display_mode = (
            self.auto_display_mode + 1
        ) % len(presets)

        self.display.invalidate_base()
        self.redraw(False)

    def fit_image(self):
        self.update_idletasks()

        cw = max(1, self.canvas.winfo_width())
        ch = max(1, self.canvas.winfo_height())

        self.zoom = float(
            np.clip(
                min(cw / self.w, ch / self.h),
                0.05,
                20.0
            )
        )

        self.offset_x = (
            cw - self.w * self.zoom
        ) / 2

        self.offset_y = (
            ch - self.h * self.zoom
        ) / 2

        self.display.invalidate_resize()
        self.redraw(True)
        self.update_status()

    # ---------------- StarDist ----------------

    def get_detection_parameters(self):
        p = {
            "pixel_size_um": float(self.pixel_size.get()),
            "prob_thresh": float(self.prob_thresh.get()),
            "nms_thresh": float(self.nms_thresh.get()),
            "min_area_um2": float(self.min_area_um2.get()),
            "max_area_um2": float(self.max_area_um2.get()),
            "normalize_pmin": float(self.normalize_pmin.get()),
            "normalize_pmax": float(self.normalize_pmax.get()),
        }

        if p["pixel_size_um"] <= 0:
            raise ValueError("1 pixel (µm)は0より大きくしてください。")

        if not 0 < p["prob_thresh"] < 1:
            raise ValueError("Probability thresholdは0～1です。")

        if not 0 < p["nms_thresh"] < 1:
            raise ValueError("NMS thresholdは0～1です。")

        if (
            p["min_area_um2"] < 0
            or p["max_area_um2"] < 0
        ):
            raise ValueError("核面積は0以上にしてください。")

        if (
            p["max_area_um2"] > 0
            and p["max_area_um2"] <= p["min_area_um2"]
        ):
            raise ValueError(
                "Max nucleus areaはMin nucleus areaより大きくしてください。"
            )

        if not (
            0 <= p["normalize_pmin"]
            < p["normalize_pmax"]
            <= 100
        ):
            raise ValueError(
                "Normalize p-min < p-max、かつ0～100にしてください。"
            )

        return p

    def run_stardist(self, confirm_replace=True):
        if self.processing_running:
            return

        try:
            params = self.get_detection_parameters()
        except Exception as exc:
            messagebox.showerror(
                "Parameter error",
                str(exc),
                parent=self
            )
            return

        callback = getattr(
            self.master,
            "persist_stardist_settings",
            None
        )
        if callback is not None:
            callback(params)

        if (
            confirm_replace
            and self.manager.count() > 0
        ):
            if not messagebox.askyesno(
                "Replace detections",
                "現在の検出結果と手動修正を"
                "新しいStarDist結果で置き換えますか？",
                parent=self
            ):
                return

        self.processing_running = True

        self.run_button.config(
            state=tk.DISABLED
        )

        self.detect_status.config(
            text="StarDist processing...",
            bg="orange"
        )

        threading.Thread(
            target=self._stardist_worker,
            args=(params,),
            daemon=True
        ).start()

        self.after(
            100,
            self.poll_processing_queue
        )

    def _stardist_worker(self, params):
        try:
            labels = StarDistNucleiSegmenter.segment(
                self.raw,
                **params
            )

            self.processing_queue.put(
                ("done", labels)
            )

        except Exception:
            self.processing_queue.put(
                ("error", traceback.format_exc())
            )

    def poll_processing_queue(self):
        try:
            typ, data = self.processing_queue.get_nowait()
        except queue.Empty:
            if self.processing_running:
                self.after(
                    100,
                    self.poll_processing_queue
                )
            return

        self.processing_running = False

        self.run_button.config(
            state=tk.NORMAL
        )

        if typ == "error":
            self.detect_status.config(
                text="StarDist error",
                bg="red"
            )

            messagebox.showerror(
                "StarDist error",
                data,
                parent=self
            )
            return

        self.manager.push_undo()
        self.manager.set_label_mask(data)

        self.display.invalidate_all()

        self.unsaved_changes = True

        self.detect_status.config(
            text=(
                f"Detected: {self.manager.count()} nuclei "
                "(not saved)"
            ),
            bg="dark orange"
        )

        self.redraw(True)
        self.update_status()

    # ---------------- mouse ----------------

    def screen_to_image(self, sx, sy):
        ix = int(
            (sx - self.offset_x)
            / max(self.zoom, 1e-9)
        )

        iy = int(
            (sy - self.offset_y)
            / max(self.zoom, 1e-9)
        )

        return (
            int(np.clip(ix, 0, self.w - 1)),
            int(np.clip(iy, 0, self.h - 1)),
        )

    @staticmethod
    def is_ctrl(event):
        return bool(event.state & 0x0004)

    def on_left_press(self, event):
        if self.processing_running:
            return

        self.left_pressed = True
        self.left_started = False
        self.left_press_event = event

        if self.left_press_after_id is not None:
            self.after_cancel(
                self.left_press_after_id
            )

        self.left_press_after_id = self.after(
            self.double_click_delay_ms,
            self.start_left_action
        )

    def start_left_action(self):
        if (
            not self.left_pressed
            or self.left_press_event is None
        ):
            return

        self.left_started = True
        event = self.left_press_event

        self.manager.push_undo()

        if self.is_ctrl(event):
            self.drag_mode = "eraser"
            self.manager.eraser_touched_labels.clear()

            self.draw_brush(
                event.x,
                event.y,
                add=False
            )
        else:
            self.drag_mode = "pen"

            self.draw_brush(
                event.x,
                event.y,
                add=True
            )

    def on_left_motion(self, event):
        if (
            self.left_pressed
            and not self.left_started
        ):
            if self.left_press_after_id is not None:
                self.after_cancel(
                    self.left_press_after_id
                )

                self.left_press_after_id = None

            self.left_press_event = event
            self.start_left_action()

        if self.drag_mode == "pen":
            self.draw_brush(
                event.x,
                event.y,
                add=True
            )

        elif self.drag_mode == "eraser":
            self.draw_brush(
                event.x,
                event.y,
                add=False
            )

    def on_left_release(self, event):
        self.left_pressed = False

        if self.left_press_after_id is not None:
            self.after_cancel(
                self.left_press_after_id
            )

            self.left_press_after_id = None

        if not self.left_started:
            self.left_press_event = None
            return

        if self.drag_mode == "pen":
            self.manager.finish_pen()

        elif self.drag_mode == "eraser":
            self.manager.finish_eraser()

        if self.drag_mode in {"pen", "eraser"}:
            # Only ROI layer changed. Keep brightness/base cache.
            self.display.invalidate_roi_layer()

            self.unsaved_changes = True

            self.detect_status.config(
                text="Manually edited (not saved)",
                bg="dark orange"
            )

            self.redraw(False)

        self.drag_mode = None
        self.left_press_event = None

        self.update_status()

    def draw_brush(self, sx, sy, add=True):
        ix, iy = self.screen_to_image(
            sx,
            sy
        )

        self.manager.draw_brush(
            ix,
            iy,
            self.brush_screen_radius,
            self.zoom,
            add=add,
        )

        # During dragging, only show the brush marker.
        # Full image redraw occurs on mouse release.
        color = "red" if add else "blue"
        r = self.brush_screen_radius

        self.canvas.create_oval(
            sx - r,
            sy - r,
            sx + r,
            sy + r,
            outline=color,
            fill=color,
            stipple="gray50",
            tags="brush_preview"
        )

    def on_left_double_click(self, event):
        if self.left_press_after_id is not None:
            self.after_cancel(
                self.left_press_after_id
            )

            self.left_press_after_id = None

        self.left_pressed = False
        self.left_started = False
        self.drag_mode = None

        ix, iy = self.screen_to_image(
            event.x,
            event.y
        )

        lab = int(
            self.manager.label_mask[iy, ix]
        )

        if lab <= 0:
            self.manager.selected_labels.clear()
        elif lab in self.manager.selected_labels:
            self.manager.selected_labels.remove(lab)
        else:
            self.manager.selected_labels.add(lab)

        # Contours themselves are cached; selection only changes color.
        self.redraw(False)
        self.update_status()

    def on_right_down(self, event):
        self.last_drag = (
            event.x,
            event.y
        )
        self.drag_mode = "pan"

    def on_right_drag(self, event):
        if (
            self.drag_mode != "pan"
            or self.last_drag is None
        ):
            return

        dx = event.x - self.last_drag[0]
        dy = event.y - self.last_drag[1]

        # Ignore tiny jitter as in reference code.
        if abs(dx) + abs(dy) < 3:
            return

        self.offset_x += dx
        self.offset_y += dy

        self.last_drag = (
            event.x,
            event.y
        )

        self.redraw(False)

    def on_right_up(self, event):
        if self.drag_mode == "pan":
            self.drag_mode = None
            self.last_drag = None

    def on_wheel(self, event):
        factor = (
            1.15
            if event.delta > 0
            else 1 / 1.15
        )

        self.zoom_at(
            event.x,
            event.y,
            factor
        )

    def on_wheel_linux(self, event):
        factor = (
            1.15
            if event.num == 4
            else 1 / 1.15
        )

        self.zoom_at(
            event.x,
            event.y,
            factor
        )

    def zoom_at(self, sx, sy, factor):
        old_zoom = self.zoom

        new_zoom = float(
            np.clip(
                self.zoom * factor,
                0.1,
                20.0
            )
        )

        scale = new_zoom / old_zoom

        self.offset_x = (
            sx
            - scale * (sx - self.offset_x)
        )

        self.offset_y = (
            sy
            - scale * (sy - self.offset_y)
        )

        self.zoom = new_zoom

        self.display.invalidate_resize()
        self.redraw(True)
        self.update_status()

    def change_brush_size(self, delta):
        self.brush_screen_radius = int(
            np.clip(
                self.brush_screen_radius + delta,
                2,
                80
            )
        )

        self.update_status()

    # ---------------- commands ----------------

    def mark_edited(self, full=False):
        if full:
            self.display.invalidate_all()
        else:
            self.display.invalidate_roi_layer()

        self.unsaved_changes = True

        self.detect_status.config(
            text="Manually edited (not saved)",
            bg="dark orange"
        )

        self.redraw(False)
        self.update_status()

    def delete_selected(self):
        if self.manager.delete_selected():
            self.mark_edited(False)

    def merge_selected(self):
        if self.manager.merge_selected():
            self.mark_edited(False)

    def fill_holes_selected(self):
        if not self.manager.selected_labels:
            messagebox.showinfo(
                "Fill hole",
                "Fill holeする核を選択してください。",
                parent=self
            )
            return

        try:
            changed, added_pixels = (
                self.manager.fill_holes_selected()
            )
        except Exception:
            messagebox.showerror(
                "Fill hole",
                traceback.format_exc(),
                parent=self
            )
            return

        if not changed:
            messagebox.showinfo(
                "Fill hole",
                "選択核に埋める内部holeはありませんでした。",
                parent=self
            )
            return

        self.mark_edited(False)

        self.status_label.config(
            text=(
                self.manager.get_status_text(
                    self.brush_screen_radius
                )
                + f"\nFill hole: +{added_pixels} px"
                + f"\nZoom: {self.zoom:.2f}x"
            )
        )

    def merge_touching_nuclei(self):
        groups = self.manager.find_touching_groups()

        if not groups:
            messagebox.showinfo(
                "Merge touching nuclei",
                "直接接している別ラベル核はありませんでした。",
                parent=self
            )
            return

        n_groups = len(groups)
        n_labels = sum(len(g) for g in groups)
        n_removed = sum(len(g) - 1 for g in groups)

        preview = "\n".join(
            " + ".join(map(str, g))
            for g in groups[:8]
        )

        if len(groups) > 8:
            preview += f"\n... and {len(groups) - 8} more groups"

        ok = messagebox.askyesno(
            "Merge touching nuclei",
            (
                "8近傍で直接接している核ラベルをMergeします。\n"
                "距離thresholdやdilationは使いません。\n\n"
                f"Touching groups: {n_groups}\n"
                f"Labels involved: {n_labels}\n"
                f"Labels removed by merge: {n_removed}\n\n"
                f"{preview}\n\n"
                "実行しますか？"
            ),
            parent=self
        )

        if not ok:
            return

        changed, merged_groups, removed = (
            self.manager.merge_touching_groups(
                groups
            )
        )

        if changed:
            # Old label contour cache includes labels that disappeared.
            # invalidate_roi_layer + update_dirty_roi_meta will prune them.
            self.mark_edited(False)

            self.status_label.config(
                text=(
                    self.manager.get_status_text(
                        self.brush_screen_radius
                    )
                    + f"\nMerged touching: "
                    f"{merged_groups} groups / "
                    f"{removed} labels"
                    + f"\nZoom: {self.zoom:.2f}x"
                )
            )

    def split_selected(self):
        ok, msg = (
            self.manager.split_selected_disconnected()
        )

        if not ok:
            messagebox.showinfo(
                "Split",
                msg,
                parent=self
            )
            return

        self.mark_edited(False)

    def undo(self):
        if self.manager.undo():
            self.mark_edited(False)

    def clear_all(self):
        if not messagebox.askyesno(
            "Clear",
            "すべての核検出を削除しますか？",
            parent=self
        ):
            return

        if self.manager.clear_all():
            self.display.roi_meta_cache.clear()
            self.mark_edited(False)

    def update_status(self):
        text = self.manager.get_status_text(
            self.brush_screen_radius
        )

        self.status_label.config(
            text=(
                text
                + f"\nZoom: {self.zoom:.2f}x"
            )
        )

    # ---------------- save ----------------

    def save_overlay(self):
        # Fast debug overlay: label transition boundaries.
        rgb = cv2.cvtColor(
            self.gray8,
            cv2.COLOR_GRAY2RGB
        )

        L = self.manager.label_mask

        boundary = np.zeros(
            L.shape,
            dtype=bool
        )

        boundary[1:, :] |= (
            L[1:, :] != L[:-1, :]
        ) & (
            (L[1:, :] > 0)
            | (L[:-1, :] > 0)
        )

        boundary[:, 1:] |= (
            L[:, 1:] != L[:, :-1]
        ) & (
            (L[:, 1:] > 0)
            | (L[:, :-1] > 0)
        )

        rgb[boundary] = (
            0,
            255,
            0
        )

        imwrite_unicode(
            str(self.paths["overlay"]),
            cv2.cvtColor(
                rgb,
                cv2.COLOR_RGB2BGR
            )
        )

    def save_annotation(self):
        if self.processing_running:
            messagebox.showwarning(
                "Processing",
                "StarDist処理中は保存できません。",
                parent=self
            )
            return False

        try:
            params = self.get_detection_parameters()
        except Exception as exc:
            messagebox.showerror(
                "Parameter error",
                str(exc),
                parent=self
            )
            return False

        try:
            self.paths["dir"].mkdir(
                parents=True,
                exist_ok=True
            )

            self.manager.save(
                self.paths["npz"],
                self.sample_name,
                self.dapi_path.name,
                self.zo1_path.name,
                params
            )

            tifffile.imwrite(
                str(self.paths["labels"]),
                self.manager.label_mask.astype(
                    np.uint16
                )
            )

            tifffile.imwrite(
                str(self.paths["binary"]),
                (
                    self.manager.label_mask > 0
                ).astype(np.uint8) * 255
            )

            self.save_overlay()

            # Human review is now finalized, so the automatic draft is obsolete.
            if self.paths["draft_npz"].exists():
                try:
                    self.paths["draft_npz"].unlink()
                except Exception:
                    pass

            self.loaded_from_draft = False
            self.unsaved_changes = False

            self.detect_status.config(
                text=(
                    f"Saved: "
                    f"{self.manager.count()} nuclei"
                ),
                bg="green"
            )

            if self.saved_callback:
                self.saved_callback(
                    self.sample,
                    params
                )

            return True

        except Exception:
            messagebox.showerror(
                "Save error",
                traceback.format_exc(),
                parent=self
            )
            return False

    def save_and_close(self):
        if self.save_annotation():
            self.destroy()

    def save_and_next(self):
        if not self.save_annotation():
            return

        sample = self.sample
        callback = self.next_callback
        parent = self.master

        self.destroy()

        if callback:
            parent.after(
                50,
                lambda: callback(sample)
            )

    def on_s(self, event):
        widget = self.focus_get()

        if isinstance(
            widget,
            (tk.Entry, ttk.Entry, tk.Text)
        ):
            return

        self.save_annotation()
        return "break"

    def on_n(self, event):
        widget = self.focus_get()

        if isinstance(
            widget,
            (tk.Entry, ttk.Entry, tk.Text)
        ):
            return

        self.save_and_next()
        return "break"

    def on_window_close(self):
        self._persist_detection_settings_if_valid()

        if self.unsaved_changes:
            ans = messagebox.askyesnocancel(
                "Close DAPI annotation",
                "DAPI核検出・修正結果を保存して閉じますか？",
                parent=self
            )

            if ans is None:
                return

            if ans:
                if not self.save_annotation():
                    return

        self.destroy()


# ============================================================
# ZO-1 processing QC / parameter tuning window
# ============================================================

class ZO1ProcessingQCWindow(tk.Toplevel):
    """
    Representative-image GUI for checking ZO-1 preprocessing one step at a time.

    Stages:
        1. Raw ZO-1
        2. Rolling-ball background
        3. Background-subtracted
        4. White top-hat
        5. Threshold binary
        6. Skeleton
        7. Raw + skeleton
        8. Top-hat + skeleton
        9. DAPI + reviewed nuclei + skeleton

    The window does not change analysis files by itself.
    "この設定を本解析へ反映" sends the tuned parameters back to MainApp.
    """

    STAGES = [
        "1. Raw ZO-1",
        "2. Rolling-ball background",
        "3. Background-subtracted",
        "4. White top-hat",
        "5. Threshold binary",
        "6. Skeleton",
        "7. Raw + skeleton",
        "8. Top-hat + skeleton",
        "9. DAPI + nuclei + skeleton",
    ]

    def __init__(
        self,
        master,
        sample,
        initial_params,
        apply_callback=None,
    ):
        super().__init__(master)

        self.master_app = master
        self.sample = sample
        self.folder = Path(sample["folder"])
        self.sample_name = sample_display_name(sample)
        self.zo1_path = Path(sample["zo1"])
        self.dapi_path = (
            Path(sample["dapi"])
            if sample.get("dapi") is not None
            else None
        )
        self.apply_callback = apply_callback

        self.title(
            f"ZO-1 preprocessing QC - {self.sample_name}"
        )
        self.geometry("1450x900")
        self.minsize(1000, 650)
        self.panel_width = 360

        # Load raw ZO-1 once. This is the representative image.
        self.zo1_raw = read_tiff_gray(
            self.zo1_path
        ).astype(np.float32)

        self.h, self.w = self.zo1_raw.shape

        # Optional DAPI/reviewed labels for final overlay stage.
        self.dapi_raw = None
        self.nuclei_labels = None
        self.nuclei_regions = None

        if (
            self.dapi_path is not None
            and self.dapi_path.exists()
        ):
            try:
                d = read_tiff_gray(
                    self.dapi_path
                )

                if d.shape == self.zo1_raw.shape:
                    self.dapi_raw = d
            except Exception:
                self.dapi_raw = None

        ann_path = sample_output_paths(
            self.sample
        )["npz"]

        if ann_path.exists():
            try:
                ann = np.load(
                    ann_path,
                    allow_pickle=True
                )
                if "nuclei_labels" in ann:
                    labs = ann[
                        "nuclei_labels"
                    ].astype(np.uint16)
                    if labs.shape == self.zo1_raw.shape:
                        self.nuclei_labels = labs
                        self.nuclei_regions = (
                            _nucleus_regions_from_labels(
                                labs
                            )
                        )
            except Exception:
                self.nuclei_labels = None
                self.nuclei_regions = None

        p = dict(initial_params or {})

        self.rolling_ball_radius_var = tk.IntVar(
            value=int(
                p.get(
                    "zo1_rolling_ball_radius_px",
                    DEFAULT_ZO1_ROLLING_BALL_RADIUS_PX
                )
            )
        )

        self.tophat_radius_var = tk.IntVar(
            value=int(
                p.get(
                    "zo1_tophat_radius_px",
                    DEFAULT_ZO1_TOPHAT_RADIUS_PX
                )
            )
        )

        self.threshold_var = tk.DoubleVar(
            value=float(
                p.get(
                    "zo1_tophat_threshold",
                    DEFAULT_ZO1_TOPHAT_THRESHOLD
                )
            )
        )

        self.min_object_size_var = tk.IntVar(
            value=int(
                p.get(
                    "zo1_min_object_size_px",
                    DEFAULT_ZO1_MIN_OBJECT_SIZE_PX
                )
            )
        )

        self.stage_var = tk.StringVar(
            value=self.STAGES[0]
        )

        self.stage_images = {}
        self.stage_stats = {}

        # Numeric arrays for mouse-position inspection.
        # These retain the actual analysis-scale values, not display-normalized values.
        self.qc_arrays = {}

        self.zoom = 1.0
        self.offset_x = 0.0
        self.offset_y = 0.0
        self.last_drag = None
        self.tk_img = None

        self.processing_queue = queue.Queue()
        self.processing_running = False
        self.compute_generation = 0

        self.build_ui()

        self.protocol(
            "WM_DELETE_WINDOW",
            self.destroy
        )

        # Initial compute after widgets exist.
        self.after(
            150,
            self.recompute
        )

    # --------------------------------------------------------
    # UI
    # --------------------------------------------------------

    def build_ui(self):
        main = ttk.Frame(self)
        main.pack(
            fill=tk.BOTH,
            expand=True
        )

        canvas_frame = ttk.Frame(main)
        canvas_frame.pack(
            side=tk.LEFT,
            fill=tk.BOTH,
            expand=True
        )

        self.canvas = tk.Canvas(
            canvas_frame,
            bg="black",
            cursor="cross"
        )
        self.canvas.pack(
            fill=tk.BOTH,
            expand=True
        )

        panel_container = ttk.Frame(
            main,
            width=self.panel_width
        )
        panel_container.pack(
            side=tk.RIGHT,
            fill=tk.Y,
            padx=4,
            pady=4
        )
        panel_container.pack_propagate(False)

        panel_canvas = tk.Canvas(
            panel_container,
            width=self.panel_width,
            highlightthickness=0
        )
        panel_canvas.pack(
            side=tk.LEFT,
            fill=tk.BOTH,
            expand=True
        )

        scroll = ttk.Scrollbar(
            panel_container,
            orient=tk.VERTICAL,
            command=panel_canvas.yview
        )
        scroll.pack(
            side=tk.RIGHT,
            fill=tk.Y
        )

        panel_canvas.configure(
            yscrollcommand=scroll.set
        )

        panel = ttk.Frame(panel_canvas)
        pw = panel_canvas.create_window(
            (0, 0),
            window=panel,
            anchor="nw"
        )

        panel.bind(
            "<Configure>",
            lambda e: panel_canvas.configure(
                scrollregion=panel_canvas.bbox("all")
            )
        )
        panel_canvas.bind(
            "<Configure>",
            lambda e: panel_canvas.itemconfigure(
                pw,
                width=e.width
            )
        )

        # Sample info
        info = ttk.LabelFrame(
            panel,
            text="Representative image"
        )
        info.pack(
            fill=tk.X,
            pady=5
        )

        ttk.Label(
            info,
            text=f"Folder: {self.sample_name}",
            wraplength=330,
            justify=tk.LEFT
        ).pack(
            anchor="w",
            padx=5,
            pady=2
        )

        ttk.Label(
            info,
            text=f"ZO-1: {self.zo1_path.name}",
            wraplength=330,
            justify=tk.LEFT
        ).pack(
            anchor="w",
            padx=5,
            pady=2
        )

        nuclei_msg = (
            f"Reviewed nuclei: "
            f"{len(self.nuclei_regions)}"
            if self.nuclei_regions is not None
            else "Reviewed nuclei: not loaded"
        )

        ttk.Label(
            info,
            text=nuclei_msg
        ).pack(
            anchor="w",
            padx=5,
            pady=(2, 5)
        )

        # Parameters
        prm = ttk.LabelFrame(
            panel,
            text="ZO-1 preprocessing parameters"
        )
        prm.pack(
            fill=tk.X,
            pady=5
        )

        fields = [
            (
                "Rolling ball radius",
                self.rolling_ball_radius_var,
                "px"
            ),
            (
                "White top-hat radius",
                self.tophat_radius_var,
                "px"
            ),
            (
                "Top-hat threshold",
                self.threshold_var,
                "raw"
            ),
            (
                "Min object size",
                self.min_object_size_var,
                "px"
            ),
        ]

        for r, (name, var, unit) in enumerate(fields):
            ttk.Label(
                prm,
                text=name
            ).grid(
                row=r,
                column=0,
                sticky="w",
                padx=5,
                pady=3
            )

            ttk.Entry(
                prm,
                textvariable=var,
                width=10
            ).grid(
                row=r,
                column=1,
                sticky="ew",
                padx=5,
                pady=3
            )

            ttk.Label(
                prm,
                text=unit
            ).grid(
                row=r,
                column=2,
                sticky="w",
                padx=(0, 5),
                pady=3
            )

        prm.columnconfigure(
            1,
            weight=1
        )

        self.recompute_button = ttk.Button(
            prm,
            text="再計算 / Apply parameters",
            command=self.recompute
        )
        self.recompute_button.grid(
            row=len(fields),
            column=0,
            columnspan=3,
            sticky="ew",
            padx=5,
            pady=(6, 4)
        )

        ttk.Label(
            prm,
            text=(
                "Thresholdはwhite top-hat後の"
                "raw輝度スケールに適用します。"
            ),
            wraplength=325,
            justify=tk.LEFT
        ).grid(
            row=len(fields) + 1,
            column=0,
            columnspan=3,
            sticky="w",
            padx=5,
            pady=(0, 5)
        )

        # Stage
        stage_box = ttk.LabelFrame(
            panel,
            text="Processing step"
        )
        stage_box.pack(
            fill=tk.X,
            pady=5
        )

        self.stage_combo = ttk.Combobox(
            stage_box,
            textvariable=self.stage_var,
            values=self.STAGES,
            state="readonly"
        )
        self.stage_combo.pack(
            fill=tk.X,
            padx=5,
            pady=5
        )
        self.stage_combo.bind(
            "<<ComboboxSelected>>",
            lambda e: self.redraw()
        )

        nav = ttk.Frame(stage_box)
        nav.pack(
            fill=tk.X,
            padx=5,
            pady=(0, 5)
        )

        ttk.Button(
            nav,
            text="← 前の処理",
            command=lambda: self.change_stage(-1)
        ).pack(
            side=tk.LEFT,
            fill=tk.X,
            expand=True
        )

        ttk.Button(
            nav,
            text="次の処理 →",
            command=lambda: self.change_stage(1)
        ).pack(
            side=tk.LEFT,
            fill=tk.X,
            expand=True,
            padx=(4, 0)
        )

        ttk.Button(
            stage_box,
            text="Fit image",
            command=self.fit_image
        ).pack(
            fill=tk.X,
            padx=5,
            pady=(0, 5)
        )

        # Stats
        stat = ttk.LabelFrame(
            panel,
            text="QC / statistics"
        )
        stat.pack(
            fill=tk.X,
            pady=5
        )

        self.status_label = ttk.Label(
            stat,
            text="未計算",
            justify=tk.LEFT,
            anchor="nw",
            wraplength=325
        )
        self.status_label.pack(
            fill=tk.X,
            padx=5,
            pady=5
        )

        ttk.Separator(
            stat
        ).pack(
            fill=tk.X,
            padx=5,
            pady=4
        )

        self.mouse_value_label = ttk.Label(
            stat,
            text="Mouse: 画像上にカーソルを置いてください",
            justify=tk.LEFT,
            anchor="nw",
            wraplength=325
        )
        self.mouse_value_label.pack(
            fill=tk.X,
            padx=5,
            pady=(2, 6)
        )

        # Apply to analysis
        apply_box = ttk.LabelFrame(
            panel,
            text="Use settings"
        )
        apply_box.pack(
            fill=tk.X,
            pady=5
        )

        ttk.Button(
            apply_box,
            text="この設定を本解析へ反映",
            command=self.apply_to_main
        ).pack(
            fill=tk.X,
            padx=5,
            pady=5
        )

        ttk.Label(
            apply_box,
            text=(
                "反映後、親GUIの「選択サンプルを解析」または"
                "「全保存済みサンプルを解析」を実行してください。"
            ),
            wraplength=325,
            justify=tk.LEFT
        ).pack(
            fill=tk.X,
            padx=5,
            pady=(0, 5)
        )

        # Mouse
        help_box = ttk.LabelFrame(
            panel,
            text="View"
        )
        help_box.pack(
            fill=tk.X,
            pady=5
        )

        ttk.Label(
            help_box,
            text=(
                "右ドラッグ: pan\n"
                "Mouse wheel: zoom\n"
                "各処理は同じviewportで比較できます。"
            ),
            justify=tk.LEFT
        ).pack(
            fill=tk.X,
            padx=5,
            pady=5
        )

        self.canvas.bind(
            "<Configure>",
            lambda e: self.redraw()
        )
        self.canvas.bind(
            "<MouseWheel>",
            self.on_wheel
        )
        self.canvas.bind(
            "<Button-4>",
            self.on_wheel_linux
        )
        self.canvas.bind(
            "<Button-5>",
            self.on_wheel_linux
        )
        self.canvas.bind(
            "<ButtonPress-3>",
            self.on_right_down
        )
        self.canvas.bind(
            "<B3-Motion>",
            self.on_right_drag
        )
        self.canvas.bind(
            "<ButtonRelease-3>",
            self.on_right_up
        )

        self.canvas.bind(
            "<Motion>",
            self.on_mouse_move
        )

        self.canvas.bind(
            "<Leave>",
            self.on_mouse_leave
        )

    # --------------------------------------------------------
    # Parameters / processing
    # --------------------------------------------------------

    def get_params(self):
        try:
            p = {
                "zo1_rolling_ball_radius_px":
                    int(
                        self.rolling_ball_radius_var.get()
                    ),
                "zo1_tophat_radius_px":
                    int(
                        self.tophat_radius_var.get()
                    ),
                "zo1_tophat_threshold":
                    float(
                        self.threshold_var.get()
                    ),
                "zo1_min_object_size_px":
                    int(
                        self.min_object_size_var.get()
                    ),
            }
        except Exception as exc:
            raise ValueError(
                f"パラメータを確認してください: {exc}"
            ) from exc

        if p["zo1_rolling_ball_radius_px"] < 1:
            raise ValueError(
                "Rolling ball radiusは1 px以上にしてください。"
            )

        if p["zo1_tophat_radius_px"] < 1:
            raise ValueError(
                "White top-hat radiusは1 px以上にしてください。"
            )

        if p["zo1_tophat_threshold"] < 0:
            raise ValueError(
                "Top-hat thresholdは0以上にしてください。"
            )

        if p["zo1_min_object_size_px"] < 1:
            raise ValueError(
                "Min object sizeは1 px以上にしてください。"
            )

        return p

    def recompute(self):
        if self.processing_running:
            return

        try:
            p = self.get_params()
        except Exception as exc:
            messagebox.showerror(
                "ZO-1 QC",
                str(exc),
                parent=self
            )
            return

        self.processing_running = True
        self.compute_generation += 1
        generation = self.compute_generation

        self.recompute_button.config(
            state=tk.DISABLED
        )
        self.status_label.config(
            text="処理中..."
        )

        thread = threading.Thread(
            target=self._compute_worker,
            args=(p, generation),
            daemon=True
        )
        thread.start()

        self.after(
            100,
            self.poll_processing
        )

    def _compute_worker(
        self,
        p,
        generation
    ):
        try:
            raw = self.zo1_raw

            background = (
                _rolling_ball_background_fast(
                    raw,
                    radius_px=p[
                        "zo1_rolling_ball_radius_px"
                    ]
                )
            )

            bg_sub = (
                raw.astype(np.float32)
                - background.astype(np.float32)
            )
            np.maximum(
                bg_sub,
                0,
                out=bg_sub
            )

            top_hat = _white_tophat_raw(
                bg_sub,
                radius_px=p[
                    "zo1_tophat_radius_px"
                ]
            )

            binary = (
                top_hat
                >= p["zo1_tophat_threshold"]
            )

            if remove_small_objects is None:
                raise RuntimeError(
                    "scikit-imageが必要です。\n"
                    "pip install scikit-image"
                )

            binary = remove_small_objects(
                binary,
                min_size=max(
                    1,
                    p["zo1_min_object_size_px"]
                )
            )

            if skeletonize is None:
                raise RuntimeError(
                    "scikit-imageが必要です。\n"
                    "pip install scikit-image"
                )

            skel = skeletonize(binary)

            raw8 = normalize_8bit(raw)
            background8 = normalize_8bit(
                background
            )
            bg_sub8 = normalize_8bit(
                bg_sub
            )
            top_hat8 = normalize_8bit(
                top_hat
            )

            binary8 = (
                binary.astype(np.uint8)
                * 255
            )

            skeleton8 = (
                skel.astype(np.uint8)
                * 255
            )

            # Slightly dilated skeleton ONLY for display overlays.
            skel_vis = cv2.dilate(
                skel.astype(np.uint8),
                np.ones(
                    (3, 3),
                    dtype=np.uint8
                ),
                iterations=1
            ).astype(bool)

            raw_skel = cv2.cvtColor(
                raw8,
                cv2.COLOR_GRAY2RGB
            )
            raw_skel[
                skel_vis
            ] = (
                255,
                0,
                0
            )

            top_skel = cv2.cvtColor(
                top_hat8,
                cv2.COLOR_GRAY2RGB
            )
            top_skel[
                skel_vis
            ] = (
                255,
                0,
                0
            )

            # DAPI + nuclei + skeleton
            if self.dapi_raw is not None:
                d8 = normalize_8bit(
                    self.dapi_raw
                )

                dapi_overlay = np.zeros(
                    (
                        self.h,
                        self.w,
                        3
                    ),
                    dtype=np.uint8
                )

                # RGB: blue-dominant DAPI
                dapi_overlay[..., 2] = d8
                dapi_overlay[..., 1] = (
                    0.25 * d8
                ).astype(np.uint8)
                dapi_overlay[..., 0] = (
                    0.10 * d8
                ).astype(np.uint8)

                # ZO-1 skeleton red
                dapi_overlay[
                    skel_vis
                ] = (
                    255,
                    0,
                    0
                )

                if (
                    self.nuclei_labels is not None
                    and self.nuclei_regions is not None
                ):
                    for lab, region in (
                        self.nuclei_regions.items()
                    ):
                        x1, y1, x2, y2 = (
                            region["bbox"]
                        )

                        local = (
                            self.nuclei_labels[
                                y1:y2,
                                x1:x2
                            ]
                            == lab
                        ).astype(np.uint8)

                        contours, _ = (
                            cv2.findContours(
                                local,
                                cv2.RETR_EXTERNAL,
                                cv2.CHAIN_APPROX_SIMPLE
                            )
                        )

                        shifted = []

                        for cnt in contours:
                            c = cnt.copy()
                            c[:, 0, 0] += x1
                            c[:, 0, 1] += y1
                            shifted.append(c)

                        # RGB green
                        cv2.drawContours(
                            dapi_overlay,
                            shifted,
                            -1,
                            (0, 255, 0),
                            1
                        )

                        cx, cy = (
                            region["centroid_xy"]
                        )

                        # RGB yellow
                        cv2.circle(
                            dapi_overlay,
                            (
                                int(round(cx)),
                                int(round(cy))
                            ),
                            2,
                            (255, 255, 0),
                            -1,
                            cv2.LINE_AA
                        )
            else:
                dapi_overlay = raw_skel.copy()

            stages = {
                self.STAGES[0]:
                    cv2.cvtColor(
                        raw8,
                        cv2.COLOR_GRAY2RGB
                    ),
                self.STAGES[1]:
                    cv2.cvtColor(
                        background8,
                        cv2.COLOR_GRAY2RGB
                    ),
                self.STAGES[2]:
                    cv2.cvtColor(
                        bg_sub8,
                        cv2.COLOR_GRAY2RGB
                    ),
                self.STAGES[3]:
                    cv2.cvtColor(
                        top_hat8,
                        cv2.COLOR_GRAY2RGB
                    ),
                self.STAGES[4]:
                    cv2.cvtColor(
                        binary8,
                        cv2.COLOR_GRAY2RGB
                    ),
                self.STAGES[5]:
                    cv2.cvtColor(
                        skeleton8,
                        cv2.COLOR_GRAY2RGB
                    ),
                self.STAGES[6]:
                    raw_skel,
                self.STAGES[7]:
                    top_skel,
                self.STAGES[8]:
                    dapi_overlay,
            }

            stats = {
                "raw_min":
                    float(np.min(raw)),
                "raw_max":
                    float(np.max(raw)),
                "background_min":
                    float(np.min(background)),
                "background_max":
                    float(np.max(background)),
                "bg_sub_max":
                    float(np.max(bg_sub)),
                "tophat_max":
                    float(np.max(top_hat)),
                "threshold":
                    float(
                        p["zo1_tophat_threshold"]
                    ),
                "binary_pixels":
                    int(np.sum(binary)),
                "binary_fraction":
                    float(np.mean(binary)),
                "skeleton_pixels":
                    int(np.sum(skel)),
            }

            qc_arrays = {
                "raw_zo1": raw.astype(np.float32, copy=False),
                "background": background.astype(np.float32, copy=False),
                "background_subtracted": bg_sub.astype(np.float32, copy=False),
                "white_tophat": top_hat.astype(np.float32, copy=False),
                "binary": binary.astype(bool, copy=False),
                "skeleton": skel.astype(bool, copy=False),
            }

            if self.dapi_raw is not None:
                qc_arrays["dapi"] = np.asarray(
                    self.dapi_raw,
                    dtype=np.float32
                )

            self.processing_queue.put({
                "type": "done",
                "generation": generation,
                "stages": stages,
                "stats": stats,
                "params": p,
                "qc_arrays": qc_arrays,
            })

        except Exception:
            self.processing_queue.put({
                "type": "error",
                "generation": generation,
                "error": traceback.format_exc(),
            })

    def poll_processing(self):
        try:
            msg = (
                self.processing_queue
                .get_nowait()
            )
        except queue.Empty:
            if self.processing_running:
                self.after(
                    100,
                    self.poll_processing
                )
            return

        # Ignore stale result if a newer compute was launched.
        if msg.get(
            "generation"
        ) != self.compute_generation:
            if self.processing_running:
                self.after(
                    100,
                    self.poll_processing
                )
            return

        self.processing_running = False

        self.recompute_button.config(
            state=tk.NORMAL
        )

        if msg["type"] == "error":
            self.status_label.config(
                text="処理エラー"
            )
            messagebox.showerror(
                "ZO-1 QC",
                msg["error"],
                parent=self
            )
            return

        self.stage_images = msg["stages"]
        self.stage_stats = msg["stats"]
        self.qc_arrays = msg.get(
            "qc_arrays",
            {}
        )
        self.last_computed_params = msg["params"]

        self.update_stats_label()
        self.fit_image()

    def update_stats_label(self):
        s = self.stage_stats

        if not s:
            self.status_label.config(
                text="未計算"
            )
            return

        self.status_label.config(
            text=(
                f"Raw min/max: "
                f"{s['raw_min']:.1f} / "
                f"{s['raw_max']:.1f}\n"
                f"Background min/max: "
                f"{s['background_min']:.1f} / "
                f"{s['background_max']:.1f}\n"
                f"BG-sub max: "
                f"{s['bg_sub_max']:.1f}\n"
                f"Top-hat max: "
                f"{s['tophat_max']:.1f}\n"
                f"Threshold: "
                f"{s['threshold']:.2f}\n"
                f"Binary area: "
                f"{100.0 * s['binary_fraction']:.3f}% "
                f"({s['binary_pixels']} px)\n"
                f"Skeleton: "
                f"{s['skeleton_pixels']} px"
            )
        )

    # --------------------------------------------------------
    # Stage / view
    # --------------------------------------------------------

    def change_stage(self, delta):
        current = self.stage_var.get()

        try:
            idx = self.STAGES.index(
                current
            )
        except ValueError:
            idx = 0

        idx = int(
            np.clip(
                idx + int(delta),
                0,
                len(self.STAGES) - 1
            )
        )

        self.stage_var.set(
            self.STAGES[idx]
        )
        self.redraw()

    def get_current_image(self):
        return self.stage_images.get(
            self.stage_var.get()
        )

    def fit_image(self):
        img = self.get_current_image()

        if img is None:
            return

        self.update_idletasks()

        cw = max(
            1,
            self.canvas.winfo_width()
        )
        ch = max(
            1,
            self.canvas.winfo_height()
        )

        ih, iw = img.shape[:2]

        self.zoom = float(
            np.clip(
                min(
                    cw / max(1, iw),
                    ch / max(1, ih)
                ),
                0.05,
                20.0
            )
        )

        self.offset_x = (
            cw - iw * self.zoom
        ) / 2.0

        self.offset_y = (
            ch - ih * self.zoom
        ) / 2.0

        self.redraw()

    def redraw(self):
        img = self.get_current_image()

        self.canvas.delete("all")

        if img is None:
            return

        cw = max(
            1,
            self.canvas.winfo_width()
        )
        ch = max(
            1,
            self.canvas.winfo_height()
        )

        ih, iw = img.shape[:2]

        # visible image coordinates
        x1 = int(
            np.floor(
                (0 - self.offset_x)
                / max(self.zoom, 1e-9)
            )
        )
        y1 = int(
            np.floor(
                (0 - self.offset_y)
                / max(self.zoom, 1e-9)
            )
        )
        x2 = int(
            np.ceil(
                (cw - self.offset_x)
                / max(self.zoom, 1e-9)
            )
        )
        y2 = int(
            np.ceil(
                (ch - self.offset_y)
                / max(self.zoom, 1e-9)
            )
        )

        x1 = max(0, x1)
        y1 = max(0, y1)
        x2 = min(iw, x2)
        y2 = min(ih, y2)

        if x2 <= x1 or y2 <= y1:
            return

        crop = img[
            y1:y2,
            x1:x2
        ]

        dw = max(
            1,
            int(round(
                crop.shape[1]
                * self.zoom
            ))
        )
        dh = max(
            1,
            int(round(
                crop.shape[0]
                * self.zoom
            ))
        )

        resized = cv2.resize(
            crop,
            (dw, dh),
            interpolation=cv2.INTER_NEAREST
        )

        self.tk_img = ImageTk.PhotoImage(
            Image.fromarray(resized)
        )

        canvas_x = (
            self.offset_x
            + x1 * self.zoom
        )
        canvas_y = (
            self.offset_y
            + y1 * self.zoom
        )

        self.canvas.create_image(
            canvas_x,
            canvas_y,
            anchor=tk.NW,
            image=self.tk_img
        )

    def canvas_to_image_xy(
        self,
        sx,
        sy
    ):
        if self.zoom <= 0:
            return None

        ix = int(
            np.floor(
                (float(sx) - self.offset_x)
                / self.zoom
            )
        )
        iy = int(
            np.floor(
                (float(sy) - self.offset_y)
                / self.zoom
            )
        )

        if (
            ix < 0
            or iy < 0
            or ix >= self.w
            or iy >= self.h
        ):
            return None

        return ix, iy

    @staticmethod
    def _fmt_mouse_value(v):
        try:
            fv = float(v)
        except Exception:
            return str(v)

        if not np.isfinite(fv):
            return "nan"

        # Keep enough precision for raw fluorescence values while avoiding noise.
        if abs(fv) >= 100:
            return f"{fv:.1f}"
        if abs(fv) >= 10:
            return f"{fv:.2f}"
        return f"{fv:.3f}"

    def on_mouse_move(self, event):
        pos = self.canvas_to_image_xy(
            event.x,
            event.y
        )

        if pos is None:
            self.mouse_value_label.config(
                text="Mouse: 画像外"
            )
            return

        x, y = pos

        if not self.qc_arrays:
            self.mouse_value_label.config(
                text=f"Mouse: x={x}, y={y}  /  未計算"
            )
            return

        q = self.qc_arrays

        def val(name):
            arr = q.get(name)
            if arr is None:
                return None
            try:
                return arr[y, x]
            except Exception:
                return None

        raw = val("raw_zo1")
        bg = val("background")
        sub = val("background_subtracted")
        top = val("white_tophat")
        binary = val("binary")
        skel = val("skeleton")
        dapi = val("dapi")

        lines = [
            f"Mouse x={x}, y={y}",
            (
                "ZO-1  Raw="
                f"{self._fmt_mouse_value(raw)}"
                "   Background="
                f"{self._fmt_mouse_value(bg)}"
            ),
            (
                "BG-sub="
                f"{self._fmt_mouse_value(sub)}"
                "   White top-hat="
                f"{self._fmt_mouse_value(top)}"
            ),
            (
                "Binary="
                f"{1 if bool(binary) else 0}"
                "   Skeleton="
                f"{1 if bool(skel) else 0}"
            ),
        ]

        if dapi is not None:
            lines.append(
                "DAPI="
                f"{self._fmt_mouse_value(dapi)}"
            )

        self.mouse_value_label.config(
            text="\n".join(lines)
        )

    def on_mouse_leave(self, event):
        self.mouse_value_label.config(
            text="Mouse: 画像上にカーソルを置いてください"
        )

    def on_wheel(self, event):
        factor = (
            1.15
            if event.delta > 0
            else 1 / 1.15
        )
        self.zoom_at(
            event.x,
            event.y,
            factor
        )

    def on_wheel_linux(self, event):
        factor = (
            1.15
            if event.num == 4
            else 1 / 1.15
        )
        self.zoom_at(
            event.x,
            event.y,
            factor
        )

    def zoom_at(
        self,
        sx,
        sy,
        factor
    ):
        old_zoom = self.zoom

        new_zoom = float(
            np.clip(
                old_zoom * factor,
                0.05,
                30.0
            )
        )

        scale = (
            new_zoom
            / max(old_zoom, 1e-9)
        )

        self.offset_x = (
            sx
            - scale
            * (sx - self.offset_x)
        )

        self.offset_y = (
            sy
            - scale
            * (sy - self.offset_y)
        )

        self.zoom = new_zoom
        self.redraw()

    def on_right_down(self, event):
        self.last_drag = (
            event.x,
            event.y
        )

    def on_right_drag(self, event):
        if self.last_drag is None:
            return

        dx = (
            event.x
            - self.last_drag[0]
        )
        dy = (
            event.y
            - self.last_drag[1]
        )

        self.offset_x += dx
        self.offset_y += dy

        self.last_drag = (
            event.x,
            event.y
        )

        self.redraw()

    def on_right_up(self, event):
        self.last_drag = None

    # --------------------------------------------------------
    # Apply tuned params
    # --------------------------------------------------------

    def apply_to_main(self):
        try:
            p = self.get_params()
        except Exception as exc:
            messagebox.showerror(
                "ZO-1 QC",
                str(exc),
                parent=self
            )
            return

        if self.apply_callback is not None:
            self.apply_callback(p)

        messagebox.showinfo(
            "ZO-1 QC",
            (
                "現在のZO-1前処理設定を"
                "本解析パラメータへ反映しました。"
            ),
            parent=self
        )


# ============================================================
# Main sample browser
# ============================================================

class ZO1AnalysisTab(ttk.Frame):
    def __init__(self, parent, host_app=None):
        super().__init__(parent)
        self.host_app = host_app

        self.root_folder = None
        self.samples = []
        self.sample_groups = {}
        self.group_name_var = tk.StringVar(value="")

        self.global_settings_path = (
            Path.home()
            / GLOBAL_SETTINGS_FILENAME
        )

        self._global_settings_existed_at_start = (
            self.global_settings_path.exists()
        )
        self._legacy_settings_migration_checked = False
        self._analysis_settings_after_id = None

        try:
            self.settings = (
                load_json(
                    self.global_settings_path
                )
                if self.global_settings_path.exists()
                else {}
            )
        except Exception:
            self.settings = {}

        self.editor = None

        self.analysis_queue = queue.Queue()
        self.analysis_running = False

        # Batch StarDist runs sequentially in one worker thread so the model is
        # loaded only once and RAM usage stays bounded.
        self.batch_stardist_queue = queue.Queue()
        self.batch_stardist_running = False

        _as = self.settings.get(
            "__zo1_analysis__",
            {}
        )

        self.dapi_postfix_var = tk.StringVar(
            value=str(_as.get("dapi_postfix", "_DAPI"))
        )
        self.zo1_postfix_var = tk.StringVar(
            value=str(_as.get("zo1_postfix", "_ZO1"))
        )

        self.junction_lanes_var = tk.IntVar(
            value=int(
                _as.get(
                    "junction_lanes",
                    DEFAULT_JUNCTION_LANES
                )
            )
        )
        self.pair_distance_factor_var = tk.DoubleVar(
            value=float(
                _as.get(
                    "pair_distance_factor",
                    DEFAULT_PAIR_DISTANCE_FACTOR
                )
            )
        )
        self.zo1_rolling_ball_radius_var = tk.IntVar(
            value=int(
                _as.get(
                    "zo1_rolling_ball_radius_px",
                    DEFAULT_ZO1_ROLLING_BALL_RADIUS_PX
                )
            )
        )

        self.zo1_tophat_radius_var = tk.IntVar(
            value=int(
                _as.get(
                    "zo1_tophat_radius_px",
                    DEFAULT_ZO1_TOPHAT_RADIUS_PX
                )
            )
        )

        self.zo1_tophat_threshold_var = tk.DoubleVar(
            value=float(
                _as.get(
                    "zo1_tophat_threshold",
                    DEFAULT_ZO1_TOPHAT_THRESHOLD
                )
            )
        )

        self.zo1_min_object_size_var = tk.IntVar(
            value=int(
                _as.get(
                    "zo1_min_object_size_px",
                    DEFAULT_ZO1_MIN_OBJECT_SIZE_PX
                )
            )
        )

        self.intact_cutoff_var = tk.DoubleVar(
            value=float(
                _as.get(
                    "intact_coverage_threshold_pct",
                    DEFAULT_INTACT_COVERAGE_THRESHOLD_PCT
                )
            )
        )

        self.lane_peak_top_n_var = tk.IntVar(
            value=int(
                _as.get(
                    "lane_peak_top_n_pixels",
                    DEFAULT_LANE_PEAK_TOP_N_PIXELS
                )
            )
        )

        self.build_ui()

        for _var in (
            self.dapi_postfix_var,
            self.zo1_postfix_var,
            self.junction_lanes_var,
            self.pair_distance_factor_var,
            self.zo1_rolling_ball_radius_var,
            self.zo1_tophat_radius_var,
            self.zo1_tophat_threshold_var,
            self.zo1_min_object_size_var,
            self.intact_cutoff_var,
            self.lane_peak_top_n_var,
        ):
            _var.trace_add(
                "write",
                self._schedule_analysis_settings_save
            )


    def build_ui(self):
        top = ttk.Frame(self)
        top.pack(
            fill=tk.X,
            padx=10,
            pady=10
        )

        ttk.Button(
            top,
            text="親フォルダを選択",
            command=self.choose_parent_folder
        ).pack(side=tk.LEFT)

        self.parent_label = ttk.Label(
            top,
            text="親フォルダ未選択"
        )

        self.parent_label.pack(
            side=tk.LEFT,
            fill=tk.X,
            expand=True,
            padx=10
        )

        postfix_box = ttk.LabelFrame(
            self,
            text="入力ファイル postfix / Channel postfix",
            padding=(10, 6)
        )
        postfix_box.pack(
            fill=tk.X,
            padx=10,
            pady=(0, 6)
        )

        ttk.Label(postfix_box, text="DAPI postfix:").pack(
            side=tk.LEFT
        )
        dapi_entry = ttk.Entry(
            postfix_box,
            textvariable=self.dapi_postfix_var,
            width=18
        )
        dapi_entry.pack(side=tk.LEFT, padx=(5, 14))

        ttk.Label(postfix_box, text="ZO-1 postfix:").pack(
            side=tk.LEFT
        )
        zo1_entry = ttk.Entry(
            postfix_box,
            textvariable=self.zo1_postfix_var,
            width=18
        )
        zo1_entry.pack(side=tk.LEFT, padx=(5, 14))

        ttk.Button(
            postfix_box,
            text="再スキャン / Rescan",
            command=self.rescan_with_current_postfixes
        ).pack(side=tk.LEFT, padx=(0, 8))

        if self.host_app is not None:
            ttk.Button(
                postfix_box,
                text="MAX設定から取得",
                command=self.copy_postfixes_from_max_tab
            ).pack(side=tk.LEFT)

        for ent in (dapi_entry, zo1_entry):
            ent.bind(
                "<Return>",
                lambda _e: self.rescan_with_current_postfixes()
            )

        self.postfix_hint_label = ttk.Label(
            self,
            text=""
        )
        self.postfix_hint_label.pack(
            anchor="w",
            padx=12,
            pady=(0, 5)
        )
        self._update_postfix_hint()

        frame = ttk.Frame(self)

        frame.pack(
            fill=tk.BOTH,
            expand=True,
            padx=10,
            pady=5
        )

        columns = (
            "no",
            "sample",
            "group",
            "dapi",
            "zo1",
            "status",
            "nuclei"
        )

        self.tree = ttk.Treeview(
            frame,
            columns=columns,
            show="headings",
            selectmode="extended"
        )

        config = {
            "no": ("#", 55),
            "sample": ("Sample", 150),
            "group": ("Group", 120),
            "dapi": ("DAPI", 210),
            "zo1": ("ZO-1", 210),
            "status": ("核アノテーション", 190),
            "nuclei": ("核数", 80),
        }

        for col in columns:
            text, width = config[col]

            self.tree.heading(
                col,
                text=text
            )

            self.tree.column(
                col,
                width=width,
                anchor=(
                    "center"
                    if col in {"no", "nuclei"}
                    else "w"
                )
            )

        ys = ttk.Scrollbar(
            frame,
            orient=tk.VERTICAL,
            command=self.tree.yview
        )

        xs = ttk.Scrollbar(
            frame,
            orient=tk.HORIZONTAL,
            command=self.tree.xview
        )

        self.tree.configure(
            yscrollcommand=ys.set,
            xscrollcommand=xs.set
        )

        self.tree.grid(
            row=0,
            column=0,
            sticky="nsew"
        )

        ys.grid(
            row=0,
            column=1,
            sticky="ns"
        )

        xs.grid(
            row=1,
            column=0,
            sticky="ew"
        )

        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)

        self.tree.bind(
            "<Double-1>",
            lambda e: self.open_selected()
        )

        buttons = ttk.Frame(self)

        buttons.pack(
            fill=tk.X,
            padx=10,
            pady=(5, 10)
        )

        ttk.Button(
            buttons,
            text="選択サンプルを開く",
            command=self.open_selected
        ).pack(
            side=tk.LEFT,
            padx=3
        )

        ttk.Button(
            buttons,
            text="最初の未保存サンプルを開く",
            command=self.open_first_unfinished
        ).pack(
            side=tk.LEFT,
            padx=3
        )

        ttk.Button(
            buttons,
            text="一覧を再読み込み",
            command=self.scan_samples
        ).pack(
            side=tk.LEFT,
            padx=3
        )

        self.batch_stardist_button = ttk.Button(
            buttons,
            text="全未確認を一括StarDist",
            command=self.run_batch_stardist
        )
        self.batch_stardist_button.pack(
            side=tk.LEFT,
            padx=(12, 3)
        )

        group_frame = ttk.Frame(self)
        group_frame.pack(
            fill=tk.X,
            padx=10,
            pady=(0, 6)
        )

        ttk.Label(
            group_frame,
            text="Group:"
        ).pack(
            side=tk.LEFT,
            padx=(3, 4)
        )

        ttk.Entry(
            group_frame,
            textvariable=self.group_name_var,
            width=18
        ).pack(
            side=tk.LEFT,
            padx=(0, 5)
        )

        ttk.Button(
            group_frame,
            text="選択サンプルにGroupを設定",
            command=self.assign_group_to_selected
        ).pack(
            side=tk.LEFT,
            padx=3
        )

        ttk.Button(
            group_frame,
            text="選択サンプルのGroupをクリア",
            command=self.clear_group_from_selected
        ).pack(
            side=tk.LEFT,
            padx=3
        )

        ttk.Label(
            group_frame,
            text=(
                "複数選択可（Ctrl/Shift）。"
                "Groupはsummary CSVに保存されます。"
            )
        ).pack(
            side=tk.LEFT,
            padx=12
        )

        self.status = ttk.Label(
            buttons,
            text=""
        )

        self.status.pack(
            side=tk.RIGHT,
            fill=tk.X,
            expand=True,
            padx=10
        )

        analysis_box = ttk.LabelFrame(
            self,
            text="Tight-junction analysis"
        )

        analysis_box.pack(
            fill=tk.X,
            padx=10,
            pady=(0, 10)
        )

        ttk.Label(
            analysis_box,
            text="台形分割数"
        ).grid(
            row=0,
            column=0,
            sticky="w",
            padx=5,
            pady=3
        )

        ttk.Entry(
            analysis_box,
            textvariable=self.junction_lanes_var,
            width=8
        ).grid(
            row=0,
            column=1,
            sticky="w",
            padx=5,
            pady=3
        )

        ttk.Label(
            analysis_box,
            text="Max pair distance factor"
        ).grid(
            row=0,
            column=2,
            sticky="w",
            padx=(14, 5),
            pady=3
        )

        ttk.Entry(
            analysis_box,
            textvariable=self.pair_distance_factor_var,
            width=8
        ).grid(
            row=0,
            column=3,
            sticky="w",
            padx=5,
            pady=3
        )

        ttk.Label(
            analysis_box,
            text="ZO-1 rolling ball radius (px)"
        ).grid(
            row=0,
            column=4,
            sticky="w",
            padx=(14, 5),
            pady=3
        )

        ttk.Entry(
            analysis_box,
            textvariable=self.zo1_rolling_ball_radius_var,
            width=8
        ).grid(
            row=0,
            column=5,
            sticky="w",
            padx=5,
            pady=3
        )

        ttk.Label(
            analysis_box,
            text="White top-hat radius (px)"
        ).grid(
            row=0,
            column=6,
            sticky="w",
            padx=(14, 5),
            pady=3
        )

        ttk.Entry(
            analysis_box,
            textvariable=self.zo1_tophat_radius_var,
            width=8
        ).grid(
            row=0,
            column=7,
            sticky="w",
            padx=5,
            pady=3
        )

        ttk.Label(
            analysis_box,
            text="Top-hat threshold"
        ).grid(
            row=1,
            column=0,
            sticky="w",
            padx=5,
            pady=3
        )

        ttk.Entry(
            analysis_box,
            textvariable=self.zo1_tophat_threshold_var,
            width=8
        ).grid(
            row=1,
            column=1,
            sticky="w",
            padx=5,
            pady=3
        )

        ttk.Label(
            analysis_box,
            text="Min object size (px)"
        ).grid(
            row=1,
            column=2,
            sticky="w",
            padx=(14, 5),
            pady=3
        )

        ttk.Entry(
            analysis_box,
            textvariable=self.zo1_min_object_size_var,
            width=8
        ).grid(
            row=1,
            column=3,
            sticky="w",
            padx=5,
            pady=3
        )

        ttk.Label(
            analysis_box,
            text="Intact cutoff (%)"
        ).grid(
            row=1,
            column=4,
            sticky="w",
            padx=(14, 5),
            pady=3
        )

        ttk.Entry(
            analysis_box,
            textvariable=self.intact_cutoff_var,
            width=8
        ).grid(
            row=1,
            column=5,
            sticky="w",
            padx=5,
            pady=3
        )

        ttk.Label(
            analysis_box,
            text="Lane peak top N pixels"
        ).grid(
            row=1,
            column=6,
            sticky="w",
            padx=(14, 5),
            pady=3
        )

        ttk.Entry(
            analysis_box,
            textvariable=self.lane_peak_top_n_var,
            width=8
        ).grid(
            row=1,
            column=7,
            sticky="w",
            padx=5,
            pady=3
        )

        self.qc_button = ttk.Button(
            analysis_box,
            text="代表画像でZO-1処理を確認・調整",
            command=self.open_zo1_qc
        )

        self.qc_button.grid(
            row=2,
            column=0,
            columnspan=3,
            sticky="ew",
            padx=5,
            pady=(5, 6)
        )

        self.analyze_selected_button = ttk.Button(
            analysis_box,
            text="選択サンプルを解析",
            command=self.analyze_selected_sample
        )

        self.analyze_selected_button.grid(
            row=3,
            column=0,
            columnspan=2,
            sticky="ew",
            padx=5,
            pady=(5, 6)
        )

        self.analyze_all_button = ttk.Button(
            analysis_box,
            text="全保存済みを解析＋Group summary",
            command=self.analyze_all_saved_samples
        )

        self.analyze_all_button.grid(
            row=3,
            column=2,
            columnspan=3,
            sticky="ew",
            padx=5,
            pady=(5, 6)
        )

        self.analysis_status = ttk.Label(
            analysis_box,
            text="解析待機",
            anchor="w"
        )

        self.analysis_status.grid(
            row=3,
            column=5,
            columnspan=3,
            sticky="ew",
            padx=8,
            pady=(5, 6)
        )

        for c in range(8):
            analysis_box.columnconfigure(
                c,
                weight=1 if c in (0, 2, 4, 6) else 0
            )

    # ---------------- folder/list ----------------

    def settings_path(self):
        # ONE shared settings file for every dataset/folder.
        return self.global_settings_path

    def _save_global_settings(self):
        try:
            save_json(
                self.global_settings_path,
                self.settings
            )
        except Exception:
            pass

    def _apply_analysis_settings_dict(
        self,
        saved
    ):
        if not isinstance(saved, dict):
            return

        mapping = (
            (
                self.dapi_postfix_var,
                "dapi_postfix",
                str,
            ),
            (
                self.zo1_postfix_var,
                "zo1_postfix",
                str,
            ),
            (
                self.junction_lanes_var,
                "junction_lanes",
                int,
            ),
            (
                self.pair_distance_factor_var,
                "pair_distance_factor",
                float,
            ),
            (
                self.zo1_rolling_ball_radius_var,
                "zo1_rolling_ball_radius_px",
                int,
            ),
            (
                self.zo1_tophat_radius_var,
                "zo1_tophat_radius_px",
                int,
            ),
            (
                self.zo1_tophat_threshold_var,
                "zo1_tophat_threshold",
                float,
            ),
            (
                self.zo1_min_object_size_var,
                "zo1_min_object_size_px",
                int,
            ),
            (
                self.intact_cutoff_var,
                "intact_coverage_threshold_pct",
                float,
            ),
            (
                self.lane_peak_top_n_var,
                "lane_peak_top_n_pixels",
                int,
            ),
        )

        for var, key, caster in mapping:
            if key not in saved:
                continue

            try:
                var.set(
                    caster(saved[key])
                )
            except Exception:
                pass

    def _schedule_analysis_settings_save(
        self,
        *_args
    ):
        if self._analysis_settings_after_id is not None:
            try:
                self.after_cancel(
                    self._analysis_settings_after_id
                )
            except Exception:
                pass

        self._analysis_settings_after_id = self.after(
            500,
            self._persist_analysis_settings_from_vars
        )

    def _persist_analysis_settings_from_vars(self):
        self._analysis_settings_after_id = None

        try:
            vals = {
                "dapi_postfix":
                    str(self.dapi_postfix_var.get()).strip(),
                "zo1_postfix":
                    str(self.zo1_postfix_var.get()).strip(),
                "junction_lanes":
                    int(
                        self.junction_lanes_var.get()
                    ),
                "pair_distance_factor":
                    float(
                        self.pair_distance_factor_var.get()
                    ),
                "zo1_rolling_ball_radius_px":
                    int(
                        self.zo1_rolling_ball_radius_var.get()
                    ),
                "zo1_tophat_radius_px":
                    int(
                        self.zo1_tophat_radius_var.get()
                    ),
                "zo1_tophat_threshold":
                    float(
                        self.zo1_tophat_threshold_var.get()
                    ),
                "zo1_min_object_size_px":
                    int(
                        self.zo1_min_object_size_var.get()
                    ),
                "intact_coverage_threshold_pct":
                    float(
                        self.intact_cutoff_var.get()
                    ),
                "lane_peak_top_n_pixels":
                    int(
                        self.lane_peak_top_n_var.get()
                    ),
            }
        except Exception:
            # User may be midway through editing an Entry.
            return

        self.settings[
            "__zo1_analysis__"
        ] = vals

        self._save_global_settings()

    def persist_stardist_settings(
        self,
        params
    ):
        self.settings[
            "__stardist__"
        ] = {
            k: float(params[k])
            for k in (
                "pixel_size_um",
                "prob_thresh",
                "nms_thresh",
                "min_area_um2",
                "max_area_um2",
                "normalize_pmin",
                "normalize_pmax",
            )
        }

        self._save_global_settings()

    def _maybe_migrate_legacy_parent_settings(
        self
    ):
        """
        First run of this version only:
        if no global settings file existed, import the old per-parent JSON once.
        Later folder changes never reload settings.
        """
        if self._legacy_settings_migration_checked:
            return

        self._legacy_settings_migration_checked = True

        if self._global_settings_existed_at_start:
            return

        if self.root_folder is None:
            return

        legacy_path = (
            self.root_folder
            / SETTINGS_FILENAME
        )

        if not legacy_path.exists():
            self._persist_analysis_settings_from_vars()
            return

        try:
            legacy = load_json(
                legacy_path
            )
        except Exception:
            legacy = {}

        if isinstance(legacy, dict):
            # Keep current structure but import old StarDist/ZO-1 settings.
            for key in (
                "__stardist__",
                "__zo1_analysis__",
            ):
                if key in legacy:
                    self.settings[key] = legacy[key]

            self._apply_analysis_settings_dict(
                self.settings.get(
                    "__zo1_analysis__",
                    {}
                )
            )

        self._save_global_settings()

    def on_main_close(self):
        self._persist_analysis_settings_from_vars()

    def error_log_path(self):
        if self.root_folder is None:
            return None

        return (
            self.root_folder
            / ERROR_LOG_FILENAME
        )

    def choose_parent_folder(self):
        folder = filedialog.askdirectory(
            title="親フォルダを選択"
        )

        if not folder:
            return

        self.set_parent_folder(folder)

    def set_parent_folder(self, folder):
        """Load a folder directly, including hand-off from the Mica MAX tab."""
        if not folder:
            return

        self.root_folder = Path(folder)

        self.parent_label.config(
            text=str(self.root_folder)
        )

        # IMPORTANT:
        # changing parent folders does NOT reload or replace current parameters.
        # Only on the first run of this version can an old per-parent settings
        # JSON be migrated once into the new global settings file.
        self._maybe_migrate_legacy_parent_settings()

        self.scan_samples()

    def sample_groups_path(self):
        if self.root_folder is None:
            return None

        return (
            self.root_folder
            / SAMPLE_GROUPS_FILENAME
        )

    def _load_sample_groups(self):
        p = self.sample_groups_path()

        if p is None or not p.exists():
            self.sample_groups = {}
            return

        x = load_json(p)

        self.sample_groups = {
            str(k): str(v)
            for k, v in x.items()
        }

    def _save_sample_groups(self):
        p = self.sample_groups_path()

        if p is None:
            return

        save_json(
            p,
            self.sample_groups
        )

    def selected_indices(self):
        out = []

        for iid in self.tree.selection():
            try:
                out.append(int(iid))
            except Exception:
                pass

        return sorted(set(out))

    def assign_group_to_selected(self):
        indices = self.selected_indices()

        if not indices:
            messagebox.showinfo(
                "Group",
                "Groupを設定するサンプルを選択してください。",
                parent=self
            )
            return

        group = self.group_name_var.get().strip()

        if not group:
            messagebox.showwarning(
                "Group",
                "Group名を入力してください。",
                parent=self
            )
            return

        selected_names = []

        for idx in indices:
            if 0 <= idx < len(self.samples):
                name = sample_display_name(
                    self.samples[idx]
                )
                self.sample_groups[name] = group
                selected_names.append(name)

        self._save_sample_groups()
        self.scan_samples(
            preserve_selection=selected_names
        )

    def clear_group_from_selected(self):
        indices = self.selected_indices()

        if not indices:
            messagebox.showinfo(
                "Group",
                "Groupをクリアするサンプルを選択してください。",
                parent=self
            )
            return

        selected_names = []

        for idx in indices:
            if 0 <= idx < len(self.samples):
                name = sample_display_name(
                    self.samples[idx]
                )
                self.sample_groups.pop(
                    name,
                    None
                )
                selected_names.append(name)

        self._save_sample_groups()
        self.scan_samples(
            preserve_selection=selected_names
        )

    def _current_channel_postfixes(self):
        dapi = str(self.dapi_postfix_var.get()).strip()
        zo1 = str(self.zo1_postfix_var.get()).strip()
        return dapi, zo1

    def _update_postfix_hint(self):
        if not hasattr(self, "postfix_hint_label"):
            return
        dapi, zo1 = self._current_channel_postfixes()
        self.postfix_hint_label.configure(
            text=(
                f"自動ペアリング: <sample>{dapi}.tif(f) / "
                f"<sample>{zo1}.tif(f)"
            )
        )

    def rescan_with_current_postfixes(self):
        dapi, zo1 = self._current_channel_postfixes()

        if not dapi or not zo1:
            messagebox.showwarning(
                "Postfix",
                "DAPI postfix と ZO-1 postfix の両方を入力してください。",
                parent=self
            )
            return

        if dapi.casefold() == zo1.casefold():
            messagebox.showwarning(
                "Postfix",
                "DAPI postfix と ZO-1 postfix は異なる文字列にしてください。",
                parent=self
            )
            return

        self._update_postfix_hint()
        self._persist_analysis_settings_from_vars()

        if self.root_folder is not None:
            self.scan_samples()

    def copy_postfixes_from_max_tab(self):
        host = self.host_app
        mica_tab = getattr(host, "mica_tab", None) if host is not None else None

        if mica_tab is None:
            return

        try:
            n_channels = int(mica_tab.n_channels.get())
        except Exception:
            n_channels = 2

        if n_channels < 2:
            return

        self.dapi_postfix_var.set(
            str(mica_tab.postfix_vars[0].get()).strip()
        )
        self.zo1_postfix_var.set(
            str(mica_tab.postfix_vars[1].get()).strip()
        )
        self.rescan_with_current_postfixes()

    def scan_samples(
        self,
        preserve_selection=None
    ):
        """
        Pair TIFF files in ROOT using the user-defined DAPI / ZO-1 postfixes.

        Example postfixes:
            DAPI: _DAPI_MAX
            ZO-1: _ZO1_MAX

        Then:
            con1_DAPI_MAX.tiff + con1_ZO1_MAX.tiff -> sample 'con1'
        """
        if self.root_folder is None:
            return

        dapi_postfix, zo1_postfix = self._current_channel_postfixes()

        if (
            not dapi_postfix
            or not zo1_postfix
            or dapi_postfix.casefold() == zo1_postfix.casefold()
        ):
            return

        self._update_postfix_hint()
        self._load_sample_groups()

        paired = {}

        files = sorted(
            [
                p
                for p in self.root_folder.iterdir()
                if (
                    p.is_file()
                    and p.suffix.lower()
                    in {".tif", ".tiff"}
                )
            ],
            key=lambda p: natural_sort_key(
                p.name
            )
        )

        for p in files:
            matched = split_sample_channel_by_postfix(
                p.stem,
                dapi_postfix,
                zo1_postfix
            )

            if matched is None:
                continue

            sample_name, channel = matched
            key = sample_name.casefold()

            if key not in paired:
                paired[key] = {
                    "name": sample_name,
                    "dapi": None,
                    "zo1": None,
                }

            # If capitalization differs between files,
            # preserve the name from the first encountered file.
            if channel == "DAPI":
                paired[key]["dapi"] = p
            elif channel == "ZO1":
                paired[key]["zo1"] = p

        samples = []

        for key in sorted(
            paired,
            key=lambda k: natural_sort_key(
                paired[k]["name"]
            )
        ):
            item = paired[key]
            name = item["name"]

            sample = {
                "folder": self.root_folder,
                "name": name,
                "group": self.sample_groups.get(
                    name,
                    ""
                ),
                "flat_input": True,
                "dapi": item["dapi"],
                "zo1": item["zo1"],
            }

            paths = sample_output_paths(
                sample
            )

            saved = paths["npz"].exists()
            draft = (
                paths["draft_npz"].exists()
                and not saved
            )

            count = None

            count_path = (
                paths["npz"]
                if saved
                else (
                    paths["draft_npz"]
                    if draft
                    else None
                )
            )

            if count_path is not None:
                try:
                    d = np.load(
                        count_path,
                        allow_pickle=True
                    )

                    if "nuclei_labels" in d:
                        labels = d["nuclei_labels"]

                        if find_objects is not None:
                            count = sum(
                                s is not None
                                for s in find_objects(
                                    labels
                                )
                            )
                        else:
                            count = int(
                                len(
                                    np.unique(labels)
                                )
                                - 1
                            )
                except Exception:
                    count = None

            dapi = sample["dapi"]
            zo1 = sample["zo1"]

            if dapi is None or zo1 is None:
                missing = []

                if dapi is None:
                    missing.append("DAPI")

                if zo1 is None:
                    missing.append("ZO-1")

                status = (
                    "不足: "
                    + ", ".join(missing)
                )

            else:
                if saved:
                    status = "確認済み・保存済み"
                elif draft:
                    status = "StarDist済み・未確認"
                else:
                    status = "未処理"

            sample.update({
                "saved": saved,
                "draft": draft,
                "nuclei": count,
                "status": status,
            })

            samples.append(sample)

        self.samples = samples
        self.refresh_tree()

        if preserve_selection:
            wanted = {
                str(x)
                for x in preserve_selection
            }

            for i, s in enumerate(
                self.samples
            ):
                if sample_display_name(s) in wanted:
                    self.tree.selection_add(
                        str(i)
                    )

    def refresh_tree(self):
        for iid in self.tree.get_children():
            self.tree.delete(iid)

        ready = 0
        saved = 0
        draft = 0

        for i, s in enumerate(
            self.samples,
            start=1
        ):
            if (
                s["dapi"] is not None
                and s["zo1"] is not None
            ):
                ready += 1

                if s["saved"]:
                    saved += 1
                elif s.get("draft", False):
                    draft += 1

            self.tree.insert(
                "",
                "end",
                iid=str(i - 1),
                values=(
                    i,
                    sample_display_name(s),
                    s.get("group", ""),
                    (
                        s["dapi"].name
                        if s["dapi"]
                        else "NOT FOUND"
                    ),
                    (
                        s["zo1"].name
                        if s["zo1"]
                        else "NOT FOUND"
                    ),
                    s["status"],
                    (
                        s["nuclei"]
                        if s["nuclei"] is not None
                        else ""
                    ),
                )
            )

        self.status.config(
            text=(
                f"解析可能 {ready} / "
                f"確認済み {saved} / "
                f"StarDist未確認 {draft}"
            )
        )

    # ---------------- batch StarDist ----------------

    def get_batch_stardist_params(self):
        """
        Use the latest accepted StarDist parameters.
        If none have been saved yet, use the same defaults as the editor.
        """
        s = self.settings.get(
            "__stardist__",
            {}
        )

        p = {
            "pixel_size_um": float(
                s.get(
                    "pixel_size_um",
                    DEFAULT_PIXEL_SIZE_UM
                )
            ),
            "prob_thresh": float(
                s.get(
                    "prob_thresh",
                    DEFAULT_PROB_THRESH
                )
            ),
            "nms_thresh": float(
                s.get(
                    "nms_thresh",
                    DEFAULT_NMS_THRESH
                )
            ),
            "min_area_um2": float(
                s.get(
                    "min_area_um2",
                    DEFAULT_MIN_AREA_UM2
                )
            ),
            "max_area_um2": float(
                s.get(
                    "max_area_um2",
                    DEFAULT_MAX_AREA_UM2
                )
            ),
            "normalize_pmin": float(
                s.get(
                    "normalize_pmin",
                    DEFAULT_NORMALIZE_PMIN
                )
            ),
            "normalize_pmax": float(
                s.get(
                    "normalize_pmax",
                    DEFAULT_NORMALIZE_PMAX
                )
            ),
        }

        return p

    def run_batch_stardist(self):
        if self.batch_stardist_running:
            return

        if self.root_folder is None:
            messagebox.showinfo(
                "Batch StarDist",
                "先に親フォルダを選択してください。",
                parent=self
            )
            return

        # Never overwrite manually reviewed/saved annotations.
        # Existing unreviewed drafts ARE recalculated, so changing StarDist
        # parameters and rerunning the batch is straightforward.
        targets = [
            s
            for s in self.samples
            if (
                s["dapi"] is not None
                and not s["saved"]
            )
        ]

        if not targets:
            messagebox.showinfo(
                "Batch StarDist",
                "未確認のDAPI画像はありません。",
                parent=self
            )
            return

        params = self.get_batch_stardist_params()
        self.persist_stardist_settings(
            params
        )

        n_existing_draft = sum(
            1
            for s in targets
            if s.get("draft", False)
        )

        msg = (
            f"未確認DAPI {len(targets)}枚を連続StarDistします。\n\n"
            f"確認済み・保存済み画像は上書きしません。\n"
            f"既存の未確認StarDistドラフトは再計算します: "
            f"{n_existing_draft}枚\n\n"
            "使用設定:\n"
            f"prob={params['prob_thresh']:.3f}, "
            f"NMS={params['nms_thresh']:.3f}\n"
            f"area={params['min_area_um2']:.1f}"
            f"–{params['max_area_um2']:.1f} µm²\n"
            f"normalize="
            f"{params['normalize_pmin']:.1f}"
            f"–{params['normalize_pmax']:.1f} percentile\n"
            "接触核の8近傍Auto-Merge: ON\n\n"
            "開始しますか？"
        )

        if not messagebox.askyesno(
            "Batch StarDist",
            msg,
            parent=self
        ):
            return

        self.batch_stardist_running = True
        self.batch_stardist_button.config(
            state=tk.DISABLED
        )

        self.status.config(
            text=(
                f"Batch StarDist開始: "
                f"0/{len(targets)}"
            )
        )

        threading.Thread(
            target=self._batch_stardist_worker,
            args=(targets, params),
            daemon=True
        ).start()

        self.after(
            100,
            self.poll_batch_stardist_queue
        )

    def _batch_stardist_worker(
        self,
        targets,
        params
    ):
        success = []
        errors = []

        # get_model() is lazy and cached, so one worker processes sequentially
        # without loading the model repeatedly.
        try:
            StarDistNucleiSegmenter.get_model()
        except Exception:
            self.batch_stardist_queue.put({
                "type": "fatal",
                "error": traceback.format_exc()
            })
            return

        total = len(targets)

        for i, sample in enumerate(
            targets,
            start=1
        ):
            folder = Path(sample["folder"])
            display_name = sample_display_name(
                sample
            )
            dapi_path = Path(sample["dapi"])
            zo1_path = sample.get("zo1")
            paths = sample_output_paths(
                sample
            )

            self.batch_stardist_queue.put({
                "type": "progress",
                "index": i,
                "total": total,
                "sample": display_name,
                "phase": "reading",
            })

            try:
                raw = read_tiff_gray(
                    dapi_path
                )

                self.batch_stardist_queue.put({
                    "type": "progress",
                    "index": i,
                    "total": total,
                    "sample": display_name,
                    "phase": "StarDist",
                })

                labels = (
                    StarDistNucleiSegmenter.segment(
                        raw,
                        **params
                    )
                )

                self.batch_stardist_queue.put({
                    "type": "progress",
                    "index": i,
                    "total": total,
                    "sample": display_name,
                    "phase": "merge touching nuclei",
                })

                (
                    labels,
                    merge_groups,
                    merge_removed,
                ) = merge_touching_labels_8neighbor(
                    labels
                )

                save_stardist_draft(
                    paths["draft_npz"],
                    labels=labels,
                    sample_name=display_name,
                    dapi_name=dapi_path.name,
                    zo1_name=(
                        Path(zo1_path).name
                        if zo1_path is not None
                        else ""
                    ),
                    params=params,
                    touching_merge_groups=merge_groups,
                    touching_labels_removed=merge_removed,
                )

                # Auto-merge helper relabels surviving nuclei contiguously.
                n_nuclei = int(
                    labels.max()
                )

                success.append({
                    "sample": display_name,
                    "nuclei": n_nuclei,
                    "touching_merge_groups": merge_groups,
                    "touching_labels_removed": merge_removed,
                })

                self.batch_stardist_queue.put({
                    "type": "sample_done",
                    "index": i,
                    "total": total,
                    "sample": display_name,
                    "nuclei": n_nuclei,
                    "touching_merge_groups": merge_groups,
                    "touching_labels_removed": merge_removed,
                })

                # Explicitly release image arrays before the next sample.
                del raw
                del labels

            except Exception:
                errors.append({
                    "sample": display_name,
                    "error": traceback.format_exc()
                })

                self.batch_stardist_queue.put({
                    "type": "sample_error",
                    "index": i,
                    "total": total,
                    "sample": display_name,
                })

        self.batch_stardist_queue.put({
            "type": "done",
            "success": success,
            "errors": errors,
            "total": total,
        })

    def poll_batch_stardist_queue(self):
        finished = False

        while True:
            try:
                msg = (
                    self.batch_stardist_queue
                    .get_nowait()
                )
            except queue.Empty:
                break

            typ = msg.get("type")

            if typ == "progress":
                self.status.config(
                    text=(
                        f"Batch StarDist "
                        f"{msg['index']}/{msg['total']} "
                        f"{msg['sample']} "
                        f"[{msg['phase']}]"
                    )
                )

            elif typ == "sample_done":
                self.status.config(
                    text=(
                        f"Batch StarDist "
                        f"{msg['index']}/{msg['total']} "
                        f"{msg['sample']}: "
                        f"{msg['nuclei']} nuclei / "
                        f"auto-merge "
                        f"{msg.get('touching_labels_removed', 0)}"
                    )
                )

            elif typ == "sample_error":
                self.status.config(
                    text=(
                        f"Batch StarDist "
                        f"{msg['index']}/{msg['total']} "
                        f"{msg['sample']}: ERROR"
                    )
                )

            elif typ == "fatal":
                finished = True
                self.batch_stardist_running = False
                self.batch_stardist_button.config(
                    state=tk.NORMAL
                )
                messagebox.showerror(
                    "Batch StarDist",
                    msg["error"],
                    parent=self
                )

            elif typ == "done":
                finished = True
                self.batch_stardist_running = False
                self.batch_stardist_button.config(
                    state=tk.NORMAL
                )

                success = msg["success"]
                errors = msg["errors"]

                # Refresh statuses/counts from saved draft files.
                self.scan_samples()

                if errors and self.root_folder is not None:
                    err_path = (
                        self.root_folder
                        / "batch_stardist_errors.csv"
                    )
                    try:
                        _write_csv(
                            err_path,
                            errors,
                            fieldnames=[
                                "sample",
                                "error"
                            ]
                        )
                    except Exception:
                        pass

                if errors:
                    messagebox.showwarning(
                        "Batch StarDist finished",
                        (
                            f"成功: {len(success)}\n"
                            f"エラー: {len(errors)}\n\n"
                            "エラー内容は親フォルダの "
                            "batch_stardist_errors.csv "
                            "を確認してください。"
                        ),
                        parent=self
                    )
                else:
                    messagebox.showinfo(
                        "Batch StarDist finished",
                        (
                            f"{len(success)}枚のStarDistが完了しました。\n\n"
                            "結果は「StarDist済み・未確認」として保存されています。\n"
                            "各画像を開くと再推論せず即ロードされます。"
                        ),
                        parent=self
                    )

        if (
            self.batch_stardist_running
            and not finished
        ):
            self.after(
                100,
                self.poll_batch_stardist_queue
            )

    # ---------------- open ----------------

    def selected_index(self):
        sel = self.tree.selection()

        if not sel:
            return None

        try:
            return int(sel[0])
        except Exception:
            return None

    def open_selected(self):
        idx = self.selected_index()

        if idx is None:
            messagebox.showinfo(
                "Sample",
                "サンプルを選択してください。",
                parent=self
            )
            return

        self.open_index(idx)

    def open_first_unfinished(self):
        for i, s in enumerate(self.samples):
            if (
                s["dapi"] is not None
                and s["zo1"] is not None
                and not s["saved"]
            ):
                self.tree.selection_set(str(i))
                self.tree.see(str(i))
                self.open_index(i)
                return

        messagebox.showinfo(
            "Finished",
            "DAPI/ZO-1が揃った全サンプルが確認済み・保存済みです。",
            parent=self
        )

    def open_index(self, idx):
        if (
            idx < 0
            or idx >= len(self.samples)
        ):
            return

        s = self.samples[idx]

        if (
            s["dapi"] is None
            or s["zo1"] is None
        ):
            messagebox.showerror(
                "Missing files",
                (
                    f"{sample_display_name(s)}\n\n"
                    f"DAPI: "
                    f"{s['dapi'].name if s['dapi'] else 'NOT FOUND'}\n"
                    f"ZO-1: "
                    f"{s['zo1'].name if s['zo1'] else 'NOT FOUND'}"
                ),
                parent=self
            )
            return

        if (
            self.editor is not None
            and self.editor.winfo_exists()
        ):
            self.editor.lift()
            return

        try:
            self.editor = NucleiAnnotationWindow(
                self,
                s,
                self.settings.get(
                    "__stardist__",
                    {}
                ),
                self.on_saved,
                self.on_next,
            )

        except Exception:
            tb = traceback.format_exc()

            p = self.error_log_path()

            if p is not None:
                try:
                    with open(
                        p,
                        "a",
                        encoding="utf-8"
                    ) as f:
                        f.write(
                            "\n"
                            + "=" * 80
                            + "\n"
                            + f"Sample: {sample_display_name(s)}\n"
                            + f"Root folder: {s['folder']}\n"
                            + f"DAPI: {s['dapi']}\n"
                            + f"ZO-1: {s['zo1']}\n"
                            + tb
                            + "\n"
                        )
                except Exception:
                    pass

            self.editor = None

            messagebox.showerror(
                "Open error",
                tb,
                parent=self
            )

    # ---------------- tight-junction analysis ----------------

    def open_zo1_qc(self):
        idx = self.selected_index()

        if idx is None:
            messagebox.showinfo(
                "ZO-1 QC",
                "代表画像として使うサンプルを選択してください。",
                parent=self
            )
            return

        sample = self.samples[idx]

        if sample["zo1"] is None:
            messagebox.showerror(
                "ZO-1 QC",
                "選択サンプルにZO-1画像がありません。",
                parent=self
            )
            return

        try:
            params = (
                self.get_junction_analysis_params()
            )
        except Exception as exc:
            messagebox.showerror(
                "ZO-1 QC",
                str(exc),
                parent=self
            )
            return

        try:
            win = ZO1ProcessingQCWindow(
                self,
                sample,
                params,
                apply_callback=
                    self.apply_zo1_qc_settings
            )
            win.transient(self)
        except Exception:
            messagebox.showerror(
                "ZO-1 QC open error",
                traceback.format_exc(),
                parent=self
            )

    def apply_zo1_qc_settings(
        self,
        params
    ):
        self.zo1_rolling_ball_radius_var.set(
            int(
                params[
                    "zo1_rolling_ball_radius_px"
                ]
            )
        )

        self.zo1_tophat_radius_var.set(
            int(
                params[
                    "zo1_tophat_radius_px"
                ]
            )
        )

        self.zo1_tophat_threshold_var.set(
            float(
                params[
                    "zo1_tophat_threshold"
                ]
            )
        )

        self.zo1_min_object_size_var.set(
            int(
                params[
                    "zo1_min_object_size_px"
                ]
            )
        )

        # Save preprocessing + junction analysis settings globally.
        self._persist_analysis_settings_from_vars()

        self.analysis_status.config(
            text=(
                "ZO-1 QC設定を本解析へ反映: "
                f"RB={self.zo1_rolling_ball_radius_var.get()} px, "
                f"TopHat={self.zo1_tophat_radius_var.get()} px, "
                f"Thr={self.zo1_tophat_threshold_var.get():g}, "
                f"MinObj={self.zo1_min_object_size_var.get()} px"
            )
        )

    def get_junction_analysis_params(self):
        try:
            p = {
                "junction_lanes": int(
                    self.junction_lanes_var.get()
                ),
                "pair_distance_factor": float(
                    self.pair_distance_factor_var.get()
                ),
                "zo1_rolling_ball_radius_px": int(
                    self.zo1_rolling_ball_radius_var.get()
                ),
                "zo1_tophat_radius_px": int(
                    self.zo1_tophat_radius_var.get()
                ),
                "zo1_tophat_threshold": float(
                    self.zo1_tophat_threshold_var.get()
                ),
                "intact_coverage_threshold_pct": float(
                    self.intact_cutoff_var.get()
                ),
                "lane_peak_top_n_pixels": int(
                    self.lane_peak_top_n_var.get()
                ),
                "max_other_nucleus_fraction":
                    DEFAULT_MAX_OTHER_NUCLEUS_FRACTION,
                "zo1_min_object_size_px": int(
                    self.zo1_min_object_size_var.get()
                ),
                "lane_min_skeleton_pixels":
                    DEFAULT_LANE_MIN_SKELETON_PIXELS,
            }
        except Exception as exc:
            raise ValueError(
                f"解析設定を確認してください: {exc}"
            ) from exc

        if p["junction_lanes"] < 4:
            raise ValueError(
                "台形分割数は4以上にしてください。"
            )

        if p["pair_distance_factor"] < 1.0:
            raise ValueError(
                "Neighbor distance factorは1.0以上にしてください。"
            )

        if p["zo1_rolling_ball_radius_px"] < 1:
            raise ValueError(
                "ZO-1 rolling ball radiusは1 px以上にしてください。"
            )

        if p["zo1_tophat_radius_px"] < 1:
            raise ValueError(
                "White top-hat radiusは1 px以上にしてください。"
            )

        if p["zo1_tophat_threshold"] < 0:
            raise ValueError(
                "Top-hat thresholdは0以上にしてください。"
            )

        if p["zo1_min_object_size_px"] < 1:
            raise ValueError(
                "Min object sizeは1 px以上にしてください。"
            )

        if p["lane_peak_top_n_pixels"] < 1:
            raise ValueError(
                "Lane peak top N pixelsは1以上にしてください。"
            )

        if not (
            0
            <= p["intact_coverage_threshold_pct"]
            <= 100
        ):
            raise ValueError(
                "Intact cutoffは0～100%にしてください。"
            )

        return p

    def _set_analysis_buttons(self, enabled):
        state = (
            tk.NORMAL
            if enabled
            else tk.DISABLED
        )

        self.qc_button.config(
            state=state
        )

        self.analyze_selected_button.config(
            state=state
        )

        self.analyze_all_button.config(
            state=state
        )

    def _analysis_progress_from_worker(self, text):
        self.analysis_queue.put(
            {
                "type": "progress",
                "text": str(text)
            }
        )

    def analyze_selected_sample(self):
        if self.analysis_running:
            return

        idx = self.selected_index()

        if idx is None:
            messagebox.showinfo(
                "Analysis",
                "解析するサンプルを選択してください。",
                parent=self
            )
            return

        sample = self.samples[idx]

        if (
            sample["dapi"] is None
            or sample["zo1"] is None
        ):
            messagebox.showerror(
                "Analysis",
                "DAPIまたはZO-1がありません。",
                parent=self
            )
            return

        if not sample_output_paths(
            sample
        )["npz"].exists():
            messagebox.showwarning(
                "Analysis",
                "先にDAPI核をStarDistで確認・修正して保存してください。",
                parent=self
            )
            return

        try:
            params = self.get_junction_analysis_params()
            self._persist_analysis_settings_from_vars()
        except Exception as exc:
            messagebox.showerror(
                "Analysis settings",
                str(exc),
                parent=self
            )
            return

        self._start_analysis_worker(
            [sample],
            params,
            analyze_all=False
        )

    def analyze_all_saved_samples(self):
        if self.analysis_running:
            return

        targets = [
            s
            for s in self.samples
            if (
                s["dapi"] is not None
                and s["zo1"] is not None
                and sample_output_paths(
                    s
                )["npz"].exists()
            )
        ]

        if not targets:
            messagebox.showwarning(
                "Analysis",
                "保存済みDAPI核アノテーションのあるサンプルがありません。",
                parent=self
            )
            return

        try:
            params = self.get_junction_analysis_params()
            self._persist_analysis_settings_from_vars()
        except Exception as exc:
            messagebox.showerror(
                "Analysis settings",
                str(exc),
                parent=self
            )
            return

        self._start_analysis_worker(
            targets,
            params,
            analyze_all=True
        )

    def _start_analysis_worker(
        self,
        targets,
        params,
        analyze_all
    ):
        self.analysis_running = True
        self._set_analysis_buttons(False)

        self.analysis_status.config(
            text=(
                f"解析開始: {len(targets)} sample"
            )
        )

        thread = threading.Thread(
            target=self._analysis_worker,
            args=(
                targets,
                params,
                analyze_all
            ),
            daemon=True
        )

        thread.start()

        self.after(
            100,
            self.poll_analysis_queue
        )

    def _analysis_worker(
        self,
        targets,
        params,
        analyze_all
    ):
        summaries = []
        errors = []

        for i, sample in enumerate(
            targets,
            start=1
        ):
            try:
                self._analysis_progress_from_worker(
                    f"[{i}/{len(targets)}] "
                    f"{sample_display_name(sample)}"
                )

                summary = (
                    analyze_tight_junction_sample(
                        sample,
                        params,
                        progress_callback=
                            self._analysis_progress_from_worker
                    )
                )

                summaries.append(summary)

            except Exception:
                errors.append({
                    "sample_folder":
                        sample_display_name(sample),
                    "error":
                        traceback.format_exc()
                })

        if (
            analyze_all
            and self.root_folder is not None
            and summaries
        ):
            try:
                _write_csv(
                    self.root_folder
                    / ALL_JUNCTION_SUMMARY_FILENAME,
                    summaries
                )

                _write_group_summary(
                    self.root_folder
                    / GROUP_JUNCTION_SUMMARY_FILENAME,
                    summaries
                )
            except Exception:
                errors.append({
                    "sample_folder":
                        "__ALL_SUMMARY__",
                    "error":
                        traceback.format_exc()
                })

        if (
            errors
            and self.root_folder is not None
        ):
            try:
                err_path = (
                    self.root_folder
                    / "ZO1_junction_analysis_errors.csv"
                )

                _write_csv(
                    err_path,
                    errors,
                    fieldnames=[
                        "sample_folder",
                        "error"
                    ]
                )
            except Exception:
                pass

        self.analysis_queue.put({
            "type": "done",
            "summaries": summaries,
            "errors": errors,
            "analyze_all": analyze_all,
        })

    def poll_analysis_queue(self):
        handled_done = False

        while True:
            try:
                msg = (
                    self.analysis_queue.get_nowait()
                )
            except queue.Empty:
                break

            if msg["type"] == "progress":
                self.analysis_status.config(
                    text=msg["text"]
                )

            elif msg["type"] == "done":
                handled_done = True

                self.analysis_running = False
                self._set_analysis_buttons(True)

                summaries = msg["summaries"]
                errors = msg["errors"]

                self.analysis_status.config(
                    text=(
                        f"解析終了: "
                        f"{len(summaries)} success / "
                        f"{len(errors)} error"
                    )
                )

                if errors:
                    messagebox.showwarning(
                        "Analysis finished",
                        (
                            f"成功: {len(summaries)}\n"
                            f"エラー: {len(errors)}\n\n"
                            "親フォルダの "
                            "ZO1_junction_analysis_errors.csv "
                            "を確認してください。"
                        ),
                        parent=self
                    )
                else:
                    messagebox.showinfo(
                        "Analysis finished",
                        (
                            f"{len(summaries)}サンプルの"
                            "tight-junction解析が終了しました。\n\n"
                            f"{ALL_JUNCTION_SUMMARY_FILENAME}\n"
                            f"{GROUP_JUNCTION_SUMMARY_FILENAME}\n"
                            "を親フォルダに出力しました。"
                        ),
                        parent=self
                    )

        if self.analysis_running and not handled_done:
            self.after(
                100,
                self.poll_analysis_queue
            )

    # ---------------- callbacks ----------------

    def on_saved(self, sample, params):
        self.persist_stardist_settings(
            params
        )

        self.scan_samples()

        target_name = sample_display_name(
            sample
        )

        for i, s in enumerate(self.samples):
            if sample_display_name(s) == target_name:
                self.tree.selection_set(str(i))
                self.tree.see(str(i))
                break

    def on_next(self, current_sample):
        self.editor = None

        current_idx = None

        current_name = sample_display_name(
            current_sample
        )

        for i, s in enumerate(self.samples):
            if sample_display_name(s) == current_name:
                current_idx = i
                break

        if current_idx is None:
            self.open_first_unfinished()
            return

        for j in range(
            current_idx + 1,
            len(self.samples)
        ):
            s = self.samples[j]

            if (
                s["dapi"] is not None
                and s["zo1"] is not None
            ):
                self.tree.selection_set(str(j))
                self.tree.see(str(j))

                self.after(
                    80,
                    lambda jj=j: self.open_index(jj)
                )
                return

        messagebox.showinfo(
            "Finished",
            "最後のサンプルまで到達しました。",
            parent=self
        )
# ============================================================
# Mica MAX projection tab
# ============================================================

def read_interleaved_mica_tiff(path):
    """
    Read a Mica multipage grayscale TIFF whose page order is:
        Z1-Ch1, Z1-Ch2, ... Z2-Ch1, Z2-Ch2, ...
    """
    path = Path(path)

    with tifffile.TiffFile(str(path)) as tif:
        if len(tif.pages) < 1:
            raise ValueError("TIFFにページがありません。")

        first = tif.pages[0].asarray()

        if first.ndim != 2:
            raise ValueError(
                f"各ページが2D画像ではありません (shape={first.shape})。\n"
                "この処理はMicaのmultipage grayscale TIFFを想定しています。"
            )

        h, w = first.shape
        dtype = first.dtype
        n_pages = len(tif.pages)

        stack = np.empty(
            (n_pages, h, w),
            dtype=dtype
        )
        stack[0] = first

        for i in range(1, n_pages):
            arr = tif.pages[i].asarray()

            if arr.shape != (h, w):
                raise ValueError(
                    f"ページサイズが不一致です: "
                    f"page 1={(h, w)}, page {i + 1}={arr.shape}"
                )

            stack[i] = arr

    return stack


def make_mica_max_projections(
    path,
    n_channels,
    postfixes,
    output_dir=None
):
    """
    Split an interleaved Mica Z-stack into channels and create one
    maximum-intensity projection per channel.
    """
    path = Path(path).resolve()
    stack = read_interleaved_mica_tiff(path)
    n_pages = int(stack.shape[0])

    if n_pages % int(n_channels) != 0:
        raise ValueError(
            f"総ページ数 {n_pages} がチャンネル数 {n_channels} "
            "で割り切れません。\n"
            "ページ順またはチャンネル数を確認してください。"
        )

    n_z = n_pages // int(n_channels)

    zcyx = stack.reshape(
        n_z,
        int(n_channels),
        stack.shape[1],
        stack.shape[2]
    )

    out_dir = (
        Path(output_dir).resolve()
        if output_dir
        else path.parent
    )
    out_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    results = []

    for c in range(int(n_channels)):
        ch_stack = zcyx[:, c, :, :]
        max_img = np.max(
            ch_stack,
            axis=0
        )

        postfix = str(postfixes[c]).strip()

        if not postfix:
            postfix = f"_Ch{c + 1}_MAX"

        out_path = (
            out_dir
            / f"{path.stem}{postfix}.tif"
        )

        tifffile.imwrite(
            str(out_path),
            max_img,
            photometric="minisblack",
            metadata=None
        )

        results.append({
            "channel": c + 1,
            "n_z": n_z,
            "shape": max_img.shape,
            "dtype": str(max_img.dtype),
            "output": str(out_path),
        })

    return results


class MicaMaxProjectionTab(ttk.Frame):
    def __init__(self, parent, host_app=None):
        super().__init__(parent)

        self.host_app = host_app
        self.files = []

        self.output_dir = tk.StringVar(
            value=""
        )
        self.n_channels = tk.IntVar(
            value=2
        )

        # Defaults intentionally match the analysis tab's flat-file scanner.
        self.postfix_vars = [
            tk.StringVar(value="_DAPI"),
            tk.StringVar(value="_ZO1"),
            tk.StringVar(value="_Ch3_MAX"),
            tk.StringVar(value="_Ch4_MAX"),
        ]

        self.last_handoff_folder = None

        self._build_ui()
        self._update_postfix_state()

    def _build_ui(self):
        pad = 10

        intro = ttk.LabelFrame(
            self,
            text="Step 1: Mica multipage TIFF → Channel MAX projection",
            padding=pad
        )
        intro.pack(
            fill="x",
            padx=pad,
            pady=(pad, 6)
        )

        ttk.Label(
            intro,
            text=(
                "想定ページ順: "
                "Z1-Ch1, Z1-Ch2, ... → Z2-Ch1, Z2-Ch2, ...\n"
                "2ch時の既定出力は <sample>_DAPI.tif / <sample>_ZO1.tif "
                "なので、そのまま解析タブで認識できます。"
            ),
            justify="left"
        ).pack(
            anchor="w"
        )

        top = ttk.Frame(
            self,
            padding=(pad, 4, pad, 4)
        )
        top.pack(fill="x")

        ttk.Button(
            top,
            text="TIFFファイルを選択",
            command=self.select_files
        ).pack(
            side="left",
            padx=(0, 8)
        )

        ttk.Button(
            top,
            text="フォルダ内TIFFを追加",
            command=self.select_folder
        ).pack(
            side="left",
            padx=(0, 8)
        )

        ttk.Button(
            top,
            text="リストをクリア",
            command=self.clear_files
        ).pack(
            side="left"
        )

        list_frame = ttk.LabelFrame(
            self,
            text="処理対象ファイル",
            padding=pad
        )
        list_frame.pack(
            fill="both",
            expand=True,
            padx=pad,
            pady=(0, pad)
        )

        self.listbox = tk.Listbox(
            list_frame,
            height=8
        )
        self.listbox.pack(
            side="left",
            fill="both",
            expand=True
        )

        sb = ttk.Scrollbar(
            list_frame,
            orient="vertical",
            command=self.listbox.yview
        )
        sb.pack(
            side="right",
            fill="y"
        )
        self.listbox.configure(
            yscrollcommand=sb.set
        )

        settings = ttk.LabelFrame(
            self,
            text="MAX Projection 設定",
            padding=pad
        )
        settings.pack(
            fill="x",
            padx=pad,
            pady=(0, pad)
        )

        row1 = ttk.Frame(settings)
        row1.pack(
            fill="x",
            pady=4
        )

        ttk.Label(
            row1,
            text="チャンネル数:"
        ).pack(side="left")

        ch_spin = ttk.Spinbox(
            row1,
            from_=2,
            to=4,
            width=5,
            textvariable=self.n_channels,
            command=self._update_postfix_state
        )
        ch_spin.pack(
            side="left",
            padx=(8, 20)
        )
        ch_spin.bind(
            "<KeyRelease>",
            lambda e: self._update_postfix_state()
        )

        ttk.Label(
            row1,
            text="ページ順: Z1-Ch1, Z1-Ch2, ... → Z2-Ch1, Z2-Ch2, ..."
        ).pack(side="left")

        pf_frame = ttk.Frame(settings)
        pf_frame.pack(
            fill="x",
            pady=(8, 4)
        )

        self.postfix_entries = []

        for i in range(4):
            f = ttk.Frame(pf_frame)
            f.pack(
                side="left",
                padx=(0, 12)
            )

            ttk.Label(
                f,
                text=f"Ch{i + 1} postfix"
            ).pack(anchor="w")

            ent = ttk.Entry(
                f,
                textvariable=self.postfix_vars[i],
                width=16
            )
            ent.pack()

            self.postfix_entries.append(ent)

        out_frame = ttk.Frame(settings)
        out_frame.pack(
            fill="x",
            pady=(8, 4)
        )

        ttk.Label(
            out_frame,
            text="出力先:"
        ).pack(side="left")

        ttk.Entry(
            out_frame,
            textvariable=self.output_dir
        ).pack(
            side="left",
            fill="x",
            expand=True,
            padx=8
        )

        ttk.Button(
            out_frame,
            text="選択",
            command=self.select_output_dir
        ).pack(side="left")

        ttk.Button(
            out_frame,
            text="各元ファイルと同じ場所",
            command=lambda: self.output_dir.set("")
        ).pack(
            side="left",
            padx=(8, 0)
        )

        ttk.Label(
            settings,
            text=(
                "出力先が空欄の場合は各入力TIFFと同じフォルダに保存します。"
            )
        ).pack(
            anchor="w",
            pady=(4, 0)
        )

        run_frame = ttk.Frame(
            self,
            padding=(pad, 0, pad, pad)
        )
        run_frame.pack(
            fill="x"
        )

        self.run_btn = ttk.Button(
            run_frame,
            text="MAX画像を作成",
            command=self.start_processing
        )
        self.run_btn.pack(
            side="left"
        )

        self.handoff_btn = ttk.Button(
            run_frame,
            text="出力フォルダを解析タブで開く",
            command=self.handoff_to_analysis,
            state="disabled"
        )
        self.handoff_btn.pack(
            side="left",
            padx=(10, 0)
        )

        self.progress = ttk.Progressbar(
            run_frame,
            mode="determinate"
        )
        self.progress.pack(
            side="left",
            fill="x",
            expand=True,
            padx=(12, 0)
        )

        log_frame = ttk.LabelFrame(
            self,
            text="ログ",
            padding=pad
        )
        log_frame.pack(
            fill="both",
            expand=True,
            padx=pad,
            pady=(0, pad)
        )

        self.log = tk.Text(
            log_frame,
            height=9,
            wrap="word"
        )
        self.log.pack(
            fill="both",
            expand=True
        )

    def _update_postfix_state(self):
        try:
            n = int(
                self.n_channels.get()
            )
        except Exception:
            n = 2

        n = max(
            2,
            min(4, n)
        )

        for i, ent in enumerate(
            self.postfix_entries
        ):
            ent.configure(
                state=(
                    "normal"
                    if i < n
                    else "disabled"
                )
            )

    def select_files(self):
        paths = filedialog.askopenfilenames(
            title="Mica TIFFを選択",
            filetypes=[
                ("TIFF files", "*.tif *.tiff"),
                ("All files", "*.*"),
            ],
            parent=self
        )
        self._add_files(paths)

    def select_folder(self):
        folder = filedialog.askdirectory(
            title="TIFFフォルダを選択",
            parent=self
        )

        if not folder:
            return

        paths = sorted(
            [
                str(p)
                for p in Path(folder).iterdir()
                if (
                    p.is_file()
                    and p.suffix.lower()
                    in {".tif", ".tiff"}
                )
            ],
            key=natural_sort_key
        )

        self._add_files(paths)

    def _add_files(self, paths):
        existing = set(
            self.files
        )

        for p in paths:
            p = str(
                Path(p).resolve()
            )

            if p not in existing:
                self.files.append(p)
                self.listbox.insert(
                    "end",
                    p
                )
                existing.add(p)

    def clear_files(self):
        self.files.clear()
        self.listbox.delete(
            0,
            "end"
        )
        self.log.delete(
            "1.0",
            "end"
        )
        self.progress["value"] = 0
        self.last_handoff_folder = None
        self.handoff_btn.configure(
            state="disabled"
        )

    def select_output_dir(self):
        folder = filedialog.askdirectory(
            title="出力フォルダを選択",
            parent=self
        )

        if folder:
            self.output_dir.set(
                folder
            )

    def write_log(self, msg):
        self.after(
            0,
            self._write_log_ui,
            msg
        )

    def _write_log_ui(self, msg):
        self.log.insert(
            "end",
            str(msg) + "\n"
        )
        self.log.see(
            "end"
        )

    def start_processing(self):
        if not self.files:
            messagebox.showwarning(
                "入力なし",
                "TIFFファイルを選択してください。",
                parent=self
            )
            return

        try:
            n_channels = int(
                self.n_channels.get()
            )
        except Exception:
            messagebox.showerror(
                "設定エラー",
                "チャンネル数は2〜4を指定してください。",
                parent=self
            )
            return

        if n_channels not in (2, 3, 4):
            messagebox.showerror(
                "設定エラー",
                "チャンネル数は2〜4を指定してください。",
                parent=self
            )
            return

        postfixes = [
            self.postfix_vars[i].get()
            for i in range(n_channels)
        ]

        output_dir = (
            self.output_dir.get().strip()
            or None
        )

        self.run_btn.configure(
            state="disabled"
        )
        self.handoff_btn.configure(
            state="disabled"
        )
        self.progress.configure(
            maximum=len(self.files),
            value=0
        )
        self.log.delete(
            "1.0",
            "end"
        )
        self.last_handoff_folder = None

        threading.Thread(
            target=self._process_worker,
            args=(
                n_channels,
                postfixes,
                output_dir
            ),
            daemon=True
        ).start()

    def _process_worker(
        self,
        n_channels,
        postfixes,
        output_dir
    ):
        ok = 0
        ng = 0
        produced_dirs = []

        for i, path in enumerate(
            self.files,
            start=1
        ):
            try:
                self.write_log(
                    f"[{i}/{len(self.files)}] "
                    f"{Path(path).name}"
                )

                results = make_mica_max_projections(
                    path,
                    n_channels,
                    postfixes,
                    output_dir
                )

                for r in results:
                    self.write_log(
                        f"  Ch{r['channel']}: "
                        f"Z={r['n_z']} → "
                        f"{Path(r['output']).name}"
                    )

                if results:
                    produced_dirs.append(
                        str(
                            Path(
                                results[0]["output"]
                            ).parent.resolve()
                        )
                    )

                ok += 1

            except Exception as e:
                ng += 1

                self.write_log(
                    f"  ERROR: {e}"
                )
                self.write_log(
                    traceback.format_exc()
                )

            self.after(
                0,
                lambda v=i: self.progress.configure(
                    value=v
                )
            )

        self.write_log("")
        self.write_log(
            f"完了: 成功 {ok} ファイル / "
            f"エラー {ng} ファイル"
        )

        unique_dirs = sorted(
            set(produced_dirs)
        )

        if len(unique_dirs) == 1:
            self.last_handoff_folder = (
                unique_dirs[0]
            )
            self.after(
                0,
                lambda: self.handoff_btn.configure(
                    state="normal"
                )
            )
        else:
            self.last_handoff_folder = None

        self.after(
            0,
            lambda: self.run_btn.configure(
                state="normal"
            )
        )

        if ng == 0:
            self.after(
                0,
                lambda: messagebox.showinfo(
                    "完了",
                    f"{ok}ファイルの処理が完了しました。",
                    parent=self
                )
            )
        else:
            self.after(
                0,
                lambda: messagebox.showwarning(
                    "完了",
                    f"成功 {ok}、エラー {ng}。\n"
                    "詳細はログを確認してください。",
                    parent=self
                )
            )

    def handoff_to_analysis(self):
        folder = self.last_handoff_folder

        if not folder:
            messagebox.showinfo(
                "解析フォルダ",
                "解析へ渡せる単一の出力フォルダがありません。",
                parent=self
            )
            return

        if self.host_app is not None:
            self.host_app.open_analysis_folder(
                folder
            )


# ============================================================
# Integrated application
# ============================================================

class IntegratedMicaZO1App(tk.Tk):
    def __init__(self):
        super().__init__()

        self.geometry(
            "1280x820"
        )
        self.minsize(
            1000,
            650
        )

        # Restore the last language selection. Japanese is the default.
        initial_lang = "ja"
        try:
            if _LANGUAGE_PREF_PATH.exists():
                _pref = load_json(_LANGUAGE_PREF_PATH)
                initial_lang = str(_pref.get("language", "ja"))
        except Exception:
            initial_lang = "ja"
        set_app_language(initial_lang)

        self.language_bar = ttk.Frame(self, padding=(8, 5))
        self.language_bar.pack(fill="x")

        ttk.Label(
            self.language_bar,
            text="Language / 言語:"
        ).pack(side="left", padx=(0, 8))

        self.ja_btn = ttk.Button(
            self.language_bar,
            text="日本語",
            command=lambda: self.set_language("ja")
        )
        self.ja_btn.pack(side="left", padx=(0, 4))

        self.en_btn = ttk.Button(
            self.language_bar,
            text="English",
            command=lambda: self.set_language("en")
        )
        self.en_btn.pack(side="left")

        self.notebook = ttk.Notebook(
            self
        )
        self.notebook.pack(
            fill="both",
            expand=True
        )

        self.mica_tab = MicaMaxProjectionTab(
            self.notebook,
            host_app=self
        )

        self.analysis_tab = ZO1AnalysisTab(
            self.notebook,
            host_app=self
        )

        self.notebook.add(
            self.mica_tab,
            text="1. Mica MAX Projection"
        )

        self.notebook.add(
            self.analysis_tab,
            text="2. DAPI / ZO-1 Analysis"
        )

        self.protocol(
            "WM_DELETE_WINDOW",
            self.on_close
        )

        self.set_language(initial_lang, persist=False)

        # The analysis code creates editor/QC Toplevels dynamically. A lightweight
        # refresh keeps newly created windows synchronized with the selected language.
        self.after(500, self._language_refresh_tick)

    def set_language(self, lang, persist=True):
        set_app_language(lang)
        lang = get_app_language()

        if persist:
            try:
                save_json(
                    _LANGUAGE_PREF_PATH,
                    {"language": lang}
                )
            except Exception:
                pass

        # Keep the language controls themselves unchanged and use state as feedback.
        try:
            self.ja_btn.configure(
                state=("disabled" if lang == "ja" else "normal")
            )
            self.en_btn.configure(
                state=("disabled" if lang == "en" else "normal")
            )
        except Exception:
            pass

        # Explicit titles/tab labels provide stable wording in both languages.
        if lang == "ja":
            self.title("Mica → DAPI / ZO-1 Junction解析 GUI")
            self.notebook.tab(self.mica_tab, text="1. Mica MAX Projection")
            self.notebook.tab(self.analysis_tab, text="2. DAPI / ZO-1 解析")
        else:
            self.title("Mica → DAPI / ZO-1 Junction Analysis GUI")
            self.notebook.tab(self.mica_tab, text="1. Mica MAX Projection")
            self.notebook.tab(self.analysis_tab, text="2. DAPI / ZO-1 Analysis")

        translate_widget_tree(self, lang)

        # Restore fixed language-button labels after recursive translation.
        try:
            self.ja_btn.configure(text="日本語")
            self.en_btn.configure(text="English")
        except Exception:
            pass

    def _language_refresh_tick(self):
        try:
            translate_widget_tree(self, get_app_language())
            self.ja_btn.configure(text="日本語")
            self.en_btn.configure(text="English")
        except Exception:
            pass
        try:
            self.after(1200, self._language_refresh_tick)
        except Exception:
            pass

    def open_analysis_folder(self, folder):
        # Keep the analysis pairing rule synchronized with the MAX-projection
        # postfixes used to create Ch1/Ch2 outputs.
        try:
            if int(self.mica_tab.n_channels.get()) >= 2:
                self.analysis_tab.dapi_postfix_var.set(
                    str(self.mica_tab.postfix_vars[0].get()).strip()
                )
                self.analysis_tab.zo1_postfix_var.set(
                    str(self.mica_tab.postfix_vars[1].get()).strip()
                )
                self.analysis_tab._update_postfix_hint()
                self.analysis_tab._persist_analysis_settings_from_vars()
        except Exception:
            pass

        self.notebook.select(
            self.analysis_tab
        )
        self.analysis_tab.set_parent_folder(
            folder
        )
        translate_widget_tree(self, get_app_language())

    def on_close(self):
        try:
            self.analysis_tab.on_main_close()
        except Exception:
            pass

        try:
            editor = getattr(
                self.analysis_tab,
                "editor",
                None
            )

            if editor is not None:
                editor.destroy()
        except Exception:
            pass

        self.destroy()

def main():
    app = IntegratedMicaZO1App()
    app.mainloop()


if __name__ == "__main__":
    main()
