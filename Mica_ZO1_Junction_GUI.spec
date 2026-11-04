# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_all

datas=[]; binaries=[]; hiddenimports=[]
for pkg in ["stardist", "csbdeep", "skimage", "tifffile", "tensorflow"]:
    try:
        d,b,h=collect_all(pkg); datas+=d; binaries+=b; hiddenimports+=h
    except Exception:
        pass

a=Analysis(["mica_zo1_junction_gui.py"], pathex=[], binaries=binaries, datas=datas,
           hiddenimports=hiddenimports, hookspath=[], hooksconfig={}, runtime_hooks=[],
           excludes=[], noarchive=False)
pyz=PYZ(a.pure)
exe=EXE(pyz, a.scripts, a.binaries, a.datas, [], name="Mica_ZO1_Junction_GUI",
        debug=False, bootloader_ignore_signals=False, strip=False, upx=True,
        upx_exclude=[], runtime_tmpdir=None, console=False,
        disable_windowed_traceback=False)
