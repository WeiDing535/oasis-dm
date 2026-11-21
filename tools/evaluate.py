#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
evaluate.py -- 通用、健壮的「相位改进」评估模块
================================================

比较多种密度修饰程序（自写程序 / CCP4 dm / cparrot / 原始实验相位）的相位
改进效果，按 HKL 做 dict join 对齐（不假设不同 MTZ 行顺序一致），计算一整套
指标并输出 JSON / Markdown 表 / PNG 图。

指标依据：
  - research/env_data.md §7「评估方案」
  - research/theory.md    §1.6「图相关系数 CC」

运行环境：/opt/miniconda3/bin/python（gemmi 0.7.1 + numpy 2.3 + scipy + matplotlib）

固定 CLI（主程序调用方式）：
  /opt/miniconda3/bin/python tests/evaluate.py \
      --case <case.mtz> \
      --truth-col PHIC [--truth-fcol FC] \
      --fo-col FOBS|FPHASED --fom-col FOM --free-col FreeR|FreeR_flag \
      --methods "name1=path1.mtz:PHI_COL[:FOM_COL],name2=path2.mtz:PHI_COL[:FOM_COL],..." \
      --out <prefix>

自测（产物只写 /tmp）：
  /opt/miniconda3/bin/python tests/evaluate.py --selftest
"""

import argparse
import json
import math
import os
import sys
import tempfile

import numpy as np
import gemmi

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------------
# 全局：警告收集 / 轻量日志
# ----------------------------------------------------------------------------
WARNINGS = []


def warn(msg):
    text = "WARNING: " + msg
    print(text)
    WARNINGS.append(msg)


def fail(msg):
    print("ERROR: " + msg)
    sys.exit(2)


def sanitize(o):
    """递归把 numpy 标量转 Python，并把 NaN/inf 换成 None，保证 JSON 合法。"""
    if isinstance(o, dict):
        return {k: sanitize(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [sanitize(v) for v in o]
    if isinstance(o, np.generic):
        o = o.item()
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    return o


def fmt(v, nd=4):
    """表格数值格式化，None/非有限 -> N/A。"""
    if v is None:
        return "N/A"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not math.isfinite(f):
        return "N/A"
    return f"{f:.{nd}f}"


# ----------------------------------------------------------------------------
# 数学工具
# ----------------------------------------------------------------------------
def wrap(d):
    """把度 wrap 到 (-180, 180]（与 env_data.md §7.1 一致，-180 映射到 +180）。"""
    r = (np.asarray(d, dtype=float) + 180.0) % 360.0 - 180.0
    r = np.where(r <= -180.0, r + 360.0, r)
    return r


def normalize_phase(phi):
    """
    相位单位自动判定：max(|phi|) <= 3.15 -> 弧度(转度)，否则按度。
    最终 wrap 到 (-180,180]。返回 (deg_wrapped, unit_str)；NaN 原样保留。
    """
    phi = np.asarray(phi, dtype=float)
    finite = np.isfinite(phi)
    mx = float(np.max(np.abs(phi[finite]))) if finite.any() else 0.0
    if mx <= 3.15:
        deg = np.rad2deg(phi)
        unit = "rad"
    else:
        deg = phi.copy()
        unit = "deg"
    return wrap(deg), unit


def pearson(a, b):
    """两数组（去均值归一化）的 Pearson CC；退化时返回 None。"""
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    if a.size < 2 or b.size < 2 or a.size != b.size:
        return None
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        return None
    a = a - a.mean()
    b = b - b.mean()
    den = math.sqrt(float((a * a).sum()) * float((b * b).sum()))
    if den <= 0:
        return None
    return float((a * b).sum() / den)


def wstats(idx, w, absdphi, cosdphi, Fo, dphi_rad):
    """
    在一组反射下标上计算加权/未加权相位误差指标。
    w = m_init * Fo（无 FOM 时为 Fo）；若 Σw<=0 退化为等权。
    返回 dict 或 None（空集）。
    """
    if idx is None:
        return None
    idx = np.asarray(idx, dtype=int)
    if idx.size == 0:
        return None
    w_ = w[idx]
    ad_ = absdphi[idx]
    cd_ = cosdphi[idx]
    Fo_ = Fo[idx]
    dr_ = dphi_rad[idx]
    sw = float(w_.sum())
    if sw <= 0:
        sw = float(idx.size)
        w_ = np.ones_like(w_)
    wMPE = float((w_ * ad_).sum() / sw)
    MPE = float(ad_.mean())
    mean_cos = float((w_ * cd_).sum() / sw)
    ww = w_ * Fo_ ** 2                       # σA_like_C 权重 = w · Fo²
    sww = float(ww.sum())
    C = float(np.abs((ww * np.exp(1j * dr_)).sum()) / sww) if sww > 0 else None
    return dict(wMPE=wMPE, MPE=MPE, mean_cos_dphi=mean_cos,
                sigmaA_like_C=C, n=int(idx.size))


def amplitude_bins(amp, mask, w, absdphi, dphi_rad, cosdphi, n=10):
    """相位误差按「振幅分位」的变化：等计数分成 n 个 bin。"""
    sel = mask & np.isfinite(amp)
    idx = np.where(sel)[0]
    if idx.size < n:
        return None
    a = amp[idx]
    order = idx[np.argsort(a)]
    out = []
    for chunk in np.array_split(order, n):
        st = wstats(chunk, w, absdphi, cosdphi, np.ones_like(w), dphi_rad)
        aa = amp[chunk]
        out.append(dict(amp_min=float(aa.min()), amp_max=float(aa.max()),
                        amp_center=float(aa.mean()),
                        n=int(chunk.size),
                        wMPE=st["wMPE"] if st else None,
                        MPE=st["MPE"] if st else None,
                        mean_cos_dphi=st["mean_cos_dphi"] if st else None))
    return out


def assign_shells(invd2, edges):
    """按全局 1/d² 分位边界把反射分到 10 个等体积壳（0..9）。"""
    return np.clip(np.searchsorted(edges, invd2, side="right"), 0, 9)


# ----------------------------------------------------------------------------
# MTZ 列 / 元数据读取
# ----------------------------------------------------------------------------
def resolve_column(mtz, spec):
    """
    spec 可为 'LAB' 或 'A|B|C'（候选列，取第一个存在的）。
    返回 (label, Column) 或 (None, None)。先精确匹配，再大小写不敏感匹配。
    """
    if spec is None:
        return None, None
    labels = [s.strip() for s in str(spec).split("|") if s.strip()]
    for lab in labels:
        c = mtz.column_with_label(lab)
        if c is not None:
            return c.label, c
    ci = {c.label.upper(): c for c in mtz.columns}
    for lab in labels:
        c = ci.get(lab.upper())
        if c is not None:
            return c.label, c
    return None, None


def get_col(mtz, spec):
    """返回 (label, float64 array) 或 (None, None)。"""
    lab, c = resolve_column(mtz, spec)
    if c is None:
        return None, None
    return lab, np.asarray(c.array, dtype=np.float64)


def hkl_index(mtz):
    """(h,k,l)tuple -> 行号 的 dict。"""
    ma = mtz.make_miller_array()
    d = {}
    for i in range(ma.shape[0]):
        d[(int(ma[i, 0]), int(ma[i, 1]), int(ma[i, 2]))] = i
    return d


def load_case(args):
    """读取 case 文件：真值相位、观测振幅、初始 FOM、free 标记、可选真值振幅。"""
    path = args.case
    if not path or not os.path.isfile(path):
        fail(f"case 文件不存在: {path}")
    try:
        mtz = gemmi.read_mtz_file(path)
    except Exception as e:
        fail(f"无法读取 case MTZ: {path} ({e})")
    try:
        mtz.update_reso()
    except Exception:
        pass

    truth_lab, phi_true_raw = get_col(mtz, args.truth_col)
    if truth_lab is None:
        fail(f"case 中找不到真值相位列 --truth-col {args.truth_col!r}；"
             f"可用列: {mtz.column_labels()}")
    fo_lab, Fo = get_col(mtz, args.fo_col)
    if fo_lab is None:
        fail(f"case 中找不到观测振幅列 --fo-col {args.fo_col!r}；"
             f"可用列: {mtz.column_labels()}")

    phi_true_deg, phi_true_unit = normalize_phase(phi_true_raw)

    fom_lab, m_init = get_col(mtz, args.fom_col)
    if fom_lab is None:
        warn(f"case 缺少初始 FOM 列 (--fom-col {args.fom_col!r})；权重退化为 w = Fo")
    free_lab, free = get_col(mtz, args.free_col)
    if free_lab is None:
        warn(f"case 缺少 free 列 (--free-col {args.free_col!r})；work/free 指标填 null")
    truth_f_lab, truth_f = get_col(mtz, args.truth_fcol) if args.truth_fcol else (None, None)

    ma = mtz.make_miller_array()
    H = ma[:, 0].astype(int)
    K = ma[:, 1].astype(int)
    L = ma[:, 2].astype(int)
    d = np.asarray(mtz.make_d_array(), dtype=np.float64)
    invd2 = np.asarray(mtz.make_1_d2_array(), dtype=np.float64)

    # 10 个等体积壳：用 case 全部有效反射的 1/d² 分位边界（所有方法共用）
    good = np.isfinite(invd2) & np.isfinite(Fo) & (Fo > 0)
    if good.sum() >= 10:
        qs = np.quantile(invd2[good], np.linspace(0.0, 1.0, 11))
        edges = qs[1:-1].copy()                     # 9 个内部分位边界
    else:
        edges = np.array([np.nanmin(invd2[good]) + (k + 1) / 10.0 *
                          (np.nanmax(invd2[good]) - np.nanmin(invd2[good]))
                          for k in range(9)]) if good.sum() >= 2 else np.linspace(0, 1, 9)
        warn("case 有效反射不足 10 个，分辨率壳边界退化为线性划分")
    sid_all = assign_shells(invd2, edges)
    shells_meta = []
    for s in range(10):
        sel = good & (sid_all == s)
        if sel.any():
            shells_meta.append(dict(index=s,
                                    d_center=float(np.mean(d[sel])),
                                    d_min=float(np.min(d[sel])),
                                    d_max=float(np.max(d[sel])),
                                    n_ref=int(sel.sum())))
        else:
            shells_meta.append(dict(index=s, d_center=None,
                                    d_min=None, d_max=None, n_ref=0))

    cell = mtz.cell
    dmin = float(np.nanmin(d[np.isfinite(d)])) if np.isfinite(d).any() else None
    dmax = float(np.nanmax(d[np.isfinite(d)])) if np.isfinite(d).any() else None
    fmark = getattr(args, 'free_marker', 0)
    nfree = int(np.count_nonzero(free == fmark)) if free is not None else 0
    nwork = int(np.count_nonzero(np.isfinite(free) & (free != fmark))) if free is not None else 0

    return dict(
        mtz=mtz, file=os.path.abspath(path),
        H=H, K=K, L=L, hkl_index=hkl_index(mtz),
        Fo=Fo, phi_true_deg=phi_true_deg, phi_true_unit=phi_true_unit,
        m_init=m_init, free=free, truth_f=truth_f,
        d=d, invd2=invd2, shell_edges=edges, shells_meta=shells_meta,
        cell=cell, sg=mtz.spacegroup, volume=float(cell.volume),
        labels=dict(truth_col=truth_lab, fo_col=fo_lab, fom_col=fom_lab,
                    free_col=free_lab, truth_fcol=truth_f_lab),
        meta=dict(file=os.path.abspath(path),
                  space_group=str(mtz.spacegroup.hm),
                  cell=[float(cell.a), float(cell.b), float(cell.c),
                        float(cell.alpha), float(cell.beta), float(cell.gamma)],
                  d_min=dmin, d_max=dmax, n_reflections=int(mtz.nreflections),
                  n_free=nfree, n_work=nwork,
                  truth_phase_unit=phi_true_unit,
                  columns=dict(truth_col=truth_lab, fo_col=fo_lab,
                               fom_col=fom_lab, free_col=free_lab,
                               truth_fcol=truth_f_lab)),
    )


def parse_methods(spec_str):
    """解析 --methods：name=path:PHI[:FOM],..."""
    out = []
    if not spec_str:
        fail("--methods 为空")
    for item in spec_str.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            warn(f"--methods 项缺少 '='，跳过: {item!r}")
            continue
        name, rest = item.split("=", 1)
        name, rest = name.strip(), rest.strip()
        parts = rest.split(":")
        if len(parts) == 2:
            path, phi, fom = parts[0], parts[1], None
        elif len(parts) == 3:
            path, phi, fom = parts[0], parts[1], parts[2]
        elif len(parts) > 3:                       # 路径可能含 ':'
            path, phi, fom = ":".join(parts[:-2]), parts[-2], parts[-1]
        else:
            warn(f"--methods 项格式错误（应为 name=path:PHI[:FOM]），跳过: {item!r}")
            continue
        out.append(dict(name=name, path=path.strip(),
                        phi=phi.strip(), fom=(fom.strip() if fom else None)))
    if not out:
        fail("--methods 未解析出任何方法")
    return out


def load_method(spec):
    """读取单个方法 MTZ 的相位列（+ 可选 FOM 列），失败返回 None 并给出明确提示。"""
    name, path = spec["name"], spec["path"]
    if not os.path.isfile(path):
        warn(f"方法 [{name}] 的 MTZ 不存在: {path}；已跳过该方法")
        return None
    try:
        mtz = gemmi.read_mtz_file(path)
    except Exception as e:
        warn(f"方法 [{name}] 的 MTZ 无法读取 ({e})；已跳过该方法")
        return None
    phi_lab, phi_raw = get_col(mtz, spec["phi"])
    if phi_lab is None:
        warn(f"方法 [{name}] 中找不到相位列 {spec['phi']!r}；"
             f"可用列: {mtz.column_labels()}；已跳过该方法")
        return None
    phi_deg, phi_unit = normalize_phase(phi_raw)
    fom_lab, fom = (get_col(mtz, spec["fom"]) if spec["fom"] else (None, None))
    if spec["fom"] and fom_lab is None:
        warn(f"方法 [{name}] 中找不到 FOM 列 {spec['fom']!r}；mean FOM 记为 null")
    return dict(name=name, file=os.path.abspath(path),
                phi_col=phi_lab, fom_col=fom_lab, phase_unit=phi_unit,
                phi_deg=phi_deg, fom=fom, hkl_index=hkl_index(mtz))


def common_indices(case, method):
    """HKL dict join：返回 (case 行号数组, method 行号数组)，顺序一致。"""
    c = case["hkl_index"]
    m = method["hkl_index"]
    keys = c.keys() & m.keys()
    ci = np.fromiter((c[k] for k in keys), dtype=int, count=len(keys))
    mi = np.fromiter((m[k] for k in keys), dtype=int, count=len(keys))
    return ci, mi


# ----------------------------------------------------------------------------
# map CC：gemmi 对称性展开 + numpy FFT（约定已与 gemmi 官方图核对 CC=1.0）
# ----------------------------------------------------------------------------
def build_map_mtz(H, K, L, F, P, sg, cell):
    """临时构造只含 H,K,L,F,PHI 的 gemmi.Mtz（F 的 NaN/无效值须由调用方填 0）。"""
    m = gemmi.Mtz(with_base=True)
    m.spacegroup = sg
    m.cell = cell
    m.add_dataset("map")
    m.add_column("F", "F")
    m.add_column("PHI", "P")
    m.set_data(np.column_stack([H, K, L, F, P]).astype(np.float32))
    return m


def grid_shape(case):
    """网格 size：d_min/3 采样；优先取 transform_f_phi_to_map(sample_rate=3) 的 shape。"""
    H, K, L = case["H"], case["K"], case["L"]
    Fo, phi = case["Fo"], case["phi_true_deg"]
    valid = np.isfinite(Fo) & (Fo > 0) & np.isfinite(phi)
    F = np.where(valid, Fo, 0.0)
    P = np.where(valid, phi, 0.0)
    try:
        m = build_map_mtz(H, K, L, F, P, case["sg"], case["cell"])
        g = m.transform_f_phi_to_map("F", "PHI", sample_rate=3)
        return [int(g.nu), int(g.nv), int(g.nw)]
    except Exception as e:
        warn(f"transform_f_phi_to_map 定网格失败 ({e})，改用 d_min/3 + 2/3/5 因子")
        cell = case["cell"]
        dmin = case["meta"]["d_min"] or 3.0
        step = dmin / 3.0

        def nice(n):
            best = None
            for a2 in range(0, 14):
                for b3 in range(0, 9):
                    for c5 in range(0, 7):
                        v = (2 ** a2) * (3 ** b3) * (5 ** c5)
                        if v >= n and (best is None or v < best):
                            best = v
            return best or int(n)

        return [nice(math.ceil(cell.a / step)),
                nice(math.ceil(cell.b / step)),
                nice(math.ceil(cell.c / step))]


def density(H, K, L, F, P, shape, sg, cell, N, V):
    """np.fft.ifftn(grid.conj()).real * (N/V) 得实空间密度。"""
    m = build_map_mtz(H, K, L, F, P, sg, cell)
    gr = m.get_f_phi_on_grid("F", "PHI", size=shape)
    arr = np.asarray(gr)                          # ReciprocalComplexGrid -> complex64
    return np.fft.ifftn(arr.conj()).real * (N / V)


def map_cc(case, ctx, H, K, L, Fo, phi_true, phi_out, valid):
    """
    参考图 = Fo·exp(i φ_true)，方法图 = Fo·exp(i φ_out)。两者用同一 HKL 公共集、
    同一 F 数组、同一网格；F 无效处填 0（避免 0*exp(i·NaN)=NaN 污染 FFT）。
    返回 (全胞 CC, 高密度区 CC, 高密度区体积分数, 阈值)。
    """
    shape, N, sg, cell, V = ctx
    vm = valid & np.isfinite(Fo) & (Fo > 0) & np.isfinite(phi_true) & np.isfinite(phi_out)
    F = np.where(vm, Fo, 0.0)
    P_ref = np.where(vm, phi_true, 0.0)
    P_out = np.where(vm, phi_out, 0.0)
    rho_ref = density(H, K, L, F, P_ref, shape, sg, cell, N, V)
    rho_out = density(H, K, L, F, P_out, shape, sg, cell, N, V)
    cc = pearson(rho_ref, rho_out)
    thr = float(np.quantile(rho_ref, 0.75))
    hi = rho_ref > thr
    cc_hi = pearson(rho_ref[hi], rho_out[hi])
    return cc, cc_hi, float(hi.mean()), thr


# ----------------------------------------------------------------------------
# 单方法评估
# ----------------------------------------------------------------------------
def evaluate_method(case, method, ctx, args):
    ci, mi = common_indices(case, method)
    ncase = len(case["H"])
    if ci.size == 0:
        warn(f"方法 [{method['name']}] 与 case 无公共 HKL；已跳过该方法")
        return None
    if ci.size < 0.5 * ncase:
        warn(f"方法 [{method['name']}] 与 case 仅 {ci.size}/{ncase} 条 HKL 匹配"
             f"（可能重索引/ASU 不一致），指标仅基于公共集")

    Fo = case["Fo"][ci]
    phi_true = case["phi_true_deg"][ci]
    phi_out = method["phi_deg"][mi]
    H, K, L = case["H"][ci], case["K"][ci], case["L"][ci]
    invd2 = case["invd2"][ci]
    m_init = case["m_init"][ci] if case["m_init"] is not None else None
    free = case["free"][ci] if case["free"] is not None else None
    truth_f = case["truth_f"][ci] if case["truth_f"] is not None else None
    m_out = method["fom"][mi] if method["fom"] is not None else None

    valid = (np.isfinite(Fo) & (Fo > 0) &
             np.isfinite(phi_true) & np.isfinite(phi_out))
    if valid.sum() == 0:
        warn(f"方法 [{method['name']}] 公共集无有效反射（全为 NaN/无效）；已跳过")
        return None

    w_base = (np.ones_like(Fo) if m_init is None
              else np.where(np.isfinite(m_init) & (m_init >= 0), m_init, 1.0))
    w = w_base * Fo
    dphi = wrap(phi_out - phi_true)               # 有符号 Δφ，(-180,180]
    absdphi = np.abs(dphi)
    dphi_rad = np.deg2rad(dphi)
    cosdphi = np.cos(dphi_rad)

    g = wstats(np.where(valid)[0], w, absdphi, cosdphi, Fo, dphi_rad)

    # 方法自带 FOM 的均值
    mean_fom = None
    if m_out is not None:
        x = m_out[valid & np.isfinite(m_out)]
        mean_fom = float(x.mean()) if x.size else None

    # work / free（free 标记值由 --free-marker 指定，CCP4 默认 0，二值 FreeR 常用 1）
    fmark = getattr(args, 'free_marker', 0)
    wf = dict(wMPE_work=None, wMPE_free=None, n_work=0, n_free=0)
    if free is not None:
        fmask = valid & (free == fmark)
        wmask = valid & np.isfinite(free) & (free != fmark)
        sf = wstats(np.where(fmask)[0], w, absdphi, cosdphi, Fo, dphi_rad)
        sw = wstats(np.where(wmask)[0], w, absdphi, cosdphi, Fo, dphi_rad)
        wf["wMPE_free"] = sf["wMPE"] if sf else None
        wf["wMPE_work"] = sw["wMPE"] if sw else None
        wf["n_free"] = int(fmask.sum())
        wf["n_work"] = int(wmask.sum())
        if wf["n_free"] == 0:
            warn(f"方法 [{method['name']}] free 集为空（free 列无 {fmark} 值，"
                 f"可试 --free-marker 1）；wMPE_free 填 null")
    # 若 free 列缺失，load_case 已给出全局提示

    # 10 个等体积分辨率壳
    sid = assign_shells(invd2, case["shell_edges"])
    shells = []
    for s in range(10):
        smask = valid & (sid == s)
        st = wstats(np.where(smask)[0], w, absdphi, cosdphi, Fo, dphi_rad)
        mfs = None
        if m_out is not None:
            x = m_out[smask & np.isfinite(m_out)]
            mfs = float(x.mean()) if x.size else None
        shells.append(dict(index=s,
                           d_center=case["shells_meta"][s]["d_center"],
                           wMPE=st["wMPE"] if st else None,
                           MPE=st["MPE"] if st else None,
                           mean_cos_dphi=st["mean_cos_dphi"] if st else None,
                           mean_FOM=mfs,
                           n=int(smask.sum())))

    # Δφ 直方图（有符号，用于 _phase_error.png）
    counts, edges = np.histogram(dphi[valid], bins=36, range=(-180.0, 180.0))
    hist = dict(bin_edges=edges.tolist(), counts=counts.astype(int).tolist())

    # 可选：相位误差按振幅分位变化（--truth-fcol 提供时）
    amp = None
    if case["truth_f"] is not None:
        amp = dict(by_fobs=amplitude_bins(Fo, valid, w, absdphi, dphi_rad, cosdphi),
                   by_fc=amplitude_bins(truth_f, valid, w, absdphi, dphi_rad, cosdphi))

    # map CC
    cc, cc_hi, hi_frac, hi_thr = map_cc(case, ctx, H, K, L, Fo, phi_true, phi_out, valid)

    return dict(
        file=method["file"], phi_col=method["phi_col"],
        fom_col=method["fom_col"], phase_unit=method["phase_unit"],
        n_common=int(ci.size), n_valid=int(valid.sum()),
        global_=dict(wMPE=g["wMPE"], MPE=g["MPE"],
                     mean_cos_dphi=g["mean_cos_dphi"],
                     sigmaA_like_C=g["sigmaA_like_C"],
                     mean_FOM=mean_fom, n=g["n"]),
        work_free=wf,
        map_CC=dict(whole_cell=cc, high_density=cc_hi,
                    high_density_threshold=hi_thr,
                    high_density_fraction=hi_frac),
        shells=shells,
        phase_error_hist=hist,
        phase_error_vs_amplitude=amp,
    )


# ----------------------------------------------------------------------------
# 输出：JSON / Markdown 表 / 3 张 PNG
# ----------------------------------------------------------------------------
def method_colors(names, best):
    """best（wMPE 最小）用红色（红=好），其余用调色板。"""
    palette = ["#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e",
               "#17becf", "#8c564b", "#7f7f7f", "#bcbd22"]
    colors = {}
    k = 0
    for n in names:
        if n == best:
            colors[n] = "#d62728"
        else:
            colors[n] = palette[k % len(palette)]
            k += 1
    return colors


def mark_better(ax, higher_is_better):
    """在子图角落画红色箭头 + 英文 'better'，标注数值改进方向（红=好）。"""
    if higher_is_better:
        y0, y1 = 0.80, 0.93
    else:
        y0, y1 = 0.93, 0.80
    ax.annotate("", xy=(0.965, y1), xytext=(0.965, y0),
                xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#d62728", lw=1.5))
    ax.text(0.945, (y0 + y1) / 2.0, "better", color="#d62728",
            transform=ax.transAxes, ha="right", va="center", fontsize=8)


def style_axes(ax):
    ax.set_facecolor("white")
    for sp in ax.spines.values():
        sp.set_color("black")
        sp.set_linewidth(0.8)
    ax.tick_params(colors="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.title.set_color("black")
    ax.grid(True, color="#dddddd", linewidth=0.6, linestyle="-")


def write_json(result, prefix):
    with open(prefix + "_metrics.json", "w", encoding="utf-8") as f:
        json.dump(sanitize(result), f, indent=2, ensure_ascii=False)


def write_table(result, prefix):
    methods = result["methods"]
    meta = result["case"]
    lines = []
    lines.append("# 相位改进评估汇总")
    lines.append("")
    lines.append(f"- case: `{meta['file']}`")
    lines.append(f"- 空间群: {meta['space_group']}  晶胞: "
                 + ", ".join(f"{x:.3f}" for x in meta["cell"]))
    lines.append(f"- 分辨率: d_min = {fmt(meta['d_min'], 3)} Å, "
                 f"d_max = {fmt(meta['d_max'], 3)} Å")
    lines.append(f"- 反射数: {meta['n_reflections']}  "
                 f"(work {meta['n_work']} / free {meta['n_free']})")
    lines.append(f"- 网格: {result['grid']}  分辨率壳数: {meta['n_shells']}")
    lines.append("")
    lines.append("> 改进方向：wMPE / MPE / wMPE_work / wMPE_free **越小越好**；"
                 "mean cosΔφ / σA C / map CC / mean FOM **越大越好**（图中红色=更好）。")
    lines.append("> 所有指标均按 HKL dict join 对齐；权重 w = m_init·Fo（无 FOM 列时用 Fo）。")
    lines.append("")
    header = ("| 方法 | wMPE | MPE | mean cosΔφ | σA C | map CC 全胞 | "
              "map CC 高密度区 | wMPE_work | wMPE_free | mean FOM |")
    sep = "|---|---|---|---|---|---|---|---|---|---|"
    lines.append(header)
    lines.append(sep)
    if not methods:
        lines.append("| _(无有效方法)_ | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A | N/A |")
    else:
        for name, r in methods.items():
            g = r["global_"]
            cc = r["map_CC"]
            wf = r["work_free"]
            lines.append(
                f"| {name} | {fmt(g['wMPE'])} | {fmt(g['MPE'])} | "
                f"{fmt(g['mean_cos_dphi'])} | {fmt(g['sigmaA_like_C'])} | "
                f"{fmt(cc['whole_cell'])} | {fmt(cc['high_density'])} | "
                f"{fmt(wf['wMPE_work'])} | {fmt(wf['wMPE_free'])} | "
                f"{fmt(g['mean_FOM'])} |")
    skipped = result.get("skipped_methods", [])
    if skipped:
        lines.append("")
        lines.append("## 已跳过的方法")
        for s in skipped:
            lines.append(f"- `{s['name']}`: {s['reason']}")
    if result.get("warnings"):
        lines.append("")
        lines.append("## 警告")
        for w in result["warnings"]:
            lines.append(f"- {w}")
    lines.append("")
    with open(prefix + "_table.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def plot_shells(result, prefix, colors, best):
    methods = result["methods"]
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    style_axes(ax)
    if methods:
        for name, r in methods.items():
            xs = [s["d_center"] for s in r["shells"] if s["d_center"] is not None]
            ys = [s["wMPE"] for s in r["shells"] if s["d_center"] is not None]
            ax.plot(xs, ys, "o-", color=colors[name],
                    lw=2.6 if name == best else 1.5,
                    ms=6 if name == best else 4,
                    alpha=1.0 if name == best else 0.85, label=name)
        ax.invert_xaxis()                          # 低分辨(大 d) -> 高分辨(小 d)
        ax.legend(frameon=True, facecolor="white", edgecolor="#cccccc",
                  fontsize=9, loc="best")
    else:
        ax.text(0.5, 0.5, "No valid methods", ha="center", va="center",
                transform=ax.transAxes)
    ax.set_xlabel("Resolution d (A)   [high resolution ->]")
    ax.set_ylabel("wMPE (deg)")
    ax.set_title("wMPE vs resolution shell (10 equal-volume shells)")
    mark_better(ax, higher_is_better=False)
    fig.tight_layout()
    fig.savefig(prefix + "_wmpe_shells.png")
    plt.close(fig)


def plot_summary(result, prefix, colors, best):
    methods = result["methods"]
    names = list(methods.keys())
    panels = [
        ("wMPE (deg)", "wMPE", False,
         lambda r: r["global_"]["wMPE"]),
        ("map CC (whole cell)", "map CC", True,
         lambda r: r["map_CC"]["whole_cell"]),
        ("mean cos(dphi)", "mean cos", True,
         lambda r: r["global_"]["mean_cos_dphi"]),
        ("wMPE_free (deg)", "wMPE_free", False,
         lambda r: r["work_free"]["wMPE_free"]),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, (title, short, higher, getter) in zip(axes.ravel(), panels):
        style_axes(ax)
        vals = []
        for n in names:
            v = getter(methods[n])
            vals.append(v if (v is not None and math.isfinite(v)) else None)
        xpos = np.arange(len(names)) if names else np.array([])
        finite_vals = [v for v in vals if v is not None]
        if finite_vals:
            bar_colors = [colors[n] for n in names]
            heights = [0.0 if v is None else v for v in vals]
            ax.bar(xpos, heights, color=bar_colors, edgecolor="#333333", linewidth=0.6)
            lo = min(finite_vals)
            hi = max(finite_vals)
            span = hi - lo if hi > lo else (abs(hi) if hi != 0 else 1.0)
            ax.set_ylim(min(0, lo - 0.15 * span), hi + 0.25 * span)
        else:
            ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                    transform=ax.transAxes, color="#888888")
        if names:
            ax.set_xticks(xpos)
            ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
        ax.set_title(title)
        ax.set_ylabel(short)
        mark_better(ax, higher_is_better=higher)
    fig.suptitle("Method comparison  (red = best / improvement direction)",
                 color="black", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(prefix + "_summary.png")
    plt.close(fig)


def plot_phase_error(result, prefix, colors, best):
    methods = result["methods"]
    names = list(methods.keys())[:4]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    axes = axes.ravel()
    for i, ax in enumerate(axes):
        style_axes(ax)
        if i >= len(names):
            ax.axis("off")
            continue
        n = names[i]
        h = methods[n]["phase_error_hist"]
        edges = np.asarray(h["bin_edges"], dtype=float)
        counts = np.asarray(h["counts"], dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        width = edges[1] - edges[0]
        ax.bar(centers, counts, width=width * 0.92, color=colors[n],
               edgecolor="#333333", linewidth=0.4)
        ax.axvline(0.0, color="#d62728", linestyle="--", linewidth=1.0)
        wm = methods[n]["global_"]["wMPE"]
        ax.set_title(f"{n}   wMPE={fmt(wm, 2)} deg")
        ax.set_xlabel("phase error dphi (deg)")
        ax.set_ylabel("count")
        ax.set_xlim(-180, 180)
    fig.suptitle("Phase error distribution  (red dashed = zero error)",
                 color="black", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(prefix + "_phase_error.png")
    plt.close(fig)


def write_plots(result, prefix):
    methods = result["methods"]
    names = list(methods.keys())
    best = None
    cand = [(r["global_"]["wMPE"], n) for n, r in methods.items()
            if r["global_"]["wMPE"] is not None and math.isfinite(r["global_"]["wMPE"])]
    if cand:
        best = min(cand)[1]
    colors = method_colors(names, best)
    plot_shells(result, prefix, colors, best)
    plot_summary(result, prefix, colors, best)
    plot_phase_error(result, prefix, colors, best)


def build_result(case, shape, results, skipped):
    meta = dict(case["meta"])
    meta["n_shells"] = 10
    meta["shells"] = case["shells_meta"]
    meta["grid"] = list(shape)
    return dict(
        case=meta,
        grid=list(shape),
        methods=results,
        skipped_methods=skipped,
        warnings=list(WARNINGS),
        notes=dict(
            wMPE="Sum(w*|wrap(dphi)|)/Sum(w), w = m_init*Fo (Fo if no FOM)",
            MPE="unweighted mean |dphi|",
            mean_cos_dphi="Sum(w*cos(dphi))/Sum(w)  [true-FOM proxy]",
            sigmaA_like_C="|Sum(w*Fo^2*exp(i*dphi))|/Sum(w*Fo^2)",
            map_CC="Fo*exp(i*phi_true) vs Fo*exp(i*phi_out); "
                   "gemmi.get_f_phi_on_grid + ifftn(arr.conj()).real*(N/V)",
            work_free="free_col == 0 is the free set",
            phase_unit="max(|phi|)<=3.15 -> radians (converted to deg), else deg; wrap (-180,180]",
        ),
    )


# ----------------------------------------------------------------------------
# 主流程
# ----------------------------------------------------------------------------
def evaluate(args):
    for key in ("case", "truth_col", "fo_col", "methods", "out"):
        if not getattr(args, key, None):
            fail(f"缺少必需参数 --{key.replace('_', '-')}")
    case = load_case(args)
    shape = grid_shape(case)
    N = int(np.prod(shape))
    ctx = (shape, N, case["sg"], case["cell"], case["volume"])
    print(f"[evaluate] case={case['file']} SG={case['meta']['space_group']} "
          f"grid={shape} nref={case['meta']['n_reflections']}")

    specs = parse_methods(args.methods)
    results = {}
    skipped = []
    for spec in specs:
        m = load_method(spec)
        if m is None:
            skipped.append(dict(name=spec["name"],
                                reason=f"MTZ 或相位列缺失/不可读 ({spec['path']}:{spec['phi']})"))
            continue
        r = evaluate_method(case, m, ctx, args)
        if r is None:
            skipped.append(dict(name=spec["name"], reason="无公共/有效反射"))
            continue
        results[spec["name"]] = r
        print(f"[evaluate]   {spec['name']:<12s} wMPE={fmt(r['global_']['wMPE'],3):>8s} "
              f"CC={fmt(r['map_CC']['whole_cell'],4):>8s} "
              f"CC_hi={fmt(r['map_CC']['high_density'],4):>8s} "
              f"unit={r['phase_unit']}")

    if not results:
        warn("没有任何方法成功评估；仍将生成 5 个输出文件（内容为空/标注）")

    result = build_result(case, shape, results, skipped)
    prefix = args.out
    outdir = os.path.dirname(os.path.abspath(prefix))
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)
    write_json(result, prefix)
    write_table(result, prefix)
    write_plots(result, prefix)
    print(f"[evaluate] 输出: {prefix}_metrics.json / _table.md / "
          f"_wmpe_shells.png / _summary.png / _phase_error.png")
    print(f"[evaluate] 警告数: {len(result['warnings'])}, 跳过方法: {len(skipped)}")
    return result


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="评估多种密度修饰程序的相位改进效果（gemmi + numpy + matplotlib）")
    p.add_argument("--case", help="case MTZ（含真值相位/观测振幅/初始FOM/free）")
    p.add_argument("--truth-col", default=None, help="case 中真值相位列，如 PHIC（可用 | 列多个候选）")
    p.add_argument("--truth-fcol", default=None, help="可选：case 中真值振幅列，如 FC")
    p.add_argument("--fo-col", default=None, help="观测振幅列，如 FOBS / FPHASED")
    p.add_argument("--fom-col", default=None, help="初始 FOM 列（权重 w = m_init*Fo），如 FOM")
    p.add_argument("--free-col", default=None, help="free 标记列，如 FreeR / FreeR_flag")
    p.add_argument("--free-marker", type=float, default=0,
                   help="free 标记值（CCP4 FreeR_flag 用 0；二值 FreeR 列常用 1）")
    p.add_argument("--methods", default=None,
                   help='"name1=path1.mtz:PHI_COL[:FOM_COL],name2=..."')
    p.add_argument("--out", default=None, help="输出前缀")
    p.add_argument("--selftest", action="store_true", help="合成 MTZ 自测（产物只写 /tmp）")
    return p.parse_args(argv)


# ----------------------------------------------------------------------------
# 自测：/tmp 合成 MTZ，验证 5 个输出 + 真值方法 wMPE≈0 / map CC≈1
# ----------------------------------------------------------------------------
def _write_mtz(path, sg, cell, columns, arrays):
    m = gemmi.Mtz(with_base=True)
    m.spacegroup = sg
    m.cell = cell
    m.add_dataset("selftest")
    for lab, typ in columns:
        m.add_column(lab, typ)
    m.set_data(np.column_stack(arrays).astype(np.float32))
    m.write_to_file(path)


def selftest():
    keep = os.environ.get("EVAL_SELFTEST_KEEP") == "1"
    d = tempfile.mkdtemp(prefix="evaluate_selftest_")
    try:
        rng = np.random.default_rng(2026)
        cell = gemmi.UnitCell(30.0, 40.0, 50.0, 90.0, 90.0, 90.0)
        sg = gemmi.find_spacegroup_by_name("P 1")

        # 由若干高斯"原子"合成真实密度，FFT 得真值结构因子（保证密度是实的、有峰）
        nx, ny, nz = 32, 40, 48
        rho = np.zeros((nx, ny, nz))
        x = np.arange(nx)[:, None, None]
        y = np.arange(ny)[None, :, None]
        z = np.arange(nz)[None, None, :]
        for _ in range(12):
            cx, cy, cz = rng.uniform(0, nx), rng.uniform(0, ny), rng.uniform(0, nz)
            amp = rng.uniform(5, 20)
            dx = np.minimum(np.abs(x - cx), nx - np.abs(x - cx))
            dy = np.minimum(np.abs(y - cy), ny - np.abs(y - cy))
            dz = np.minimum(np.abs(z - cz), nz - np.abs(z - cz))
            rho += amp * np.exp(-(dx ** 2 + dy ** 2 + dz ** 2) / (2.0 * 3.0 ** 2))
        Fgrid = np.fft.fftn(rho)

        H, K, L, Fo, PHI = [], [], [], [], []
        for h in range(-6, 7):
            for k in range(-6, 7):
                for l in range(-6, 7):
                    if h == 0 and k == 0 and l == 0:
                        continue
                    if not (l > 0 or (l == 0 and k > 0) or (l == 0 and k == 0 and h > 0)):
                        continue                       # P1 ASU（每对 Friedel 取一个）
                    fval = Fgrid[h % nx, k % ny, l % nz]
                    H.append(h); K.append(k); L.append(l)
                    Fo.append(abs(fval))
                    PHI.append(np.degrees(np.angle(fval)))
        H = np.asarray(H, dtype=int)
        K = np.asarray(K, dtype=int)
        L = np.asarray(L, dtype=int)
        Fo = np.asarray(Fo, dtype=float)
        Fo = Fo / max(Fo.max(), 1e-9) * 120.0        # 归一到合理振幅
        PHI = np.asarray(PHI, dtype=float)
        n = len(H)
        fom = rng.uniform(0.35, 0.95, n)
        free = np.where(rng.random(n) < 0.07, 0, rng.integers(1, 20, n)).astype(float)
        FC = Fo * rng.uniform(0.8, 1.2, n)
        PHIB = PHI + rng.normal(0, 55, n)            # 初始实验相位（仅占位）

        case_path = os.path.join(d, "case.mtz")
        _write_mtz(case_path, sg, cell,
                   [("FOBS", "F"), ("FOM", "W"), ("FreeR", "I"),
                    ("PHIC", "P"), ("FC", "F"), ("PHIB", "P")],
                   [H, K, L, Fo, fom, free, PHI, FC, PHIB])

        # 方法 1：真值，但相位用【弧度】且行顺序【打乱】——测单位判定 + HKL dict join
        perm = rng.permutation(n)
        p1 = os.path.join(d, "truth.mtz")
        _write_mtz(p1, sg, cell, [("PHIT", "P"), ("FOMT", "W")],
                   [H[perm], K[perm], L[perm], np.radians(PHI)[perm], fom[perm]])

        # 方法 2/3：加噪相位，带 / 不带 FOM 列
        p2 = os.path.join(d, "err20.mtz")
        _write_mtz(p2, sg, cell, [("PHIE20", "P"), ("FOME", "W")],
                   [H, K, L, PHI + rng.normal(0, 20, n), rng.uniform(0.4, 0.8, n)])
        p3 = os.path.join(d, "err60.mtz")
        _write_mtz(p3, sg, cell, [("PHIE60", "P")],
                   [H, K, L, PHI + rng.normal(0, 60, n)])

        # 方法 4：随机相位
        p4 = os.path.join(d, "random.mtz")
        _write_mtz(p4, sg, cell, [("PHIR", "P")],
                   [H, K, L, rng.uniform(-180, 180, n)])

        # 方法 5：故意缺失相位列 -> 验证优雅跳过
        p5 = os.path.join(d, "broken.mtz")
        _write_mtz(p5, sg, cell, [("PHIX", "P")], [H, K, L, PHI])

        methods = (
            f"truth={p1}:PHIT:FOMT,"
            f"err20={p2}:PHIE20:FOME,"
            f"err60={p3}:PHIE60,"
            f"random={p4}:PHIR,"
            f"broken={p5}:NOPE"
        )
        out = os.path.join(d, "selftest")
        argv = ["--case", case_path, "--truth-col", "PHIC", "--truth-fcol", "FC",
                "--fo-col", "FOBS", "--fom-col", "FOM", "--free-col", "FreeR",
                "--methods", methods, "--out", out]
        args = parse_args(argv)
        evaluate(args)

        # 校验 5 个输出文件
        files = [out + s for s in ("_metrics.json", "_table.md",
                                   "_wmpe_shells.png", "_summary.png",
                                   "_phase_error.png")]
        for f in files:
            assert os.path.isfile(f) and os.path.getsize(f) > 0, f"缺少/空输出: {f}"

        with open(out + "_metrics.json", encoding="utf-8") as f:
            res = json.load(f)
        assert "truth" in res["methods"], "真值方法未被评估"
        assert len(res["methods"]) == 4, f"应有 4 个成功方法, 实得 {len(res['methods'])}"
        assert any(s["name"] == "broken" for s in res["skipped_methods"]), \
            "缺失列的方法未被优雅跳过"
        t = res["methods"]["truth"]
        wmpe = t["global_"]["wMPE"]
        cc = t["map_CC"]["whole_cell"]
        cchi = t["map_CC"]["high_density"]

        print("\n===== SELFTEST 结果 =====")
        for name, r in res["methods"].items():
            g = r["global_"]
            print(f"  {name:<8s} wMPE={fmt(g['wMPE'],3):>8s}  MPE={fmt(g['MPE'],3):>8s}  "
                  f"cosDphi={fmt(g['mean_cos_dphi'],3):>7s}  "
                  f"CC={fmt(r['map_CC']['whole_cell'],4):>7s}  "
                  f"CC_hi={fmt(r['map_CC']['high_density'],4):>7s}  "
                  f"wMPE_free={fmt(r['work_free']['wMPE_free'],3):>8s}  "
                  f"unit={r['phase_unit']}")
        print(f"  skipped: {[s['name'] for s in res['skipped_methods']]}")
        print(f"  warnings: {len(res['warnings'])}")
        print(f"  grid: {res['grid']}, shells: {res['case']['n_shells']}")

        assert wmpe is not None and wmpe < 0.5, f"真值方法 wMPE 应≈0, 实得 {wmpe}"
        assert cc is not None and cc > 0.999, f"真值方法 map CC 应≈1, 实得 {cc}"
        assert cchi is not None and cchi > 0.999, f"真值方法高密度区 CC 应≈1, 实得 {cchi}"
        print(f"\nSELFTEST OK: truth wMPE={wmpe:.6f} deg, "
              f"map CC={cc:.6f}, CC_high={cchi:.6f}")
        print("输出目录:", d)
        return 0
    finally:
        if not keep:
            import shutil
            shutil.rmtree(d, ignore_errors=True)


def main(argv=None):
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    evaluate(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
