#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import math
import numpy as np

import matplotlib
matplotlib.use("Agg")  # 安全 headless
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter


class DiffChartGenerator:
    """
    專門畫「diff 圖」的產圖器。
    diff_mode:
      - "norm": (ours - other) / (UB - other)
      - "gap" : (UB - other) - (UB - ours) == (ours - other)
    只會畫「other」演算法（ours 與 UB 不會畫），
    並在輸出檔名自動加上 `_diff` 後綴。
    """
    # ---- 你可以依需要調整這些樣式/表 ---- #
    THEME = {
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "axes.labelsize": 20,
        "axes.titlesize": 20,
        "font.family": "DejaVu Serif",
        "mathtext.fontset": "dejavuserif",
        "text.usetex": False,
    }

    # 圖例錨點（仿你的設定）
    RIGHT_TOP = (0.7, 0.87)
    LEFT_TOP  = (0.3, 0.87)
    RIGHT_DOWN= (0.7, 0.13)
    LEFT_DOWN = (0.3, 0.13)

    BBOX_POS = {
        'fidelity_gain':{
            'request_cnt': LEFT_TOP,
            'tao': RIGHT_TOP,
            'time_limit': LEFT_TOP,
            'avg_memory': LEFT_TOP,
            'min_fidelity': LEFT_TOP,
            'fidelity_threshold': RIGHT_TOP,
            'swap_prob': LEFT_TOP,
            'entangle_time': RIGHT_DOWN,
            'entangle_prob': LEFT_DOWN,
            'Zmin': LEFT_TOP,
            'time_eta': LEFT_TOP,
            'bucket_eps': LEFT_TOP
        },
        'succ_request_cnt':{
            'request_cnt': LEFT_TOP,
            'tao': RIGHT_TOP,
            'time_limit': LEFT_TOP,
            'avg_memory': LEFT_TOP,
            'min_fidelity': LEFT_TOP,
            'fidelity_threshold': RIGHT_TOP,
            'swap_prob': LEFT_TOP,
            'entangle_time': RIGHT_DOWN,
            'entangle_prob': LEFT_DOWN,
            'Zmin': LEFT_TOP,
            'time_eta': LEFT_TOP,
            'bucket_eps': LEFT_TOP
        },
        'pure_fidelity':{
            'request_cnt': LEFT_TOP,
            'tao': RIGHT_TOP,
            'time_limit': LEFT_TOP,
            'avg_memory': LEFT_TOP,
            'min_fidelity': LEFT_TOP,
            'fidelity_threshold': RIGHT_TOP,
            'swap_prob': LEFT_TOP,
            'entangle_time': RIGHT_DOWN,
            'entangle_prob': LEFT_DOWN,
            'Zmin': LEFT_TOP,
            'time_eta': LEFT_TOP,
            'bucket_eps': LEFT_TOP
        }
    }

    # x/y 顯示名（沿用你的 pretty 命名）
    PRETTY = {
        "request_cnt": r"$\#$Requests",
        "time_limit": "$|T|$",
        "avg_memory": "Average Memory Limit",
        "tao": r"$\it{\tau}$",
        "swap_prob": "Swapping Probability",
        "fidelity_gain": "Expected Fidelity Sum",
        "succ_request_cnt": r"$\#$Accepted Requests",
        "pure_fidelity": "Average Final Fidelity",
        "fidelity_threshold": "Fidelity Threshold",
        "min_fidelity": "Minimum Initial Fidelity",
        "entangle_time": "Entangling Time",
        "entangle_prob": "Entangling Probability",
        "Zmin": r"$Z_{min}$",
        "time_eta": r"$\eta$",
        "bucket_eps": r"bucket $\epsilon$",
    }

    # 顏色/標記池（diff 會只畫 other，所以先準備足夠多）
    COLORS  = ["#FF0000","#00AA00","#0044FF","#000000","#900321",
               "#FF00FF","#00CCCC","#FF8800","#8A2BE2","#2E8B57"]
    MARKERS = ['s','v','o','^','x','1','D','P','*','h']

    # ---- 主要 API ---- #
    def __init__(self,
                 our_algo="G-FNPR",
                 ub_algo="G-UB",
                 diff_mode="norm",
                 y_fixed=None,  # 例如 (-1, 1, 0.2)；None 表示自動
                 label_every=None):  # 主要刻度每幾格標註一次（None 用自動/1）
        self.our_algo = our_algo
        self.ub_algo = ub_algo
        assert diff_mode in ("norm", "gap")
        self.diff_mode = diff_mode
        self.y_fixed = y_fixed
        self.label_every = label_every

        matplotlib.rcParams.update(self.THEME)

    def render(self, dataName, Xlabel_key, Ylabel_key, out_dir="../data/pdf"):
        """
        讀 ../data/ans/<dataName>，畫 diff 圖，輸出到 out_dir，
        檔名自動加 `_diff`。
        """
        filename = os.path.join("../data/ans", dataName)
        if not os.path.exists(filename):
            print(f"[WARN] file doesn't exist: {os.path.abspath(filename)}")
            return

        with open(filename, 'r', encoding='utf-8') as f:
            raw_lines = f.readlines()
        print("start generate (diff)", filename)

        # 讀資料（穩健）
        x, y, algo_names = self._read_ans(raw_lines)

        if x is None:  # 無有效資料
            print(f"[ERROR] no valid data rows in {filename}")
            return

        # X 軸工程縮放（仿你的行為）
        Xpow = -3 if Xlabel_key in ("tao","entangle_time","entangle_prob") else 0
        Xdiv = float(10**Xpow)
        for i in range(len(x)):
            if Xpow != 0:
                xv = float(x[i]) / Xdiv
                if abs(xv - int(xv)) <= 1e-6:
                    xv = int(round(xv * 10) / 10)
                x[i] = xv

        # 做 diff 轉換：產出 new_y 與對應的 other 演算法名稱
        new_y, other_names = self._diff_transform(y, algo_names)
        if not new_y:
            print(f"[WARN] nothing to plot after diff in {filename}")
            return

        # y 範圍
        minv = min(min(row) for row in new_y)
        maxv = max(max(row) for row in new_y)
        if self.y_fixed is None:
            Ystart, Yend, Yinterval = self._auto_y_range(minv, maxv, target_ticks=7, padding=0.05)
            label_step_factor = max(1, int(math.ceil((Yend - Ystart) / max(1e-12, Yinterval) / 7)))
        else:
            Ystart, Yend, Yinterval = self.y_fixed
            label_step_factor = 1
        if self.label_every is not None:
            label_step_factor = max(1, int(self.label_every))

        # 圖面
        fig, ax1 = plt.subplots(figsize=(6.8,5.4), dpi=600, constrained_layout=True)
        ax1.tick_params(direction="in", bottom=True, top=True, left=True, right=True, pad=10)

        # 畫線
        colors  = (self.COLORS  * ((len(new_y)//len(self.COLORS))+1))[:len(new_y)]
        markers = (self.MARKERS * ((len(new_y)//len(self.MARKERS))+1))[:len(new_y)]

        for i, row in enumerate(new_y):
            ax1.plot(
                range(len(x)), row,
                color=colors[i], lw=1, ls="-",
                marker=markers[i], markersize=8,
                markerfacecolor='None', markeredgewidth=2, zorder=-i
            )

        # 圖例
        bbox = self.BBOX_POS.get(Ylabel_key, {}).get(Xlabel_key, self.LEFT_TOP)
        plt.legend(
            other_names,
            loc='center',
            bbox_to_anchor=bbox,
            prop={"size": 16, "family": "DejaVu Serif"},
            frameon=False, handletextpad=0.2, handlelength=1,
            labelspacing=0.2, columnspacing=0.6, ncol=2, facecolor='None',
        )

        # 標籤
        unit = " s" if Xlabel_key in ("tao", "entangle_time") else ""
        Xlabel_text = self._pretty(Xlabel_key) + self._scale_suffix(Xpow, unit)
        baseY = self._pretty(Ylabel_key)
        if self.diff_mode == "norm":
            Ylabel_text = baseY + " (diff = (ours-other)/(UB-other))"
        else:
            Ylabel_text = baseY + " (diff = ours - other)"

        # y locators
        major_step = max(1e-12, Yinterval * label_step_factor)
        minor_step = max(1e-12, Yinterval)
        span = max(1e-12, Yend - Ystart)
        approx_labels = span / major_step
        if approx_labels > 1000:
            factor = max(1, int(math.ceil(approx_labels / 12)))
            major_step *= factor

        ax1.set_ylim(Ystart, Yend)
        ax1.yaxis.set_major_locator(MultipleLocator(major_step))
        ax1.yaxis.set_minor_locator(MultipleLocator(minor_step))
        ax1.tick_params(axis='y', which='major', labelsize=20)
        ax1.tick_params(axis='y', which='minor', length=4, labelsize=0)

        def _fmt_y(val, pos):
            if abs(val) >= 100 or abs(val) == int(val):
                return f"{val:.0f}"
            elif abs(val) >= 10:
                return f"{val:.1f}"
            else:
                return f"{val:.2f}"
        ax1.yaxis.set_major_formatter(FuncFormatter(_fmt_y))

        plt.xticks(ticks=range(len(x)), labels=x, fontsize=20)
        plt.ylabel(Ylabel_text, fontsize=22)
        plt.xlabel(Xlabel_text, fontsize=22, labelpad=16)
        ax1.yaxis.set_label_coords(-0.13, 0.5)
        ax1.xaxis.set_label_coords(0.45, -0.13)

        # 輸出
        os.makedirs(out_dir, exist_ok=True)
        base = dataName[:-4].replace('#','') + "_diff"
        plt.savefig(os.path.join(out_dir, f"Fie4_{base}.eps"), bbox_inches="tight", pad_inches=0.3)
        plt.savefig(os.path.join(out_dir, f"{base}.jpg"),     bbox_inches="tight", pad_inches=0.3)
        plt.close()

    # --------- internal utils ---------- #
    def _read_ans(self, lines):
        """
        讀 .ans：每行 'x y1 y2 ... yk'
        回傳 x(list of str/float), y(list[list[float]] with shape [k][n]), algo_names(依你資料順序推測)
        """
        x = []
        rows = []
        num_cols = None

        for raw in lines:
            s = raw.strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            if num_cols is None:
                num_cols = len(parts)
            elif len(parts) != num_cols:
                print(f"[WARN] skip row (expected {num_cols} cols, got {len(parts)}): {s}")
                continue
            x.append(parts[0])
            rows.append([float(v) for v in parts[1:]])

        if num_cols is None or not rows:
            return None, None, None

        data = np.array(rows, dtype=float)  # shape [n][k]
        y = data.T.tolist()                 # shape [k][n]

        # 你的演算法名稱固定順序（比 y 條數多時截取）
        per_algo_name = ["G-Werner","G-FNPR","G-UB","G-FLTO","G-Nesting","G-Linear","G-ASAP"]
        k = len(y)
        algo_names = per_algo_name[:k]
        return x, y, algo_names

    def _diff_transform(self, y, algo_names):
        """
        y: [algo][n]
        回傳 new_y(list of rows) 與 other_names（不含 ours/UB）
        """
        name_to_idx = {name:i for i, name in enumerate(algo_names)}
        if self.our_algo not in name_to_idx or self.ub_algo not in name_to_idx:
            print(f"[ERROR] cannot find algo index: ours={self.our_algo}, ub={self.ub_algo}, available={algo_names}")
            return [], []

        io = name_to_idx[self.our_algo]
        iu = name_to_idx[self.ub_algo]
        y_ours = np.array(y[io])
        y_ub   = np.array(y[iu])

        new_y = []
        names = []
        for i, name in enumerate(algo_names):
            if i in (io, iu):
                continue
            other = np.array(y[i])
            if self.diff_mode == "norm":
                denom = (y_ub - other)
                denom = np.where(np.abs(denom) < 1e-12, np.nan, denom)  # 避免除零
                vals = (y_ours - other) / denom
            else:
                vals = (y_ours - other)
            new_y.append(vals.tolist())
            names.append(name)
        return new_y, names

    def _auto_y_range(self, ymin, ymax, target_ticks=6, padding=0.05):
        if ymin == math.inf or ymax == -math.inf:
            return (0.0, 1.0, 0.2)
        if ymin == ymax:
            if ymin == 0:
                return (0.0, 1.0, 0.2)
            pad = abs(ymin) * max(padding, 0.02)
            ymin, ymax = ymin - pad, ymax + pad
        span = ymax - ymin
        rough = span / max(2, target_ticks - 1)
        step = self._nice_step(rough)
        start = math.floor(ymin / step) * step
        end   = math.ceil (ymax / step) * step
        return (start, end, step)

    def _nice_step(self, span):
        if span <= 0:
            return 1.0
        exp = math.floor(math.log10(span))
        base = span / (10**exp)
        if base <= 1.5: m = 1
        elif base <= 3: m = 2
        elif base <= 7: m = 5
        else: m = 10
        return m * (10**exp)

    def _pretty(self, key):
        return self.PRETTY.get(key, key)

    def _scale_suffix(self, multiple, unit):
        if multiple == 0:
            return ""
        return f" ($10^{{{multiple}}}$" + (unit if unit else "") + ")"


if __name__ == "__main__":
    # ===== 你可以在這裡調整要跑的檔名模式 =====
    Xlabels = ["entangle_prob","swap_prob","request_cnt", "time_limit", "tao",
               "avg_memory", "fidelity_threshold", "min_fidelity", "entangle_time",
               "Zmin", "time_eta", "bucket_eps"]
    Ylabels = ["fidelity_gain", "succ_request_cnt","pure_fidelity"]
    PathNames = ["Greedy", "QCAST", "REPS"]  # 依你的資料前綴

    # ===== diff 圖設定 =====
    DIFF_MODE = "norm"          # "norm" 或 "gap"
    OUR_ALGO  = "G-FNPR"
    UB_ALGO   = "G-UB"
    # 想固定 y 範圍（例如常見的 -1~1）：把 y_fixed 換成 (-1.0, 1.0, 0.2)
    Y_FIXED = None
    LABEL_EVERY = None  # 例如 2 -> 每兩格顯示一次主刻度標籤

    gen = DiffChartGenerator(
        our_algo=OUR_ALGO,
        ub_algo=UB_ALGO,
        diff_mode=DIFF_MODE,
        y_fixed=Y_FIXED,
        label_every=LABEL_EVERY
    )

    for Path in PathNames:
        for X in Xlabels:
            for Y in Ylabels:
                dataFileName = f"{Path}_{X}_{Y}.ans"
                gen.render(dataFileName, X, Y, out_dir="../data/pdf")
