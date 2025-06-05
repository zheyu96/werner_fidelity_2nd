import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

class ChartGenerator:
    def __init__(self, dataName, _Xlabel, _Ylabel, Xlabel, Ylabel, Xpow, Ypow, Ystart, Yend, Yinterval):
        filename = '../data/ans/' + dataName

        if not os.path.exists(filename):
            print("file doesn't exist")
            return
        
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
        print("start generate", filename)
        
        fontsize = 22
        Xlabel_fontsize = fontsize
        Ylabel_fontsize = fontsize
        Xticks_fontsize = fontsize
        Yticks_fontsize = fontsize
            
        andy_theme = {
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
            "axes.labelsize": 20,
            "axes.titlesize": 20,
            "mathtext.fontset": "custom",
            "font.family": "Times New Roman",
            "mathtext.default": "default",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.cal": "Times New Roman:italic",
            "text.usetex": True,
        }
        
        matplotlib.rcParams.update(andy_theme)        
        fig, ax1 = plt.subplots(figsize=(6, 5), dpi=600)
   
        ax1.tick_params(direction="in", bottom=True, top=True, left=True, right=True, pad=10)

        RIGHT_TOP = (0.7, 0.87)
        LEFT_TOP = (0.3, 0.87)
        bbox_pos_settings = {
            'fidelity_gain_log':{
                'request_cnt': LEFT_TOP,
                'tao': RIGHT_TOP,
                'time_limit': LEFT_TOP,
                'avg_memory': LEFT_TOP,
                'min_fidelity': LEFT_TOP,
                'fidelity_threshold': RIGHT_TOP,
                'swap_prob': LEFT_TOP,
                'entangle_time': RIGHT_TOP,
                'entangle_prob': LEFT_TOP
            },
            'succ_request_cnt_log':{
                'request_cnt': LEFT_TOP,
                'tao': RIGHT_TOP,
                'time_limit': LEFT_TOP,
                'avg_memory': LEFT_TOP,
                'min_fidelity': LEFT_TOP,
                'fidelity_threshold': RIGHT_TOP,
                'swap_prob': LEFT_TOP,
                'entangle_time': RIGHT_TOP,
                'entangle_prob': LEFT_TOP
            },
        }

        Y_interval_settings = {
            'fidelity_gain_log':{
                'entangle_prob': (-10, 2, 3)
            },
            'succ_request_cnt_log':{
                'entangle_prob': (-10, 2, 3)
            }
        }

        x = []
        y = []
        num_of_data = 0
        for line in lines:
            line = line.strip()
            data = line.split(' ')
            num_of_line = len(data)
            num_of_data += 1

            for i in range(num_of_line):
                if i == 0:
                    x.append(data[i])
                else:
                    y.append(data[i])
        num_of_algo = num_of_line - 1
        y = np.array(y).reshape(-1, num_of_algo).T.tolist()

        max_data = 0
        min_data = math.inf
        Ydiv = float(10 ** Ypow)
        Xdiv = float(10 ** Xpow)
        
        for i in range(num_of_data):
            if Xpow != 0:
                x[i] = float(x[i]) / Xdiv
                x[i] = int(round(x[i] * 10) / 10)
            x[i] = int(x[i])
            x[i] = str("$10^{" + str((-x[i])) + "}$")
        for i in range(num_of_algo):
            for j in range(num_of_data):
                y[i][j] = float(y[i][j]) / Ydiv
                max_data = max(max_data, y[i][j])
                min_data = min(min_data, y[i][j])

        per_algo_name = ["G-FNPR", "G-UB", "G-FLTO", "G-Nesting", "G-Linear", "G-ASAP"]
        algo_name = []
        per = [0, 2, 1, 3, 4, 5]
        marker = ['s', 'v', 'o', '^', 'x', '1']
        color = [
            "#FF0000",  # red
            "#00FF00",  # lime
            "#0000FF",  # blue
            "#000000",  # black
            "#900321",  # brown
            "#FF00FF",  # green
        ]

        for idx in range(num_of_algo):
            i = per[idx]
            if _Ylabel == "succ_request_cnt_log" and per_algo_name[i][-2:] == "UB":
                continue
            ax1.semilogy(
                range(len(x)), 
                y[i],
                color = color[i], 
                lw = 1, 
                ls = "-", 
                marker = marker[i], 
                markersize = 8, 
                markerfacecolor = 'None', 
                markeredgewidth = 2, 
                zorder = -idx
            )
  
        for idx in range(num_of_algo):
            i = per[idx]
            if _Ylabel == "succ_request_cnt_log" and per_algo_name[i][-2:] == "UB":
                continue
            algo_name.append(per_algo_name[i])
        
        bbox_pos = bbox_pos_settings[_Ylabel][_Xlabel]
        (Ystart, Yend, Yinterval) = Y_interval_settings[_Ylabel][_Xlabel]
        Ystart /= Ydiv
        Yend /= Ydiv
        Yinterval /= Ydiv
        plt.legend(
            algo_name,
            loc = 'center',
            bbox_to_anchor = bbox_pos,
            prop = {"size": fontsize-6, "family": "Times New Roman"},
            frameon = False,
            handletextpad = 0.2,
            handlelength = 1,
            labelspacing = 0.2,
            columnspacing = 0.6,
            ncol = 2,
            facecolor = 'None',
        )
        
        unit = " s" if _Xlabel in ["tao", "entangle_time"] else ""
        Xlabel += self.genMultiName(Xpow, unit)
        Ylabel += self.genMultiName(Ypow, unit)
        plt.subplots_adjust(top = 0.97, left = 0.2, right = 0.975, bottom = 0.18)
        
        plt.yticks(np.logspace(Ystart, Yend, base = 10, num = 5), fontsize=Yticks_fontsize)
        ax1.tick_params(axis='y', pad=47)
        plt.xticks(ticks=range(len(x)), labels=x, fontsize=Xticks_fontsize)
        plt.ylabel(Ylabel, fontsize=Ylabel_fontsize)
        plt.xlabel(Xlabel, fontsize=Xlabel_fontsize, labelpad=500)
        ax1.yaxis.set_label_coords(-0.175, 0.5)
        ax1.xaxis.set_label_coords(0.45, -0.13)
        for label in ax1.get_yticklabels():
            label.set_horizontalalignment('left')
        pdfName = dataName[:-4].replace('#', '')
        plt.savefig('../data/pdf/Fie4_{}.eps'.format(pdfName)) 
        plt.savefig('../data/pdf/{}.jpg'.format(pdfName)) 
        plt.close()

    def genMultiName(self, multiple, unit):
        return "" if multiple == 0 else " ($" + "10" + "^{" + str(multiple) + "}$" + unit + ")"

if __name__ == "__main__":
    Xlabels = ["entangle_prob"]
    Ylabels = ["fidelity_gain_log", "succ_request_cnt_log"]
    PathNames = ["Greedy"]

    LabelsName = {
        "request_cnt": "$\#$Requests",
        "time_limit": "$|T|$",
        "avg_memory": "Average Memory Limit",
        "tao": "$\\it{\\tau}$",
        "swap_prob": "Swapping Probability",
        "fidelity_gain_log": "Expected Fidelity Sum",
        "succ_request_cnt_log": "$\#$Accepted Requests",
        "fidelity_threshold": "Fidelity Threshold",
        "min_fidelity": "Minimum Initial Fidelity",
        "entangle_time": "Entangling Time",
        "entangle_prob": "Entangling Probablity"
    }

    for Path in PathNames:
        for Xlabel in Xlabels:
            for Ylabel in Ylabels:
                dataFileName = Path + '_' + Xlabel + '_' + Ylabel + '.ans'
                ChartGenerator(dataFileName, Xlabel, Ylabel, LabelsName[Xlabel], LabelsName[Ylabel], 0, 0, 0, 2e-4, 4e-5)
