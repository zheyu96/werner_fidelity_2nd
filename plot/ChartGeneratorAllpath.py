import numpy as np
import math
import os
import latex
import matplotlib.pyplot as plt
import matplotlib.transforms
import matplotlib
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

class ChartGenerator:
    # data檔名 Y軸名稱 X軸名稱 Y軸要除多少(10的多少次方) Y軸起始座標 Y軸終止座標 Y軸座標間的間隔
    def __init__(self, pathNames, _Xlabel, _Ylabel, Xlabel, Ylabel, Xpow, Ypow, Ystart, Yend, Yinterval):

        
        fontsize = 22
        Xlabel_fontsize = fontsize
        Ylabel_fontsize = fontsize
        Xticks_fontsize = fontsize
        Yticks_fontsize = fontsize
            
        andy_theme = {
        # "axes.grid": True,
        # "grid.linestyle": "--",
        # "legend.framealpha": 1,
        # "legend.facecolor": "white",
        # "legend.shadow": True,
        # "legend.fontsize": 14,
        # "legend.title_fontsize": 16,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "axes.labelsize": 20,
        "axes.titlesize": 20,
        "mathtext.fontset": "custom",
        "font.family": "Times New Roman",
        "mathtext.default": "default",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.cal": "Times New Roman:italic",
        # "mathtext.fontset": "regular",
        # "figure.autolayout": True,
        "text.usetex": True,
        # "figure.dpi": 800
        }
        
        matplotlib.rcParams.update(andy_theme)        
        fig, ax1 = plt.subplots(figsize=(6, 5), dpi=600)
   
        # axis x, y ticks settings
        ax1.tick_params(direction="in")
        ax1.tick_params(bottom=True, top=True, left=True, right=True)
        ax1.tick_params(pad=10)

        total_algo_name = []

        for path in pathNames:
            filename = '../data/ans/' + path  + '_' + _Xlabel + '_' + _Ylabel + '.ans'
            print(filename)
            if not os.path.exists(filename):
                print("file doesn't exist")
                return
            
            with open(filename, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                
            print("start generate", filename)

            """
            Read Data to y

            """
            x = []
            y = []
            num_of_data = 0
            for line in lines:
                line = line.replace('\n','')
                if line[-1] == ' ':
                    line = line[0:-1]
                data = line.split(' ')
                num_of_line = len(data)
                num_of_data += 1

                for i in range(num_of_line):
                    if i == 0:
                        x.append(data[i])
                    else:
                        if Ylabel.endswith("(%)"):
                            y.append(str(float(data[i]) * 100))        
                        else:
                            y.append(data[i])
            num_of_algo = num_of_line - 1
            y = np.array(y).reshape(-1, num_of_algo).T.tolist() # transform 1-D list to 2-D list
            # print(x)
            # print(y)

            max_data = 0
            min_data = math.inf
            Ydiv = float(10 ** Ypow)
            Xdiv = float(10 ** Xpow)
            
            for i in range(num_of_data):
                x[i] = float(x[i]) / Xdiv

            for i in range(num_of_algo):
                for j in range(num_of_data):
                    y[i][j] = float(y[i][j]) / Ydiv
                    max_data = max(max_data, y[i][j])
                    min_data = min(min_data, y[i][j])

            """
            Start plotting

            """
            per_algo_name = ["FNPR", "UB", "FLTO", "Nesting", "Linear"]
            algo_name = []
            per = [0, 2, 1, 3, 4]
            marker = ['s', 'v', 'o', '^', 'x']
            color = [
                "#FF0000",  # red
                "#00FF00",  # lime
                "#0000FF",  # blue
                "#000000",  # black
                "#900321",  # brown
            ]

            line_style = {"Greedy": "-", "QCAST": ":", "REPS": "-."}
            for idx in range(num_of_algo):
                i = per[idx]
                if _Ylabel == "succ_request_cnt":
                    if per_algo_name[i] == "UB":
                        continue
                ax1.plot(
                    x, 
                    y[i], 
                    color = color[i], 
                    lw = 1.3, 
                    ls = line_style[path], 
                    marker = marker[i], 
                    markersize = 8, 
                    markerfacecolor = 'None', 
                    markeredgewidth = 2, 
                    zorder = -idx
                )
  
            for idx in range(num_of_algo):
                i = per[idx]
                if _Ylabel == "succ_request_cnt":
                    if per_algo_name[i] == "UB":
                        continue
                algo_name.append(path[0] + "-" + per_algo_name[i])    # adjust the order
            
            total_algo_name.extend(algo_name)
        
        # print(total_algo_name)
        
        # if "tau" in Xlabel:
        #     bbox_pos = (0.3, 0.15)
        # else:
        #     bbox_pos = (0.7, 0.15)

        leg = plt.legend(
            total_algo_name,
            loc = 'upper left',
            # bbox_to_anchor = (0.7, 0.8),
            # bbox_to_anchor = bbox_pos,
            prop = {"size": fontsize-6, "family": "Times New Roman"},
            frameon = False,
            handletextpad = 0.2,
            handlelength = 1,
            labelspacing = 0.25,
            columnspacing = 0.6,
            ncol = 2,
            facecolor = 'None',
        )
        
        Xlabel += self.genMultiName(Xpow)
        plt.subplots_adjust(top = 0.97)
        plt.subplots_adjust(left = 0.2)
        plt.subplots_adjust(right = 0.975)
        plt.subplots_adjust(bottom = 0.18)
        
        plt.yticks(np.arange(Ystart, Yend+Yinterval, step=Yinterval), fontsize=Yticks_fontsize)
        plt.xticks(x, fontsize=Xticks_fontsize)
        plt.ylabel(Ylabel, fontsize=Ylabel_fontsize)
        plt.xlabel(Xlabel, fontsize=Xlabel_fontsize, labelpad=500)
        ax1.yaxis.set_label_coords(-0.13, 0.5)
        ax1.xaxis.set_label_coords(0.45, -0.13)
        # plt.show()
    
        pdfName =  _Xlabel + "_" + _Ylabel
        plt.savefig('../data/pdf/Fie2_allPath_{}.eps'.format(pdfName)) 
        plt.savefig('../data/pdf/allPath_{}.jpg'.format(pdfName)) 
        # Xlabel = Xlabel.replace(' (%)','')
        # Xlabel = Xlabel.replace('# ','')
        # Ylabel = Ylabel.replace('# ','')
        plt.close()

    def genMultiName(self, multiple):
        if multiple == 0:
            return str()
        else:
            return "($" + "10" + "^{" + str(multiple) + "}" + "$)"

if __name__ == "__main__":
    # data檔名 Y軸名稱 X軸名稱 Y軸要除多少(10的多少次方) Y軸起始座標 Y軸終止座標 Y軸座標間的間隔
    # ChartGenerator("numOfnodes_waitingTime.txt", "need #round", "#Request of a round", 0, 0, 25, 5)
    Xlabels = ["request_cnt", "time_limit", "tao", "fidelity_threshold", "avg_memory", "min_fidelity"]
    Ylabels = ["fidelity_gain", "succ_request_cnt"]
    PathNames = ["Greedy", "QCAST", "REPS"]

    LabelsName = {}
    # LabelsName["num_nodes"] = "#Nodes"
    LabelsName["request_cnt"] = "$\#$Request"
    LabelsName["time_limit"] = "$|T|$"
    LabelsName["avg_memory"] = "Average Memory"
    LabelsName["tao"] = "$\\it{\\tau}$"
    LabelsName["fidelity_gain"] = "Expected Fidelity Sum"
    LabelsName["succ_request_cnt"] = "$\#$Accepted Request"
    LabelsName["fidelity_threshold"] = "Fidelity Threshold"
    LabelsName["min_fidelity"] = "Min Initial Fidelity"


    ChartGenerator(PathNames, "request_cnt", "fidelity_gain", LabelsName["request_cnt"], LabelsName["fidelity_gain"], 0, 0, 0, 90, 15)
    ChartGenerator(PathNames, "request_cnt", "succ_request_cnt", LabelsName["request_cnt"], LabelsName["succ_request_cnt"], 0, 0, 0, 105, 15)

    # for Path in PathNames:
    #     for Xlabel in Xlabels:
    #         for Ylabel in Ylabels:
    #             dataFileName = Path + '_' + Xlabel + '_' + Ylabel + '.ans'
    #             if Xlabel == "request_cnt":
    #                 (Ystart, Yend, Yinternal) = (5, 35, 5)
    #             # if Xlabel == "time_limit":
    #             #     (Ystart, Yend, Yinternal) = (2, 30, 7)
    #             # if Xlabel == "tao":
    #             #     (Ystart, Yend, Yinternal) = (6, 26, 5)
    #             # if Xlabel == "fidelity_threshold":
    #             #     (Ystart, Yend, Yinternal) = (15, 30, 3)
    #             # if Xlabel == "avg_memory":
    #             #     (Ystart, Yend, Yinternal) = (8, 33, 5)
    #             # if Xlabel == "min_fidelity":
    #             #     (Ystart, Yend, Yinternal) = (8, 33, 5)

    #             ChartGenerator(Path, dataFileName, LabelsName[Xlabel], LabelsName[Ylabel], 0, 0, Ystart, Yend, Yinternal)