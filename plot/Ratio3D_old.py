import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm  # 引入 colormap 模組

def plot_3d_bar_chart(x_values, y_layers, z_values_list):
    fig = plt.figure(figsize=(8.5, 7))
    ax = fig.add_subplot(projection='3d')
    yticks = list(range(len(y_layers)))
    xticks = list(range(len(x_values)))
    for k, z_values in zip(yticks, z_values_list):
        if len(x_values) != len(z_values):
            print(len(x_values), len(z_values))
            raise ValueError("每層的 z_values 必須和 x_values 的長度一致。")
        
        norm = plt.Normalize(1, 2)  # 標準化 z 值範圍
        colors = cm.viridis_r(norm(z_values))
        ax.bar(np.arange(len(x_values)), z_values, zs=k, zdir='y', color=colors, alpha=1, width=0.75)

    # 設定軸標籤和刻度
    ax.set_xlabel('$\\tau (10^{-3})$', fontsize=16, labelpad=16)
    ax.set_ylabel('$T$', fontsize=16, labelpad=20)
    ax.set_zlabel(r'$\frac{F_{\text{max}}}{F_{\text{min}}}$', rotation=0, labelpad=20, fontsize=16)
    
    # 設定標籤位置調整
    # ax.xaxis.label.set_position((0.5, -0.1))  # (水平位置, 垂直位置)

    ax.set_yticks(yticks)
    ax.set_xticks(xticks)
    ax.set_xticklabels(x_values)
    ax.set_yticklabels(y_layers)

    plt.subplots_adjust(bottom=0, top=1, left=0, right=1.1)

    ax.view_init(elev=30, azim=225)  # 改變視角來拉長 Z 軸
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    # plt.show()
    plt.savefig('../data/pdf/{}.eps'.format("FmaxFminRatio")) 
    plt.savefig('../data/pdf/{}.png'.format("FmaxFminRatio")) 

x_len = 8
y_len = 15
x_values = np.arange(1e-3, 4e-3, 2e-4)
x_values = x_values[0 : x_len]
x_values = [round(x * (10 ** 3), 2) for x in x_values]
y_layers = np.arange(1, y_len, 1)
z_values_list = []

for i in range(y_len):
    one_z_list = []
    for j in range(x_len):
        x, y, z = map(float, input().split(' '))
        one_z_list.append(z)
    z_values_list.append(one_z_list)

plot_3d_bar_chart(x_values, y_layers, z_values_list)
