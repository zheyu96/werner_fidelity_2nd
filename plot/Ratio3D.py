import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm  # 引入 colormap 模組

def plot_3d_bar_chart(x_values, y_layers, z_values_list):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('mathtext', fontset='stix')

    fig = plt.figure(figsize=(8.5, 7))
    ax = fig.add_subplot(projection='3d')
    yticks = list(range(len(y_layers)))
    xticks = list(range(len(x_values)))
    # for k, z_values in zip(yticks, z_values_list):
    #     if len(x_values) != len(z_values):
    #         print(len(x_values), len(z_values))
    #         raise ValueError("每層的 z_values 必須和 x_values 的長度一致。")
        
    #     norm = plt.Normalize(1, 2)  # 標準化 z 值範圍
    #     colors = cm.viridis_r(norm(z_values))
    #     ax.bar(np.arange(len(x_values)), z_values, zs=k, zdir='y', color=colors, alpha=1, width=0.75)

    X, Y = np.meshgrid(xticks, yticks)
    # print(X, Y)
    # 將 z_values_list 轉換成 numpy 陣列以便繪製
    Z = np.array(z_values_list)

    # 使用 plot_surface 來繪製連續的表面
    norm = plt.Normalize(1, 2)  # 標準化 z 值範圍
    colors = cm.viridis_r(norm(Z))

    # 繪製表面
    surf = ax.plot_surface(X, Y, Z, facecolors=colors, rstride=1, cstride=1, alpha=1, linewidth=0, antialiased=False)

    my_fontsize = 22
    ax.set_xlabel('$\\tau (10^{-3})$', fontsize=my_fontsize, labelpad=16)
    ax.set_ylabel('$|T|$', fontsize=my_fontsize, labelpad=20)
    ax.set_zlabel('$\\frac{F_{\\max}}{F_{\\min}}$', rotation=0, labelpad=25, fontsize=my_fontsize+8)
    
    # 設定標籤位置調整
    # ax.xaxis.label.set_position((0.5, -0.1))  # (水平位置, 垂直位置)



    ax.set_yticks(yticks[0::2])
    ax.set_xticks(xticks[0::2])
    ax.set_xticklabels(x_values[0::2], fontsize = my_fontsize)
    ax.set_yticklabels(y_layers[0::2], fontsize = my_fontsize)
    ax.tick_params(axis='z', labelsize = my_fontsize)  # 調整 z 軸數字字體大小



    plt.subplots_adjust(bottom=0, top=1.049, left=0, right=1.1)

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
y_layers = np.arange(1, y_len + 1, 1)
# z_values_list = [[] for _ in range(x_len)]
# for i in range(y_len):
#     for j in range(x_len):
#         x, y, z = map(float, input().split(' '))
#         z_values_list[j].append(z)
z_values_list = []
for i in range(y_len):
    one_z_list = []
    for j in range(x_len):
        x, y, z = map(float, input().split(' '))
        one_z_list.append(z)
    z_values_list.append(one_z_list)

plot_3d_bar_chart(x_values, y_layers, z_values_list)
