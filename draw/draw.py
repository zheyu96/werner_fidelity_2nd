import matplotlib.pyplot as plt
import numpy as np


# 讀檔
with open("graph.txt", "r") as f:
    lines = f.read().strip().splitlines()

n = int(lines[0])
points = [tuple(map(float, lines[i + 1].split())) for i in range(n)]
m = int(lines[n + 1])
edges = [tuple(map(int, lines[n + 2 + i].split())) for i in range(m)]

# 顏色設定
node_color = "#800080"  # 紫色
edge_color = "#0000FF"  # 藍色

# 畫圖
plt.figure(figsize=(10, 10))
plt.axis('off')

# 畫邊
for u, v in edges:
    x1, y1 = points[u]
    x2, y2 = points[v]
    plt.plot([x1, x2], [y1, y2], color=edge_color, linewidth=1.2, alpha=0.7)

# 畫點
x_coords, y_coords = zip(*points)
# 畫點的陰影（光暈感）
plt.scatter(x_coords, y_coords, 
            s=360,  # 陰影尺寸較大
            color=node_color, 
            alpha=0.1, 
            linewidths=0)

# 畫點本體（有邊框）
plt.scatter(x_coords, y_coords, 
            s=90,                  # 節點大小
            color=node_color,      # 節點填色
            edgecolors='white',    # 白色邊框
            linewidths=0.8, 
            zorder=2)

# 儲存圖片
plt.tight_layout()
plt.savefig("topo.png", dpi=300, bbox_inches='tight')
plt.savefig("topo.eps", dpi=300, bbox_inches='tight')

# 顯示圖形
# plt.show()
