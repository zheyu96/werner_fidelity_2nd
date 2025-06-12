import matplotlib.pyplot as plt

# 讀檔
with open("graph.txt", "r") as f:
    lines = f.read().strip().splitlines()

n = int(lines[0])
points = [tuple(map(float, lines[i + 1].split())) for i in range(n)]
m = int(lines[n + 1])
edges = [tuple(map(int, lines[n + 2 + i].split())) for i in range(m)]

# 顏色設定
node_color = "#800080"  # 節點主色（深紫色）
edge_color = "#0000FF"  # 邊的顏色（藍色）

# --- 畫圖 ---
plt.figure(figsize=(10, 10))
plt.axis('off')

# 畫邊
for u, v in edges:
    x1, y1 = points[u]
    x2, y2 = points[v]
    plt.plot([x1, x2], [y1, y2], color=edge_color, linewidth=1.2)

# 外層淡紫色模擬光暈（更粗！）
# 把半徑層級增大 1.3 倍
shadow_layers = [
    (460, '#f2e6f2'),  # 淺粉紫色
    (380, '#e0b0e0'),  # 淺紫色
    (300, '#c080c0'),  # 中紫色
]

for size, color in shadow_layers:
    plt.scatter(*zip(*points),
                s=size,
                color=color,
                edgecolors='none',
                zorder=1)

# 節點本體（有白色邊框）
plt.scatter(*zip(*points),
            s=90,
            color=node_color,
            edgecolors='white',
            linewidths=0.8,
            zorder=2)

plt.tight_layout()
plt.savefig("topo.eps", dpi=1600, bbox_inches='tight')
plt.savefig("topo.png", dpi=1600, bbox_inches='tight')

# plt.show()  # 可選
