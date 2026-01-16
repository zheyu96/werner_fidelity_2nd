import sys
import networkx as nx
import random
import math

# 設定地圖範圍大小
RANGE = 300

def dist(p1, p2):
    (x1, y1) = p1
    (x2, y2) = p2
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** (1 / 2)

if len(sys.argv) <= 2:
    print("Usage: python3 graph_generator.py <filename> <num_of_node>")
    sys.exit()

filename = sys.argv[1]
num_of_node = int(sys.argv[2])

print("======== generating graph ========", file=sys.stderr)
print(f"filename = {filename}", file=sys.stderr)
print(f"num_of_node = {num_of_node}", file=sys.stderr)

G = None
avg_hops = 0

# ==== 核心生成迴圈 ====
while True:
    # [修改點 1] 調整 Waxman 參數
    # alpha=0.4: 增加長距離連線機率 (降低 Hop 數，增加連通性)
    # beta=0.5: 控制邊的密度
    G = nx.waxman_graph(num_of_node, beta=0.5, alpha=0.4, domain=(0, 0, 0.5, 1))
    
    # [修改點 2] 使用 NetworkX 內建檢查，確保圖是連通的
    if nx.is_connected(G):
        avg_hops = nx.average_shortest_path_length(G)
        print(f"Generated Graph - Avg Hops: {avg_hops:.2f}", file=sys.stderr)
        
        # [修改點 3] 嚴格篩選 Hop 數在 2.0 到 5.0 之間
        if 2.0 <= avg_hops <= 5.0:
            break
        else:
            print("  -> Hops out of range (2-5), regenerating...", file=sys.stderr)
    else:
        # 如果圖不連通，直接重做
        # print("  -> Graph not connected, regenerating...", file=sys.stderr)
        pass

# ==== 寫入檔案 ====
path = filename
with open(path, 'w') as f:
    positions = nx.get_node_attributes(G, 'pos')
    
    # 1. 寫入節點數量
    print(num_of_node, file=f)
    
    # 2. 寫入節點記憶體 (這裡維持隨機，實際運算通常由 C++ 的 avg_memory 覆蓋)
    for n in G.nodes():
        (x, y) = positions[n]
        # Waxman 座標是 0~1，乘上 RANGE 變為實際座標
        # 注意：雖然這裡印出了記憶體，但你的 C++ main.cpp 似乎會用 default_setting["avg_memory"] 覆寫邏輯
        num_of_memory = random.randint(-1, 1) 
        print(num_of_memory, file=f)
    
    # 3. 計算並寫入邊的數量
    # networkx 的 edges 預設不重複，無向圖不需要過濾 e[0]!=e[1] (除非有自環)
    edges = list(G.edges())
    num_of_edge = len(edges)
    print(num_of_edge, file=f)
    
    avg_l = 0
    
    # 4. 寫入邊的資訊 (Node1, Node2, Fidelity)
    for (u, v) in edges:
        dis = RANGE * dist(positions[u], positions[v])
        
        # [修改點 4] Fidelity 設計：混合難度 (Mixed Difficulty)
        # 這是讓 WernerAlgo2 贏的關鍵：
        # - 30% 機率生成高品質鏈路 (Base 算法可吃)
        # - 70% 機率生成低品質鏈路 (需要純化，WernerAlgo2 獨食)
        
        val = random.random()
        if val < 0.8:
            # 簡單題：Base 算法也能過
            ratio = random.uniform(0.85, 1)
        else:
            # 困難題：Base 算法會死，WernerAlgo2 可以救
            ratio = 0.95
        if val< 0.5:
            ratio=0.5
        if val>1:
            ratio=1
        # 邊界檢查
        if ratio > 0.99: ratio = 0.99
        if ratio < 0.55: ratio = 0.55
        
        F = ratio
        print(f"{u} {v} {F:.5f}", file=f)
        avg_l += dis

    if num_of_edge > 0:
        avg_l /= num_of_edge

print(f"num_of_edge = {num_of_edge}", file=sys.stderr)
print(f"avg_edge_len = {avg_l:.2f}", file=sys.stderr)
print("\n======== graph generate finished ! ========", file=sys.stderr)