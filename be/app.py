from flask import Flask, request, jsonify
import scanpy as sc
import scvelo as scv
import numpy as np
import json
from flask_cors import CORS
from sklearn.neighbors import kneighbors_graph
# import arviz as az

app = Flask(__name__)
CORS(app)

# class NumpyEncoder(json.JSONEncoder):
#     def default(self, obj):
#         if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
#             np.int16, np.int32, np.int64, np.uint8,
#             np.uint16,np.uint32, np.uint64)):
#             return int(obj)
#         elif isinstance(obj, (np.float_, np.float16, np.float32,
#             np.float64)):
#             return float(obj)
#         elif isinstance(obj, (np.bool_)):
#             return bool(obj)
#         elif isinstance(obj, (np.ndarray,)):
#             return obj.tolist()
#         return json.JSONEncoder.default(self, obj)
# app.json_encoder = NumpyEncoder

def jsonify_safe(obj):
    import numpy as np

    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: jsonify_safe(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [jsonify_safe(v) for v in obj]
    elif isinstance(obj, np.generic):
        return obj.item()
    else:
        return obj

# grid
def grid_partition(X, grid_num=20, min_thres=5):

    xmin, ymin = X.min(axis=0)
    xmax, ymax = X.max(axis=0)

    x_range = xmax - xmin
    y_range = ymax - ymin

    cell_size = max(x_range, y_range) / grid_num

    x_bins = xmin + np.arange(grid_num + 1) * cell_size
    y_bins = ymin + np.arange(grid_num + 1) * cell_size

    x_idx = np.digitize(X[:, 0], x_bins) - 1
    y_idx = np.digitize(X[:, 1], y_bins) - 1

    x_idx = np.clip(x_idx, 0, grid_num - 1)
    y_idx = np.clip(y_idx, 0, grid_num - 1)

    flat_id = x_idx * grid_num + y_idx

    unique_ids, counts = np.unique(flat_id, return_counts=True)
    valid_ids = unique_ids[counts >= min_thres]

    id_map = {old: new for new, old in enumerate(valid_ids)}

    point_labels = np.full(len(X), -1, dtype=int)
    for i, fid in enumerate(flat_id):
        if fid in id_map:
            point_labels[i] = id_map[fid]

    grid_centers = []
    grid_sizes = []

    for fid in valid_ids:
        xi = fid // grid_num
        yi = fid % grid_num

        cx = 0.5 * (x_bins[xi] + x_bins[xi + 1])
        cy = 0.5 * (y_bins[yi] + y_bins[yi + 1])

        grid_centers.append([cx, cy])
        grid_sizes.append(np.sum(flat_id == fid))


    grid_xmin = x_bins[0]
    grid_xmax = x_bins[-1]
    grid_ymin = y_bins[0]
    grid_ymax = y_bins[-1]

    grid_extent = [[grid_xmin, grid_xmax], 
                    [grid_ymin, grid_ymax]]

    return (
        point_labels,
        np.asarray(grid_centers),
        cell_size,
        grid_extent
    )
    
    
# def getHDI(data, hdi_thres=0.5, filter_relS_thres=0.3):
    
#     import arviz as az
    
#     intervals = az.hdi(
#         data,
#         hdi_prob=hdi_thres,
#         multimodal=True
#     )
#     intervals_counts = np.array([
#         np.sum((data >= l) & (data <= r))
#         for l, r in intervals
#     ])
#     intervals_strength = intervals_counts / intervals_counts.sum()
    
#     strength_mask = intervals_strength >= filter_relS_thres
    
#     if not np.any(strength_mask):
#         # 保留 top-k 区间，而不是全删
        
#         k = min(4, len(intervals))
#         top_idx = np.argsort(intervals_strength)[-k:]
#         strength_mask = np.zeros_like(intervals_strength, dtype=bool)
#         strength_mask[top_idx] = True
        
    
#     # 筛选掉强度过低的区间
#     intervals = intervals[strength_mask]
#     intervals_strength = intervals_strength[strength_mask]

#     return intervals,intervals_strength

def getHDI2(data, hdi_thres=0.5, filter_relS_thres=0.3, bins=12):
    data = np.asarray(data)
    n = len(data)

    if n == 0:
        return np.empty((0, 2)), np.empty((0,))

    # 1️⃣ 构建 bins
    data_min, data_max = data.min(), data.max()
    if data_min == data_max:
        # 所有值相同，退化为单点区间
        intervals = np.array([[data_min, data_max]])
        return intervals, np.array([1.0])

    bin_edges = np.linspace(data_min, data_max, bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # 2️⃣ 统计每个 bin 的计数
    counts, _ = np.histogram(data, bins=bin_edges)
    probs = counts / counts.sum()  # 概率质量

    # 3️⃣ 按概率从大到小排序 bin
    order = np.argsort(probs)[::-1]

    selected = []
    cum_prob = 0.0
    for idx in order:
        if probs[idx] == 0:
            break
        selected.append(idx)
        cum_prob += probs[idx]
        if cum_prob >= hdi_thres:
            break

    selected = np.array(sorted(selected))

    # 4️⃣ 把相邻 bin 合并成区间
    intervals = np.array([
        [bin_edges[i], bin_edges[i + 1]]
        for i in selected
    ])
    strengths = probs[selected]

    # 5️⃣ 归一化区间强度（和你原函数一致）
    strengths = strengths / strengths.sum()

    # 6️⃣ 过滤弱区间（逻辑保持一致）
    mask = strengths >= filter_relS_thres

    # if not np.any(mask):
    #     k = min(4, len(strengths))
    #     top_idx = np.argsort(strengths)[-k:]
    #     mask = np.zeros_like(strengths, dtype=bool)
    #     mask[top_idx] = True


    # if not np.any(mask):
        # k = min(4, len(strengths))
    k = 2
    top_idx = np.argsort(strengths)[-k:]
    mask = np.zeros_like(strengths, dtype=bool)
    mask[top_idx] = True

    return intervals[mask], strengths[mask]

def dir_partition(data,bins=12):
    """
    theta: shape (n,), values in (-π, π]
    bins:  number of bins
    
    return:
        bin_ranges: shape (bins, 2)
        intensity:  shape (bins,)
    """
    
    # 1️⃣ 定义边界
    edges = np.linspace(-np.pi, np.pi, bins + 1)
    
    # 2️⃣ 统计频数
    counts, _ = np.histogram(data, bins=edges)
    
    # 3️⃣ 归一化为比例
    intensity = counts / len(data)
    
    # 4️⃣ 返回每个bin范围
    bin_ranges = np.stack([edges[:-1], edges[1:]], axis=1)
    
    return bin_ranges, intensity

def accumulate_angle_density(intervals, strengths, bins=120, total_samples=3600):
    """
    intervals: list of (l, r)
    strengths: list of s
    """
    theta = np.linspace(-np.pi, np.pi, bins + 1)
    theta_centers = (theta[:-1] + theta[1:]) / 2
    density = np.zeros(bins)

    for (l, r), s in zip(intervals, strengths):
        mask = (theta[:-1] >= l) & (theta[1:] <= r)
        if np.any(mask):
            density[mask] += s / (r - l)

    density /= density.sum()
    counts = np.round(density * total_samples).astype(int)
    if counts.sum() == 0:
        counts[np.argmax(density)] = total_samples
    samples = np.repeat(theta_centers, counts)

    return samples

def read_data(project_name):
    adata = sc.read_h5ad(f"./data/{project_name}.h5ad")
    return adata

# # 返回按基因的grid集成可视化视图数据
# @app.route("/init_plot_GridVis1", methods=["POST"])
# def init_plot_GridVis1():
    
    # 读取参数
    reqParams = json.loads(request.get_data())
    project_name = reqParams['project_name']
    adata = read_data(project_name)
    
    embeddings = adata.obsm['X_embedding']
    velo2D = adata.uns['velo2D']
    clusters_color = adata.obs['clusters_color']
    velocity_embedding = adata.obsm['velocity_embedding']
    

    # grid partition
    point_gridLabels, grid_centers, grid_size, gird_bounds = grid_partition(embeddings,grid_num=30)

    # aggregate in grid
    VeloDirCI = []
    meanVelos = []
    grid_colors = []

    valid_grids = np.unique(point_gridLabels)
    valid_grids = valid_grids[valid_grids >= 0]
    for gid in valid_grids:
        # 提取gird中的速度向量
        mask = point_gridLabels == gid
        data = velo2D[mask]  # n_points * n_genes * 2

        data = data.reshape(-1, 2)
        angles = np.arctan2(data[:, 1], data[:, 0])

        # data = data.reshape(data.shape[0],-1, 2)

        # grid_intervals = []
        # grid_intervals_strengths = []

        # for cell_velos in data: # 以细胞为单位先做HDI
        #     cell_angles =  np.arctan2(cell_velos[:, 1], cell_velos[:, 0])
        #     intervals,intervals_strength = getHDI(cell_angles,hdi_thres=0.7, filter_relS_thres=0.2)
        #     for i, s in zip(intervals,intervals_strength):
        #         grid_intervals.append(i)
        #         grid_intervals_strengths.append(s)


        # grid_samples = accumulate_angle_density(
        #     grid_intervals,
        #     grid_intervals_strengths,
        #     bins=360,
        #     total_samples=3600,
        # )

        # mean_intervals, mean_intervals_strengths = getHDI(grid_samples, hdi_thres=0.7, filter_relS_thres=0.2)    
        mean_intervals, mean_intervals_strengths = getHDI2(angles, hdi_thres=1, filter_relS_thres=0,bins=12)    
        # 计算平均值
        # meanVelo = np.mean(velocity_embedding[mask],axis=0).tolist()
        meanVelo = np.mean(data.reshape(-1,2),axis=0)
        meanVelos.append(meanVelo)
        
        # 计算grid的聚类
        grid_colors.append(
            np.unique(clusters_color[mask], return_counts=True)[0][np.unique(clusters_color[mask], return_counts=True)[1].argmax()]        
        )
            
        VeloDirCI.append([{'interval':interval,'strength':strength} for interval,strength in zip(mean_intervals,mean_intervals_strengths)])



    data = {}
    data['name'] = project_name
    data['type'] = 'GridVis1'
    data['nCells'] = adata.shape[0]
    data['nGenes'] = adata.shape[1]
    data['embeddings'] = adata.obsm['X_embedding']
    data['grid_assign'] = point_gridLabels.tolist()
    data['grid_pos'] = grid_centers.tolist()
    data['grid_size'] = grid_size
    data['grid_bounds'] = gird_bounds
    # data['GeneVelos2D'] = velo2D.tolist()
    data['VeloDirCI'] = VeloDirCI
    data['meanVelo'] = meanVelos
    data['grid_colors'] = grid_colors
    data['Genes'] = adata.var_names.tolist()

    return jsonify(jsonify_safe(data))

# # 返回按基因模块的grid集成可视化视图数据
# @app.route("/init_plot_GridVis1", methods=["POST"])
# def init_plot_GridVis1():
    
    # 读取参数
    reqParams = json.loads(request.get_data())
    project_name = reqParams['project_name']
    adata = read_data(project_name)
    
    embeddings = adata.obsm['X_embedding']
    velo2D = adata.uns['velo2D']['all']
    clusters_color = adata.obs['clusters_color']
    
    # grid partition
    point_gridLabels, grid_centers, grid_size, gird_bounds = grid_partition(embeddings,grid_num=30)

    # aggregate in grid
    VeloDirCI = []
    meanVelos = []
    grid_colors = []

    valid_grids = np.unique(point_gridLabels)
    valid_grids = valid_grids[valid_grids >= 0]
    for gid in valid_grids:
        # 提取gird中的速度向量
        mask = point_gridLabels == gid
        data = velo2D[mask]  # n_points * n_genes * 2

        data = data.reshape(-1, 2)
        angles = np.arctan2(data[:, 1], data[:, 0])

        # data = data.reshape(data.shape[0],-1, 2)

        # grid_intervals = []
        # grid_intervals_strengths = []

        # for cell_velos in data: # 以细胞为单位先做HDI
        #     cell_angles =  np.arctan2(cell_velos[:, 1], cell_velos[:, 0])
        #     intervals,intervals_strength = getHDI(cell_angles,hdi_thres=0.7, filter_relS_thres=0.2)
        #     for i, s in zip(intervals,intervals_strength):
        #         grid_intervals.append(i)
        #         grid_intervals_strengths.append(s)


        # grid_samples = accumulate_angle_density(
        #     grid_intervals,
        #     grid_intervals_strengths,
        #     bins=360,
        #     total_samples=3600,
        # )

        # mean_intervals, mean_intervals_strengths = getHDI(grid_samples, hdi_thres=0.7, filter_relS_thres=0.2)    
        # mean_intervals, mean_intervals_strengths = getHDI2(angles, hdi_thres=1, filter_relS_thres=0,bins=12)
        
        mean_intervals, mean_intervals_strengths = dir_partition(angles, bins=12)    

        # 计算平均值
        # meanVelo = np.mean(velocity_embedding[mask],axis=0).tolist()
        meanVelo = np.mean(data.reshape(-1,2),axis=0)
        meanVelos.append(meanVelo)
        
        # 计算grid的聚类
        grid_colors.append(
            np.unique(clusters_color[mask], return_counts=True)[0][np.unique(clusters_color[mask], return_counts=True)[1].argmax()]        
        )
            
        VeloDirCI.append([{'interval':interval,'strength':strength} for interval,strength in zip(mean_intervals,mean_intervals_strengths)])



    data = {}
    data['name'] = project_name
    data['type'] = 'GridVis1'
    data['nCells'] = adata.shape[0]
    data['nGenes'] = adata.shape[1]
    data['embeddings'] = adata.obsm['X_embedding']
    data['grid_assign'] = point_gridLabels.tolist()
    data['grid_pos'] = grid_centers.tolist()
    data['grid_size'] = grid_size
    data['grid_bounds'] = gird_bounds
    # data['GeneVelos2D'] = velo2D.tolist()
    data['VeloDirCI'] = VeloDirCI
    data['meanVelo'] = meanVelos
    data['grid_colors'] = grid_colors
    data['Genes'] = adata.var_names.tolist()

    return jsonify(jsonify_safe(data))
    

# # 返回按基因z组后的模块的grid集成可视化视图数据
# @app.route("/init_plot_GridVis1", methods=["POST"])
# def init_plot_GridVis1():
    
#     # 读取参数
#     reqParams = json.loads(request.get_data())
#     project_name = reqParams['project_name']
#     adata = read_data(project_name)
    
#     embeddings = adata.obsm['X_embedding']
#     velo2D = adata.uns['velo2D']
#     clusters_color = adata.obs['clusters_color']
    
#     # grid partition
#     point_gridLabels, grid_centers, grid_size, gird_bounds = grid_partition(embeddings,grid_num=30)

#     # aggregate in grid
#     VeloDirCI = []
#     meanVelos = []
#     grid_colors = []

#     valid_grids = np.unique(point_gridLabels)
#     valid_grids = valid_grids[valid_grids >= 0]
        
    
#     for gid in valid_grids:
#         # 提取gird中的速度向量
#         mask = point_gridLabels == gid

#         angles = []
#         datas = []

#         for m in velo2D.keys():
#             data = velo2D[m][mask]  # n_points * n_genes * 2
#             data = data.reshape(-1, 2)
#             angles.append(np.arctan2(data[:, 1], data[:, 0]))
#             datas.append(data)
#         angles = np.concatenate(angles, axis=0)
#         datas = np.concatenate(datas, axis=0)

#         mean_intervals, mean_intervals_strengths = dir_partition(angles, bins=12)    

#         # 计算平均值
#         # meanVelo = np.mean(velocity_embedding[mask],axis=0).tolist()
#         meanVelo = np.mean(data.reshape(-1,2),axis=0)
#         meanVelos.append(meanVelo)
        
#         # 计算grid的聚类
#         grid_colors.append(
#             np.unique(clusters_color[mask], return_counts=True)[0][np.unique(clusters_color[mask], return_counts=True)[1].argmax()]        
#         )
            
#         VeloDirCI.append([{'interval':interval,'strength':strength} for interval,strength in zip(mean_intervals,mean_intervals_strengths)])



#     data = {}
#     data['name'] = project_name
#     data['type'] = 'GridVis1'
#     data['nCells'] = adata.shape[0]
#     data['nGenes'] = adata.shape[1]
#     data['embeddings'] = adata.obsm['X_embedding']
#     data['grid_assign'] = point_gridLabels.tolist()
#     data['grid_pos'] = grid_centers.tolist()
#     data['grid_size'] = grid_size
#     data['grid_bounds'] = gird_bounds
#     # data['GeneVelos2D'] = velo2D.tolist()
#     data['VeloDirCI'] = VeloDirCI
#     data['meanVelo'] = meanVelos
#     data['grid_colors'] = grid_colors
#     data['Genes'] = adata.var_names.tolist()

#     return jsonify(jsonify_safe(data))


# 返回不同RNA速率计算方法的grid集成可视化视图数据
@app.route("/init_plot_GridVis1", methods=["POST"])
def init_plot_GridVis1():
    
    # 读取参数
    reqParams = json.loads(request.get_data())
    method_list = reqParams['project_name']
    velo2D_list = []
    for method in method_list:
        velo2D_list.append(read_data(method).uns['velo2D'])
    embeddings = read_data(method_list[0]).obsm['X_embedding']
    clusters_color = read_data(method_list[0]).obs['clusters_color']
    nCells = read_data(method_list[0]).shape[0]
    nGenes = read_data(method_list[0]).shape[1]
    
    
    # grid partition
    point_gridLabels, grid_centers, grid_size, gird_bounds = grid_partition(embeddings,grid_num=30)

    # aggregate in grid
    VeloDirCI = []
    meanVelos = []
    grid_colors = []

    valid_grids = np.unique(point_gridLabels)
    valid_grids = valid_grids[valid_grids >= 0]
        
    
    for gid in valid_grids:
        # 提取gird中的速度向量
        mask = point_gridLabels == gid
                
        data = np.vstack([velo2D[mask] for velo2D in velo2D_list])
        
        angles = np.arctan2(data[:, 1], data[:, 0])
        
        mean_intervals, mean_intervals_strengths = dir_partition(angles, bins=12)    

        # 计算平均值
        # meanVelo = np.mean(velocity_embedding[mask],axis=0).tolist()
        meanVelo = np.mean(data,axis=0)
        meanVelos.append(meanVelo)
        
        # 计算grid的聚类
        grid_colors.append(
            np.unique(clusters_color[mask], return_counts=True)[0][np.unique(clusters_color[mask], return_counts=True)[1].argmax()]        
        )
            
        VeloDirCI.append([{'interval':interval,'strength':strength} for interval,strength in zip(mean_intervals,mean_intervals_strengths)])



    data = {}
    data['name'] = 'multi method'
    data['type'] = 'GridVis1'
    data['nCells'] = nCells
    data['nGenes'] = nGenes
    data['embeddings'] = embeddings
    data['grid_assign'] = point_gridLabels.tolist()
    data['grid_pos'] = grid_centers.tolist()
    data['grid_size'] = grid_size
    data['grid_bounds'] = gird_bounds
    # data['GeneVelos2D'] = velo2D.tolist()
    data['VeloDirCI'] = VeloDirCI
    data['meanVelo'] = meanVelos
    data['grid_colors'] = grid_colors

    return jsonify(jsonify_safe(data))


# 根据给出的gene或者基因的组合，返回新的置信区间的平均速度
@app.route("/update_plot_GridVis1", methods=["POST"])
def update_plot():
    
    # 读取参数
    reqParams = json.loads(request.get_data())
    project_name = reqParams['project_name']
    gene_list = reqParams['gene_list']
    adata = read_data(project_name)

    velo2D = adata.uns['velo2D']
    
    # 按基因过滤
    velo2D = velo2D[:,adata.var_names.get_indexer(gene_list),...]
    adata = adata[:, gene_list]
    point_gridLabels = adata.obs['point_gridLabels']

    # 更新GeneVeloDirCI和meanVelo
    GeneVeloDirCI = []
    meanVelos = []
    scv.tl.velocity(adata)
    # scv.tl.velocity_graph(adata)
    del adata.obsm['velocity_embedding']
    scv.tl.velocity_embedding(adata,basis='embedding')
    velocity_embedding = adata.obsm['velocity_embedding']
    valid_grids = np.unique(point_gridLabels)
    valid_grids = valid_grids[valid_grids >= 0]
    for gid in valid_grids:
        # 提取gird中的速度向量
        mask = point_gridLabels == gid
        data = velo2D[mask]  # n_points * n_genes * 2
        data = data.reshape(-1, 2)
        angles = np.arctan2(data[:, 1], data[:, 0])
        # grid_intervals = []
        # grid_intervals_strengths = []

        # for cell_velos in data: # 以细胞为单位先做HDI
        #     cell_angles =  np.arctan2(cell_velos[:, 1], cell_velos[:, 0])
        #     intervals,intervals_strength = getHDI(cell_angles,hdi_thres=0.7, filter_relS_thres=0.2)
        #     for i, s in zip(intervals,intervals_strength):
        #         grid_intervals.append(i)
        #         grid_intervals_strengths.append(s)


        # grid_samples = accumulate_angle_density(
        #     grid_intervals,
        #     grid_intervals_strengths,
        #     bins=360,
        #     total_samples=3600,
        # )

        mean_intervals, mean_intervals_strengths = getHDI2(angles, hdi_thres=1, filter_relS_thres=0,bins=12)    
        
        # 计算平均值
        # meanVelo = np.mean(velocity_embedding[mask],axis=0).tolist()
        meanVelo = np.mean(data.reshape(-1,2),axis=0)
        meanVelos.append(meanVelo)
        
            
        GeneVeloDirCI.append([{'interval':interval,'strength':strength} for interval,strength in zip(mean_intervals,mean_intervals_strengths)])

    
    return jsonify(jsonify_safe({
        'GeneVeloDirCI': GeneVeloDirCI,
        'meanVelo': meanVelos,
    }))


# 简单的箭头指向
@app.route("/init_plot_GridVis2", methods=["POST"])
def init_plot_GridVis2():
    # 读取参数
    reqParams = json.loads(request.get_data())
    project_name = reqParams['project_name']
    adata = read_data(project_name)

    embedding = adata.obsm['X_embedding']
    cluster_color = adata.obs['clusters_color']
    velocity_embedding = adata.uns['velo2D']['all']


    return jsonify(jsonify_safe({
        'embedding': embedding.tolist(),
        'velocity_embedding': velocity_embedding.tolist(),
        'cluster_color':cluster_color.tolist(),
        'type':'GridVis2'
    }))



if __name__ == '__main__': ##!important vscode的debug不会走该路径执行该函数..

    app.run(port = 5005,debug=True)
