from cellphonedb.src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method
import scanpy as sc
import pandas as pd
import os
from tqdm import tqdm


def find_files(directory, substring, extension):
    # 存储符合条件的文件路径
    matching_files = []
    # 遍历给定目录及其所有子目录
    for root, dirs, files in os.walk(directory):
        # 检查每个文件
        for file in files:
            # 检查文件名是否包含特定子字符串并且以特定扩展名结尾
            if substring in file and file.endswith(extension):
                # 如果符合条件，添加完整路径到列表
                matching_files.append(os.path.join(root, file))
    return matching_files


save_path = '/data/twang15/spatial_protein/CODEX_data/hubmap/cci/large_intestine'
os.makedirs(save_path, exist_ok=True)
# tissues = adata.obs['Tissue'].drop_duplicates().tolist()
files = find_files('/data/twang15/spatial_protein/CODEX_data/hubmap/integrated/large_intestine', 'integrated', '.h5ad')
files.remove('/data/twang15/spatial_protein/CODEX_data/hubmap/integrated/large_intestine/combined_protein_integrated.h5ad')
for file in tqdm(files):
    tissue = file.split('/')[8].replace('_integrated_protein.h5ad', '')
    os.makedirs(f'{save_path}/{tissue}/', exist_ok=True)
    tmp = sc.read_h5ad(file)
    tmp.var_names = tmp.var['feature_name']
    # tmp.var["mt"] = tmp.var_names.str.startswith("MT-")
    tmp = tmp[:, ~tmp.var_names.str.startswith("MT-")]
    tmp = tmp[:, ~tmp.var_names.str.startswith('ENSG')]
    tmp.obs_names_make_unique()
    tmp.var_names_make_unique()
    sc.pp.filter_genes(tmp, min_cells=10)
    # tmp.obs_names = pd.Index(range(tmp.n_obs))
    tmp.write_h5ad(f'{save_path}/{tissue}/adata.h5ad')
    meta = tmp[tmp.obs['Tissue'] == tissue].obs[['cell_type']]
    meta.to_csv(f'{save_path}/{tissue}/input_meta.txt', sep='\t', index=True, index_label='Cell')
#     counts = tmp.to_df().T
    # counts = pd.DataFrame(tmp.X.toarray(), columns=tmp.var_names, index=tmp.obs_names).T
    # counts
#     counts.to_csv(f'{save_path}/{tissue}/counts.txt', sep='\t', index=True, index_label='Gene')
    cpdb_results = cpdb_statistical_analysis_method.call(
            cpdb_file_path = "/data/twang15/spatial_protein/code/cellphonedb.zip",
            meta_file_path = f'{save_path}/{tissue}/input_meta.txt',
            counts_file_path = f'{save_path}/{tissue}/adata.h5ad',
            counts_data = 'gene_name',
            score_interactions = True,
            output_path = f'{save_path}/{tissue}/')
    del tmp