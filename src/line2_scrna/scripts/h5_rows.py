"""h5_rows：h5ad 稀疏矩阵的行块流式读取器（Line 2 / Line 3 公共核心库）。

为什么需要这个库（2026-11-05 真机实测，不是推测）：

anndata 在 backed 模式下做子集物化（`adata[rows, cols].to_memory()`）时，会先把**整个**
`layers` 里的稀疏矩阵读进内存，再做切片。本基准数据 nnz = 2.1 亿：

  - 仅 indices 一个数组就是 802 MiB（numpy._core._exceptions._ArrayMemoryError 实测）
  - 叠加 data(float32) 与 indptr，单层约 2.4 GB
  - 两份 h5ad 同时加载必然超过手机 11.5 GB 的可用量

后果不是慢慢变卡，而是 Android 内核低内存杀手直接杀掉 Termux 进程：

  lowmemorykiller: Kill "com.termux" (26208), oom_score_adj 0
    to free 63172kB rss, 44028kB swap; reason: device is in direct reclaim

4 秒粒度采样显示 Python RSS 在 12 秒内从 476MB 涨到 9,063MB。换成下面的按行区间读取后，
同一任务峰值 792MB。

使用纪律（静态门禁会检查）：
  在两条线（本库服务的全部消费方）的脚本里，一旦以 backed="r" 打开 h5ad，
  就不允许对该对象调用 .to_memory() 或对它做全量 .X[...] 切片；必须走本模块的
  read_csr_rows / read_var_column。
"""

import numpy as np
import scipy.sparse as sps
import h5py


def read_var_column(h5_path, col, group="var"):
    """读 h5ad 的 var/obs 一列。

    字符串列在 h5ad 里可能是普通 dataset，也可能是 categorical 组（categories + codes），
    两种形态都要支持，否则会在「看起来能读」的地方报 TypeError。
    """
    with h5py.File(h5_path, "r") as f:
        node = f[group][col]
        if isinstance(node, h5py.Dataset):
            raw = node[...]
            return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in raw])
        cats = node["categories"][...]
        codes = node["codes"][...]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        return cats[codes]


def read_obs_index(h5_path):
    """读 obs 的行索引。索引列名存在 obs 组的 _index 属性里，不一定叫 _index。"""
    with h5py.File(h5_path, "r") as f:
        name = f["obs"].attrs.get("_index", "_index")
        if isinstance(name, bytes):
            name = name.decode()
        node = f["obs"][name]
        if isinstance(node, h5py.Dataset):
            return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in node[...]])
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][...]])
        return cats[node["codes"][...]]


def matrix_shape(h5_path, key="X"):
    """读稀疏矩阵的 (n_obs, n_vars)，不加载任何元素。"""
    with h5py.File(h5_path, "r") as f:
        return tuple(int(x) for x in f[key].attrs["shape"])


def read_csr_rows(h5_path, key, rows):
    """按行切片 CSR 稀疏矩阵，内存与选中行数成正比。

    只读 indptr/data/indices 的对应区间；相邻行合并成连续区间以减少 h5py 调用次数。
    非 CSR 编码直接抛错，不静默走回退路径（CSC 上按行切片必然错）。
    """
    with h5py.File(h5_path, "r") as f:
        g = f[key]
        enc = str(g.attrs.get("encoding-type", ""))
        if "csr" not in enc.lower():
            raise RuntimeError("read_csr_rows 只支持 CSR，实际 %s（key=%s）" % (enc, key))
        shape = tuple(int(x) for x in g.attrs["shape"])
        indptr = np.asarray(g["indptr"][...], dtype=np.int64)
        rows = np.asarray(rows, dtype=np.int64)
        if rows.size == 0:
            return sps.csr_matrix((0, shape[1]), dtype=np.float32)
        spans = []
        start = prev = int(rows[0])
        for r in rows[1:]:
            r = int(r)
            if r == prev + 1:
                prev = r
                continue
            spans.append((start, prev))
            start = prev = r
        spans.append((start, prev))
        out_d, out_i = [], []
        for lo_r, hi_r in spans:
            lo, hi = int(indptr[lo_r]), int(indptr[hi_r + 1])
            if hi > lo:
                out_d.append(np.asarray(g["data"][lo:hi]))
                out_i.append(np.asarray(g["indices"][lo:hi]))
        # indptr 必须跨区间累计：每段偏移相对本段起点，直接拼会得到负数
        out_ptr = [0]
        running = 0
        for r in rows:
            running += int(indptr[r + 1] - indptr[r])
            out_ptr.append(running)
    data = np.concatenate(out_d) if out_d else np.zeros(0, dtype=np.float32)
    indices = np.concatenate(out_i) if out_i else np.zeros(0, dtype=np.int32)
    return sps.csr_matrix((data, indices, np.asarray(out_ptr, dtype=np.int64)),
                          shape=(rows.size, shape[1]))


def iter_row_blocks(n_obs, block_rows=2000):
    """按块产出 (start, stop)，供流式累加使用。"""
    if block_rows < 1:
        raise ValueError("block_rows 必须 >= 1，实际 %r" % (block_rows,))
    for start in range(0, int(n_obs), int(block_rows)):
        yield start, min(start + int(block_rows), int(n_obs))
