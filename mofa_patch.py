"""
临时补丁：修复 MOFA+ 中的 scipy 兼容性问题
动态代理所有 numpy 函数到 scipy
"""
import scipy
import numpy as np

# 方法1: 将常用的 numpy 函数直接映射到 scipy
numpy_functions = [
    'where', 'unique', 'nan', 'inf', 'concatenate', 'array', 'asarray',
    'zeros', 'ones', 'eye', 'dot', 'transpose', 'sum', 'mean',
    'std', 'var', 'min', 'max', 'argmin', 'argmax', 'sort',
    'argsort', 'reshape', 'squeeze', 'expand_dims', 'clip',
    'abs', 'sqrt', 'exp', 'log', 'log2', 'log10', 'power',
    'isnan', 'isinf', 'isfinite', 'any', 'all', 'arange',
    'linspace', 'logspace', 'meshgrid', 'tile', 'repeat',
    'vstack', 'hstack', 'dstack', 'stack', 'split', 'vsplit', 'hsplit',
    'newaxis', 'ndarray', 'float64', 'int64', 'bool_', 'object_',
    'empty', 'full', 'zeros_like', 'ones_like', 'empty_like', 'full_like',
    'copy', 'frombuffer', 'fromfile', 'fromfunction', 'fromiter', 'fromstring',
    'prod', 'cumprod', 'cumsum', 'diff', 'ediff1d', 'gradient', 'cross',
    'inner', 'outer', 'kron', 'trace', 'diagonal', 'diag', 'diagflat',
    'tril', 'triu', 'flip', 'fliplr', 'flipud', 'roll', 'rot90',
    'ceil', 'floor', 'trunc', 'round', 'rint', 'fix', 'sign',
    'sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'arctan2',
    'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh',
    'deg2rad', 'rad2deg', 'degrees', 'radians', 'unwrap',
    'real', 'imag', 'conj', 'conjugate', 'angle', 'real_if_close',
    'percentile', 'quantile', 'nansum', 'nanmean', 'nanstd', 'nanvar',
    'nanmin', 'nanmax', 'nanargmin', 'nanargmax', 'nanmedian', 'nanquantile',
    'median', 'average', 'bincount', 'digitize', 'histogram', 'histogram2d',
    'corrcoef', 'correlate', 'cov', 'convolve'
]

patched_count = 0
for func_name in numpy_functions:
    if hasattr(np, func_name) and not hasattr(scipy, func_name):
        setattr(scipy, func_name, getattr(np, func_name))
        patched_count += 1

# 方法2: 创建动态 __getattr__ 来捕获所有未来可能缺失的函数
original_getattr = scipy.__getattr__ if hasattr(scipy, '__getattr__') else None

def dynamic_getattr(name):
    """动态从 numpy 获取缺失的 scipy 属性"""
    # 首先尝试原始的 getattr
    if original_getattr:
        try:
            return original_getattr(name)
        except AttributeError:
            pass

    # 如果在 numpy 中存在，使用它
    if hasattr(np, name):
        attr = getattr(np, name)
        # 缓存它以避免重复查找
        setattr(scipy, name, attr)
        return attr

    # 否则抛出标准错误
    raise AttributeError(f"module 'scipy' has no attribute '{name}'")

# 替换 scipy 的 __getattr__
scipy.__getattr__ = dynamic_getattr

print(f"✓ MOFA+ 兼容性补丁已加载")
print(f"  - 静态映射: {patched_count} 个函数")
print(f"  - 动态代理: 已启用 (自动处理未来的缺失函数)")



