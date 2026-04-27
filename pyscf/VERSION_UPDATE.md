# PySCF 版本更新信息

## 最新版本
- **版本**: 2.12.1
- **发布日期**: 2026-01-27
- **官方仓库**: https://github.com/pyscf/pyscf
- **文档**: http://www.pyscf.org

## 安装更新

```bash
# 最新稳定版
pip install --upgrade pyscf

# 开发版新功能
pip install pyscf-forge

# 所有扩展
pip install pyscf[all]
```

## 新增功能（v2.12.x）

### 1. GW/RPA 方法
```python
from pyscf import gw, rpa

# G0W0 近似
gw_calc = gw.G0W0(mf, freq_int='ac')
gw_result = gw_calc.kernel()

# 直接 RPA
rpa_calc = rpa.dRPA(mf)
rpa_result = rpa_calc.kernel()
```

### 2. MC-PDFT（多组态对密度泛函理论）
```python
from pyscf import mcpdft

# MC-PDFT 计算
cas = mcscf.CASSCF(mf, 6, 8)
pdft = mcpdft.MCPDFT(cas, 'tpbe')
energy = pdft.kernel()
```

### 3. ADC（代数图解构造）
```python
from pyscf import adc

# ADC(2) 激发态
adc_calc = adc.ADC(mf)
adc_calc.nstates = 5
adc_calc.kernel()
```

### 4. 性能改进

1. **密度拟合优化**：大体系内存占用减少 30%
2. **MPI 并行**：支持 >30000 基函数的计算
3. **积分缓存**：改进的磁盘/内存管理
4. **JAX 集成**：更好的自动微分支持

## 兼容性

- Python 3.8-3.12
- NumPy 1.20+
- SciPy 1.6+
- JAX 0.4+（可选）
