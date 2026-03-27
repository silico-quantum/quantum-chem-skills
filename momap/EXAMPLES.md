# MOMAP 使用示例

## 示例 1: 苯分子荧光光谱计算

### 步骤 1: 准备 Gaussian 输入文件

**gs_opt.com** (基态优化):
```bash
#P B3LYP/6-31G* Opt

Benzene - Ground State Optimization

0 1
C     0.0000    1.3970    0.0000
C     1.2098    0.6985    0.0000
C     1.2098   -0.6985    0.0000
C     0.0000   -1.3970    0.0000
C    -1.2098   -0.6985    0.0000
C    -1.2098    0.6985    0.0000
H     0.0000    2.4770    0.0000
H     2.1458    1.2385    0.0000
H     2.1458   -1.2385    0.0000
H     0.0000   -2.4770    0.0000
H    -2.1458   -1.2385    0.0000
H    -2.1458    1.2385    0.0000
```

**gs_freq.com** (基态频率):
```bash
#P B3LYP/6-31G* Freq

Benzene - Ground State Frequency

0 1
C     0.0000    1.3970    0.0000
C     1.2098    0.6985    0.0000
C     1.2098   -0.6985    0.0000
C     0.0000   -1.3970    0.0000
C    -1.2098   -0.6985    0.0000
C    -1.2098    0.6985    0.0000
H     0.0000    2.4770    0.0000
H     2.1458    1.2385    0.0000
H     2.1458   -1.2385    0.0000
H     0.0000   -2.4770    0.0000
H    -2.1458   -1.2385    0.0000
H    -2.1458    1.2385    0.0000
```

**es_opt.com** (激发态优化):
```bash
#P TD(NStates=3) B3LYP/6-31G* Opt

Benzene - Excited State (S1) Optimization

0 1
C     0.0000    1.3970    0.0000
C     1.2098    0.6985    0.0000
C     1.2098   -0.6985    0.0000
C     0.0000   -1.3970    0.0000
C    -1.2098   -0.6985    0.0000
C    -1.2098    0.6985    0.0000
H     0.0000    2.4770    0.0000
H     2.1458    1.2385    0.0000
H     2.1458   -1.2385    0.0000
H     0.0000   -2.4770    0.0000
H    -2.1458   -1.2385    0.0000
H    -2.1458    1.2385    0.0000
```

### 步骤 2: 运行 Gaussian

```bash
# 提交 Gaussian 作业
g16 < gs_opt.com > gs_opt.log
g16 < gs_freq.com > gs_freq.log
g16 < es_opt.com > es_opt.log
```

### 步骤 3: 创建 MOMAP 输入文件

**momap_fluor.inp**:
```bash
&control
    jobtype = 'fluor',
    qc_software = 'gaussian',
    gs_file = 'gs_freq.log',
    es_file = 'es_opt.log',
    temperature = 298.0,
    spectrum_range = 3.0 6.0,
    spectrum_step = 0.01,
/
```

### 步骤 4: 运行 MOMAP

```bash
# 加载 MOMAP 模块
module load momap/2024A-openmpi

# 单核运行
momap < momap_fluor.inp > fluor.out

# 或者并行运行
mpirun -np 8 momap < momap_fluor.inp > fluor.out
```

### 步骤 5: 分析结果

```bash
# 查看辐射速率
grep "Radiative rate" fluor.out

# 绘制光谱
python3 << 'EOF'
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('fluor_spectrum.dat')
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('Energy (eV)')
plt.ylabel('Intensity (a.u.)')
plt.title('Benzene Fluorescence Spectrum')
plt.savefig('benzene_fluor.png', dpi=300)
print("✅ 光谱已保存: benzene_fluor.png")
EOF
```

---

## 示例 2: 内转换 (IC) 速率计算

### 步骤 1: 准备 Gaussian 文件

需要 S0 和 S1 态的频率计算文件：
- `s0_freq.log` - S0 态频率
- `s1_freq.log` - S1 态频率

### 步骤 2: 计算非绝热耦合矩阵元 (NACME)

**nacme.com**:
```bash
#P TD(NStates=3) B3LYP/6-31G* Freq NoSymm

Benzene - NACME Calculation

0 1
[使用 S1 优化的几何结构]
```

### 步骤 3: 创建 MOMAP 输入文件

**momap_ic.inp**:
```bash
&control
    jobtype = 'ic',
    qc_software = 'gaussian',
    gs_file = 's0_freq.log',
    es_file = 's1_freq.log',
    nacme_file = 'nacme.dat',
    temperature = 298.0,
/
```

### 步骤 4: 运行 MOMAP

```bash
module load momap/2024A-openmpi
momap < momap_ic.inp > ic.out

# 提取 IC 速率
grep "IC rate" ic.out
```

---

## 示例 3: 系间交叉 (ISC) 速率计算

### 步骤 1: 准备 Gaussian 文件

需要 S0 和 T1 态的频率计算：
- `s0_freq.log` - S0 态频率
- `t1_freq.log` - T1 态频率

**t1_opt.com** (三重态优化):
```bash
#P UB3LYP/6-31G* Opt

Benzene - Triplet State (T1) Optimization

0 3
C     0.0000    1.3970    0.0000
C     1.2098    0.6985    0.0000
...
```

### 步骤 2: 创建 MOMAP 输入文件

**momap_isc.inp**:
```bash
&control
    jobtype = 'isc',
    qc_software = 'gaussian',
    gs_file = 's0_freq.log',
    es_file = 't1_freq.log',
    spin_orbit = .true.,
    temperature = 298.0,
/
```

### 步骤 3: 运行 MOMAP

```bash
module load momap/2024A-openmpi
momap < momap_isc.inp > isc.out

# 提取 ISC 速率
grep "ISC rate" isc.out
```

---

## 示例 4: 完整的 Slurm 作业脚本

**momap_job.slurm**:
```bash
#!/bin/bash
#SBATCH --job-name=momap_benzene
#SBATCH --partition=X48Cv4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --output=momap_%j.log

# 加载模块
module load gaussian/16
module load momap/2024A-openmpi

# 创建工作目录
WORK_DIR=$SLURM_SUBMIT_DIR
cd $WORK_DIR

# 步骤 1: 运行 Gaussian 计算
echo "===== Step 1: Gaussian Calculations ====="
g16 < gs_opt.com > gs_opt.log
g16 < gs_freq.com > gs_freq.log
g16 < es_opt.com > es_opt.log

# 步骤 2: 运行 MOMAP 荧光计算
echo "===== Step 2: MOMAP Fluorescence ====="
cat > momap_fluor.inp << 'EOF'
&control
    jobtype = 'fluor',
    qc_software = 'gaussian',
    gs_file = 'gs_freq.log',
    es_file = 'es_opt.log',
    temperature = 298.0,
/
EOF

mpirun -np $SLURM_NTASKS momap < momap_fluor.inp > fluor.out

# 步骤 3: 提取结果
echo "===== Step 3: Extract Results ====="
echo "Radiative rate:" > results.txt
grep "Radiative rate" fluor.out >> results.txt

echo "✅ Job completed!"
```

提交作业：
```bash
sbatch momap_job.slurm
```

---

## 常用分析脚本

### 自动提取速率常数

**extract_rates.sh**:
```bash
#!/bin/bash

echo "===== MOMAP Results Summary ====="

if [ -f "fluor.out" ]; then
    echo "Fluorescence:"
    grep "Radiative rate" fluor.out
    echo ""
fi

if [ -f "ic.out" ]; then
    echo "Internal Conversion:"
    grep "IC rate" ic.out
    echo ""
fi

if [ -f "isc.out" ]; then
    echo "Intersystem Crossing:"
    grep "ISC rate" isc.out
    echo ""
fi

# 计算量子产率 (示例)
if [ -f "fluor.out" ] && [ -f "ic.out" ] && [ -f "isc.out" ]; then
    kr=$(grep "Radiative rate" fluor.out | awk '{print $4}')
    kic=$(grep "IC rate" ic.out | awk '{print $3}')
    kisc=$(grep "ISC rate" isc.out | awk '{print $3}')
    
    python3 << EOF
kr = $kr
kic = $kic
kisc = $kisc
qy = kr / (kr + kic + kisc)
print(f"Quantum Yield: {qy:.4f} ({qy*100:.2f}%)")
EOF
fi
```

### 批量绘制光谱

**plot_spectra.py**:
```python
#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import glob

plt.figure(figsize=(12, 8))

for i, filename in enumerate(glob.glob('*_spectrum.dat')):
    data = np.loadtxt(filename)
    label = filename.replace('_spectrum.dat', '')
    plt.plot(data[:, 0], data[:, 1], label=label, linewidth=2)

plt.xlabel('Energy (eV)', fontsize=14)
plt.ylabel('Intensity (a.u.)', fontsize=14)
plt.title('MOMAP Spectra', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('all_spectra.png', dpi=300)
print("✅ 所有光谱已保存: all_spectra.png")
```

---

## 故障排除

### 问题 1: 频率计算有虚频

**症状**: Gaussian 输出中有负频率
**解决**:
```bash
# 检查虚频
grep "Frequencies --" gs_freq.log | grep "-"

# 重新优化几何结构，使用更严格的收敛标准
#P B3LYP/6-31G* Opt=VeryTight Freq
```

### 问题 2: MOMAP 内存不足

**症状**: `Insufficient memory`
**解决**:
```bash
# 增加内存限制
export MOMAP_MEM=8GB

# 或者在输入文件中指定
&control
    memory = 8GB,
/
```

### 问题 3: 相关函数不收敛

**症状**: `Correlation function not converged`
**解决**:
```bash
&control
    max_time = 2000.0,    ! 增加最大时间
    time_step = 0.2,      ! 减小时间步长
    convergence = 1.0E-6, ! 放宽容差
/
```

---

## 工作流程模板

### 完整的光物理性质计算流程

```bash
#!/bin/bash
# complete_photophysics.sh

# 1. 准备分子坐标
# 2. 运行 Gaussian (基态 + 激发态)
# 3. 运行 MOMAP (荧光 + IC + ISC)
# 4. 计算量子产率
# 5. 生成报告

echo "Starting complete photophysics calculation..."

# 步骤 1: Gaussian
g16 < gs_opt.com > gs_opt.log
g16 < gs_freq.com > gs_freq.log
g16 < es_opt.com > es_opt.log
g16 < t1_opt.com > t1_opt.log

# 步骤 2: MOMAP
module load momap/2024A-openmpi

# 荧光
cat > momap_fluor.inp << 'EOF'
&control
    jobtype = 'fluor',
    qc_software = 'gaussian',
    gs_file = 'gs_freq.log',
    es_file = 'es_opt.log',
/
EOF
momap < momap_fluor.inp > fluor.out

# IC
cat > momap_ic.inp << 'EOF'
&control
    jobtype = 'ic',
    qc_software = 'gaussian',
    gs_file = 's0_freq.log',
    es_file = 's1_freq.log',
/
EOF
momap < momap_ic.inp > ic.out

# ISC
cat > momap_isc.inp << 'EOF'
&control
    jobtype = 'isc',
    qc_software = 'gaussian',
    gs_file = 's0_freq.log',
    es_file = 't1_freq.log',
    spin_orbit = .true.,
/
EOF
momap < momap_isc.inp > isc.out

# 步骤 3: 生成报告
python3 << 'PYEOF'
import numpy as np

# 提取速率
with open('fluor.out') as f:
    for line in f:
        if 'Radiative rate' in line:
            kr = float(line.split()[4])

with open('ic.out') as f:
    for line in f:
        if 'IC rate' in line:
            kic = float(line.split()[3])

with open('isc.out') as f:
    for line in f:
        if 'ISC rate' in line:
            kisc = float(line.split()[3])

# 计算量子产率
qy = kr / (kr + kic + kisc)

# 生成报告
report = f"""
===== Photophysics Report =====

Radiative Rate (kr): {kr:.2e} s^-1
IC Rate (kic):       {kic:.2e} s^-1
ISC Rate (kisc):     {kisc:.2e} s^-1
Total Decay Rate:    {kr+kic+kisc:.2e} s^-1

Fluorescence Quantum Yield (ΦF): {qy:.4f} ({qy*100:.2f}%)
Fluorescence Lifetime (τ):       {1/(kr+kic+kisc)*1e9:.2f} ns

================================
"""

print(report)
with open('photophysics_report.txt', 'w') as f:
    f.write(report)
PYEOF

echo "✅ Complete photophysics calculation finished!"
```

运行：
```bash
chmod +x complete_photophysics.sh
./complete_photophysics.sh
```
