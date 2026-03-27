# Gaussian 高级使用技巧

## 🎯 计算类型详解

### # (Route Section)

`#` 开启 Route Section，所有关键词都在这里定义。

```
#P - 详细输出
#N - 标准输出（默认）
#T - 简洁输出
```

---

## ⚛️ 高级方法

### CASSCF (Complete Active Space SCF)

```
CASSCF(N,M)
```

- N: 活性空间电子数
- M: 活性空间轨道数

**示例:**
```
# CASSCF(6,6)/6-31G*
```

**适用:**
- 多参考态系统
- 激发态优化
- 化学键断裂

### CBS Methods (Complete Basis Set)

高精度能量计算，达到化学精度 (±1 kcal/mol)。

| 方法 | 说明 |
|------|------|
| CBS-4M | 快速，中等精度 |
| CBS-QB3 | 推荐，高精度 |
| CBS-APNO | 最高精度 |

```
# CBS-QB3
```

### Gn Methods (Gaussian-n)

| 方法 | 精度 |
|------|------|
| G2 | ~2 kcal/mol |
| G3 | ~1 kcal/mol |
| G4 | ~1 kcal/mol |

```
# G4
```

### W1 Methods

Weizmann-1，极高精度。

```
# W1U
```

### CC Methods (Coupled Cluster)

| 方法 | 说明 | 计算量 |
|------|------|--------|
| CCD | 耦合簇双激发 | N⁶ |
| CCSD | 单+双激发 | N⁶ |
| CCSD(T) | 加三重微扰 | N⁷ |
| EOMCCSD | 激发态 | N⁶ |

```
# CCSD(T)/cc-pVTZ
```

---

## 🔬 自定义基组 (Gen/GenECP)

### Gen - 通用基组输入

```
# B3LYP/Gen

Title
0 1
C 0.0 0.0 0.0
H 0.0 0.0 1.09

C 0
6-31G*
****
H 0
6-31G**
****
```

### GenECP - 带赝势

```
# B3LYP/GenECP

Title
0 1
Pt 0.0 0.0 0.0

Pt 0
SDD
****
```

### ChkBasis - 从检查点读取基组

```
# B3LYP/ChkBasis Geom=AllCheck
```

---

## 🌈 激发态高级选项

### TD (Time-Dependent)

```
TD(NStates=10,Root=1,EqSolv,Singlets)
```

| 选项 | 说明 |
|------|------|
| NStates=N | 计算N个态 |
| Root=N | 优化第N个态 |
| Singlets | 单重态 |
| Triplets | 三重态 |
| EqSolv | 平衡溶剂 |
| NonEqSolv | 非平衡溶剂 |
| IVO | 改进的虚轨道 |

### CIS (Configuration Interaction Singles)

```
CIS(NStates=10,Triplets)
```

### SAC-CI (Symmetry-Adapted Cluster CI)

```
SAC-CI(NSinglet=10,NDoublet=5)
```

---

## 💧 溶剂效应详解 (SCRF)

### 模型类型

| 模型 | 说明 |
|------|------|
| PCM | 极化连续介质模型 |
| IEFPCM | 积分方程形式 PCM |
| CPCM | 导电边界 PCM |
| SMD | 溶化模型（推荐） |
| COSMO | 导体屏蔽模型 |

### 选项

```
SCRF=(SMD,Solvent=Water,Read)
```

| 选项 | 说明 |
|------|------|
| Solvent=Name | 溶剂名称 |
| Read | 自定义参数 |
| DoVacuum | 真空计算 |
| NoDielUpdate | 不更新介电常数 |

### 显式+隐式溶剂

```
# B3LYP/6-31G* SCRF=(SMD,Solvent=Water)

Title
0 1
溶质坐标
--Link1--
# B3LYP/6-31G* SCRF=(SMD,Solvent=Water)

Title
0 1
溶质 + 显式水分子坐标
```

---

## 🎛️ SCF 高级选项

### 收敛算法

| 选项 | 说明 |
|------|------|
| Conventional | 传统SCF |
| Direct | 直接SCF（不存积分） |
| InCore | 全部存内存 |
| QC | 二次收敛 |
| XQC | 先常规后QC |
| YQC | 更激进的QC |

### 收敛标准

```
SCF=Tight      # 10⁻⁸
SCF=VeryTight  # 10⁻¹⁰
SCF=Sleazy     # 10⁻⁴ (测试用)
```

### 阻尼和DIIS

```
SCF=NoDIIS     # 禁用DIIS
SCF=DM         # 阻尼方法
SCF=DDIIS      # 阻尼+DIIS
```

### Fermi展布（金属）

```
SCF=Fermi
```

### 最大迭代次数

```
SCF=MaxCycle=500
```

---

## 🎲 初始猜测高级选项

### Guess类型

| 选项 | 说明 |
|------|------|
| Huckel | 扩展Hückel（默认） |
| Harris | Harris泛函 |
| Read | 从chk读取 |
| Mix | 混合HOMO/LUMO |
| Alter | 手动交换轨道 |
| Only | 只算猜测 |
| Cards | 从输入读取 |
| Sad | 超分子近似 |

### 轨道交换

```
Guess=Alter

Title
0 1
coordinates

8,9    ! 交换轨道8和9
0       ! 结束
```

---

## 📐 几何输入高级选项

### 坐标来源

| 选项 | 说明 |
|------|------|
| AllCheck | 从chk读取所有 |
| Check | 从chk读取坐标 |
| Modify | 修改坐标 |
| ZMatrix | Z-矩阵格式 |
| Connectivity | 指定连接性 |
| Step=n | 从第n步继续 |

### 冻结原子

```
Geom=ModRedundant

Title
0 1
coordinates

B 1 2 F     ! 冻结键1-2
A 1 2 3 F   ! 冻结角1-2-3
D 1 2 3 4 F ! 冻结二面角
```

### Z-矩阵示例

```
Geom=ZMatrix

Title
0 1
C
C 1 CC
H 1 CH 2 HCC
H 1 CH 2 HCC 3 180.0

CC=1.40
CH=1.09
HCC=120.0
```

---

## 📊 布居分析详解

### Mulliken (默认)

```
Pop=Regular
```

### Natural Population Analysis

```
Pop=NPA
Pop=NBO
Pop=NBODel
```

### ESP拟合电荷

| 方法 | 说明 |
|------|------|
| MK | Merz-Kollman |
| CHelp | CHelp |
| CHelpG | CHelpG |

```
Pop=MK IOp(6/33=2)
```

### Hirshfeld

```
Pop=Hirshfeld
```

### 原子偶极矩

```
Pop=Dipole
```

---

## 🔄 对称性

| 选项 | 说明 |
|------|------|
| Symm | 使用对称性（默认） |
| NoSymm | 禁用对称性 |
| Loose | 宽松对称性 |
| None | 不使用对称性 |
| Int | 内坐标对称性 |

**何时禁用对称性:**
- 过渡态搜索
- 频率计算出现虚频
- SCF不收敛
- 分子畸变

---

## 🧪 特殊关键词

### NMR

```
NMR=GIAO    # 规范不变原子轨道
NMR=CSGT    # 连续集规范变换
```

### EPT (Electron Propagator Theory)

```
EPT=IP   # 电离能
EPT=EA   # 电子亲和能
```

### Counterpoise (BSSE校正)

```
Counterpoise=2   # 2个片段
```

### ONIOM

```
ONIOM(B3LYP/6-31G*:HF/3-21G:PM6)
```

### ADMP/BOMD (分子动力学)

```
ADMP(B3LYP/6-31G*) Step=1000
BOMD(B3LYP/6-31G*) Step=500
```

### IRCmax

```
IRCmax(Step=20,MaxPoints=50)
```

---

## ⚡ IOp (内部选项)

常用IOp列表：

### DFT网格

```
IOp(3/76=1000000000)  # UltraFine网格
IOp(3/76=10000000)    # FineGrid
IOp(3/77=0)           # 剪裁网格
```

### 频率

```
IOp(6/7=1)    # 归一化频率
IOp(6/7=2)    # 不归一化
```

### 输出

```
IOp(7/33=2)   # 输出分子轨道图形
IOp(8/10=1)   # 输出更多信息
```

### 积分

```
IOp(3/5=1)    # 使用单电子积分近似
IOp(3/7=1)    # 使用双电子积分近似
```

---

## 📋 输出控制

### GFPrint/GFInput

```
GFPrint   # 打印基组
GFInput   # 输出基组到文件
```

### Punch

```
Punch=MO   # 输出分子轨道
Punch=Density  # 输出密度矩阵
```

### Output

```
Output=Pickett
```

---

## 🔍 稳定性测试

```
Stable=Opt   # 测试并优化波函数
Stable=NoOpt # 只测试
```

---

## 🌡️ 温度和压力

```
Temperature=298.15
Pressure=1.0
```

频率分析后自动计算热力学性质。

---

## 🧲 外场

```
Field=X+100  # X方向100 a.u.电场
Field=Z-50   # Z方向-50 a.u.电场
```

---

## 实用技巧

### 1. SCF不收敛

```
# B3LYP/6-31G* SCF=XQC Guess=Mix
```

### 2. 优化不收敛

```
# Opt=CalcFC MaxCycle=200
```

### 3. 过渡态优化

```
# Opt=TS CalcFC NoEigenTest
```

### 4. 开壳层计算

```
# UB3LYP/6-31G* Guess=Mix
```

### 5. 阴离子/激发态

```
# B3LYP/6-31+G**  # 扩散函数
```

### 6. 弱相互作用

```
# ωB97X-D/def2-TZVP  # 色散修正
```

### 7. 高精度单点

```
# CCSD(T)/cc-pVTZ Geom=AllCheck Guess=Read
```

### 8. 溶剂中优化

```
# Opt B3LYP/6-31G* SCRF=(SMD,Solvent=Water)
```

---

*基于 Gaussian 16 Rev. C.01 官方关键词列表*
*最后更新: 2026-03-10*
