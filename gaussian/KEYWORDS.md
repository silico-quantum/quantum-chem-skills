# Gaussian 关键词完整参考

## 计算类型 (Job Type)

| 关键词 | 说明 | 常用选项 |
|--------|------|----------|
| SP | 单点能量 | 默认 |
| Opt | 几何优化 | Opt=CalcFC, Opt=TS |
| Freq | 频率分析 | Freq=Raman |
| Scan | 势能面扫描 | Scan=Opt |
| IRC | 内禀反应坐标 | IRC=CalcFC |
| Polar | 极化率 | Polar=EnOnly |

## DFT 泛函

| 泛函 | 类型 | 特点 |
|------|------|------|
| B3LYP | 杂化 | 最常用 (20% HF) |
| PBE0 | 杂化 | 25% HF |
| M06-2X | Meta-GGA | 54% HF |
| CAM-B3LYP | 范围分离 | 长程修正 |
| ωB97X-D | 范围分离+色散 | 长程+色散 |

## 基组

### Pople 系列
- 6-31G: 双ζ
- 6-31G*: 双ζ+d极化
- 6-31G**: 双ζ+dp极化
- 6-31+G*: 双ζ+扩散+d
- 6-311G**: 三ζ+极化

### 相关一致
- cc-pVDZ, cc-pVTZ, cc-pVQZ
- aug-cc-pVTZ (加扩散)

### def2 系列
- def2-SVP, def2-TZVP, def2-QZVP

## TDDFT 激发态

```
TD(NStates=10,Root=1,Singlets)
```

- NStates: 激发态数量
- Root: 优化哪个态
- Singlets/Triplets: 单/三重态

## 溶剂效应 (SCRF)

```
SCRF=(SMD,Solvent=Water)
```

常用溶剂: Water, Ethanol, Acetonitrile, DMF, DMSO, Toluene

## Opt 选项

- Opt=CalcFC: 初始计算力常数
- Opt=TS: 过渡态
- Opt=Tight: 紧收敛
- Opt=MaxStep=N: 最大步长

## SCF 收敛

- SCF=QC: 二次收敛（最稳定）
- SCF=XQC: 先快后慢
- SCF=Tight: 紧收敛

## 布居分析 (Pop)

- Pop=Full: 完整分析
- Pop=NBO: Natural Bond Orbital
- Pop=NPA: Natural Population
- Pop=MK: Merz-Kollman 电荷

## Link 0 命令

- %chk=file.chk: 检查点
- %mem=60GB: 内存
- %nproc=64: 核心数

