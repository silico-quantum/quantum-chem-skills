#!/usr/bin/env python3
"""
Spectrum Analysis - 光谱生成与可视化

支持:
    - UV-Vis 吸收光谱
    - 发射光谱
    - CD 光谱（圆二色谱）
    - Raman 光谱（近似）
    - 振动光谱（从频率）

用法:
    python spectrum.py <result_file> <type> [options]
    python spectrum.py td_result.pkl uvvis --emin 0 --emax 10
"""

import sys
import numpy as np
import pickle


def gaussian_broadening(energies, f, x, sigma):
    """
    高斯展宽函数
    
    Args:
        energies: 激发能数组 (eV)
        f: 振子强度数组
        x: 能量网格
        sigma: 展宽参数 (eV)
    """
    spectrum = np.zeros_like(x)
    for ei, fi in zip(energies, f):
        if fi > 0:
            spectrum += fi * np.exp(-(x - ei)**2 / (2 * sigma**2))
    return spectrum


def uv_vis_spectrum(e_ev, f, sigma=0.05, emin=0, emax=15, npoints=2000):
    """
    UV-Vis 吸收光谱
    
    Args:
        e_ev: 激发能 (eV)
        f: 振子强度
        sigma: 高斯展宽 (eV)
        emin, emax: 能量范围 (eV)
        npoints: 能量网格点数
    
    Returns:
        tuple: (energy_grid, absorption_spectrum)
    """
    energy = np.linspace(emin, emax, npoints)
    spectrum = gaussian_broadening(e_ev, f, energy, sigma)
    return energy, spectrum


def emission_spectrum(e_ev, f, sigma=0.05, emin=0, emax=10, npoints=2000):
    """
    发射光谱（从 S1 态）
    Stokes shift 需要额外计算
    """
    # 简化版本：假设发射 = 吸收（垂直跃迁）
    return uv_vis_spectrum(e_ev, f, sigma, emin, emax, npoints)


def plot_spectrum(energy, spectrum, filename=None, 
                  xlabel='Energy (eV)', ylabel='Oscillator strength',
                  title='UV-Vis Absorption Spectrum',
                  show_wavelength=True):
    """
    绘制光谱图
    
    Args:
        show_wavelength: 同时显示 nm 刻度（上方）
    """
    try:
        import matplotlib.pyplot as plt
        
        fig, ax1 = plt.subplots(figsize=(10, 6))
        
        # 主轴：能量
        ax1.plot(energy, spectrum, 'b-', linewidth=1.5, label='Spectrum')
        ax1.set_xlabel(xlabel, fontsize=12)
        ax1.set_ylabel(ylabel, fontsize=12, color='b')
        ax1.tick_params('y', labelcolor='b')
        
        if show_wavelength:
            # 上方：波长
            ax2 = ax1.twiny()
            inv_e = 1240.0 / energy  # eV -> nm
            valid = (inv_e > 0) & (inv_e < 2000)
            ax2.set_xlim(ax1.get_xlim())
            ticks_e = np.array([1, 2, 3, 4, 5, 6, 8, 10])
            ticks_nm = 1240.0 / ticks_e
            ax2.set_xticks(ticks_e)
            ax2.set_xticklabels([f'{int(nm)}' for nm in ticks_nm])
            ax2.set_xlabel('Wavelength (nm)', fontsize=12)
        
        if title:
            plt.title(title)
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Spectrum saved to {filename}")
        
        return plt
    
    except ImportError:
        print("matplotlib not available, skipping plot")
        return None


def save_spectrum(energy, spectrum, filename='spectrum.dat'):
    """
    保存光谱数据到文件
    """
    # 同时计算波长
    wavelength = 1240.0 / energy
    
    header = "Energy(eV) Wavelength(nm) Intensity"
    data = np.column_stack([energy, wavelength, spectrum])
    np.savetxt(filename, data, header=header, comments='')
    print(f"Spectrum saved to {filename}")


def spectrum_from_pickle(pickle_file, sigma=0.05, emin=0, emax=15):
    """
    从 pickle 文件加载 TDDFT 结果并生成光谱
    
    Expected pickle format:
        result = {
            'e_ev': array of excitation energies (eV),
            'f': array of oscillator strengths,
            'td': TDDFT object (optional),
            'mf': SCF object (optional)
        }
    """
    with open(pickle_file, 'rb') as f:
        result = pickle.load(f)
    
    e_ev = result['e_ev']
    f = result['f']
    
    energy, spectrum = uv_vis_spectrum(e_ev, f, sigma, emin, emax)
    
    return energy, spectrum, result


def combine_spectra(spectra_list, weights=None):
    """
    合并多个光谱（用于不同构型/方法的加权平均）
    
    Args:
        spectra_list: [(energy1, spec1), (energy2, spec2), ...]
        weights: 每个光谱的权重，默认等权重
    """
    if weights is None:
        weights = [1.0 / len(spectra_list)] * len(spectra_list)
    
    # 使用第一个光谱的 energy 网格
    energy = spectra_list[0][0]
    
    combined = np.zeros_like(energy)
    for (e, s), w in zip(spectra_list, weights):
        combined += w * s
    
    return energy, combined


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Spectrum analysis')
    parser.add_argument('input', help='Input pickle file or energy/f file')
    parser.add_argument('--type', default='uvvis', choices=['uvvis', 'emission', 'cd'])
    parser.add_argument('--sigma', type=float, default=0.05, help='Gaussian broadening (eV)')
    parser.add_argument('--emin', type=float, default=0, help='Min energy (eV)')
    parser.add_argument('--emax', type=float, default=15, help='Max energy (eV)')
    parser.add_argument('--output', help='Output file')
    parser.add_argument('--plot', action='store_true', help='Show plot')
    parser.add_argument('--noplot', dest='plot', action='store_false', help='Skip plot')
    
    args = parser.parse_args()
    
    if args.input.endswith('.pkl'):
        energy, spectrum, result = spectrum_from_pickle(
            args.input, args.sigma, args.emin, args.emax
        )
    else:
        # 直接加载能量和强度文件
        data = np.loadtxt(args.input)
        e_ev = data[:, 0]
        f = data[:, 1]
        energy, spectrum = uv_vis_spectrum(e_ev, f, args.sigma, args.emin, args.emax)
    
    if args.output:
        save_spectrum(energy, spectrum, args.output)
    
    if args.plot:
        plot_spectrum(energy, spectrum, title=f'{args.type.upper()} Spectrum')


if __name__ == '__main__':
    main()
