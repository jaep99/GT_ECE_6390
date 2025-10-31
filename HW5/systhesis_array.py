# Disk-like planar array synthesis with radial Kaiser taper
# Three separate figures for D = 8λ, 12λ, 16λ:

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


cmap_w2b = LinearSegmentedColormap.from_list(
    'white_to_blue', ['#ffffff', '#9ec5fe', '#2b78e4'], N=256
)


def radial_kaiser(rho, R, beta):
    """Radial Kaiser–Bessel taper on a disk of radius R."""
    w = np.zeros_like(rho, dtype=float)
    inside = rho <= R + 1e-12
    nu = np.sqrt(np.maximum(0.0, 1.0 - (rho[inside] / R) ** 2))
    w[inside] = np.i0(beta * nu) / np.i0(beta)
    return w


def array_factor_cut(theta_deg, X, Y, W, k, phi, chunk=5000):
    """Array factor along a principal-plane cut at azimuth phi."""
    th = np.deg2rad(theta_deg)
    kx = k * np.cos(phi) * np.sin(th)
    ky = k * np.sin(phi) * np.sin(th)
    Nt = th.size
    AF = np.zeros(Nt, dtype=complex)
    denom = W.sum()

    for s in range(0, Nt, chunk):
        e = min(s + chunk, Nt)
        phase = np.outer(X, kx[s:e]) + np.outer(Y, ky[s:e])  # (Ne x ns)
        AF[s:e] = (W[:, None] * np.exp(1j * phase)).sum(axis=0) / denom
    return AF


def compute_SLL(P, theta_deg, D, lam):
    """Peak sidelobe level (dB) w.r.t. mainlobe, mainlobe excluded."""
    P = np.asarray(P).ravel()
    th = np.asarray(theta_deg).ravel()

    theta_null = 1.22 * lam / D * 180 / np.pi
    exclude = np.abs(th) <= max(1.5 * theta_null, 2.0)

    Q = P.copy()
    Q[exclude] = -np.inf

    loc = np.where((Q[1:-1] > Q[:-2]) & (Q[1:-1] >= Q[2:]))[0] + 1
    if loc.size == 0:
        return -np.inf
    return 10 * np.log10(Q[loc].max())


def main():
    np.random.seed(0)

    lam = 1.0
    d = lam / 2
    D_list = [8, 12, 16] # diameters in wavelengths
    targetSLL_dB = -22
    phi_deg = 0
    phi = np.deg2rad(phi_deg)
    theta_deg = np.linspace(-90, 90, 40001)
    k = 2 * np.pi / lam
    fig_size = (15, 4.8)

    for Dlam in D_list:
        D = Dlam * lam
        R = D / 2

        # Element grid (centers) and disk mask
        nspan = int(np.floor(R / d))
        g = np.arange(-nspan, nspan + 1) * d
        Xg, Yg = np.meshgrid(g, g, indexing='xy')
        Rho = np.hypot(Xg, Yg)
        mask = Rho <= R + 1e-12

        Xv = Xg[mask].astype(float)
        Yv = Yg[mask].astype(float)

        # Smallest beta that satisfies target SLL
        beta_opt = 12.0
        for b in np.arange(0.0, 12.0 + 1e-9, 0.25):
            W_full = radial_kaiser(Rho, R, b)
            Wv = W_full[mask].astype(float)
            AF = array_factor_cut(theta_deg, Xv, Yv, Wv, k, phi)
            P = np.abs(AF) ** 2
            if compute_SLL(P, theta_deg, D, lam) <= targetSLL_dB:
                beta_opt = b
                break

        # Final weights and patterns
        W_full = radial_kaiser(Rho, R, beta_opt)
        Wv = W_full[mask].astype(float)
        AF = array_factor_cut(theta_deg, Xv, Yv, Wv, k, phi)
        AF /= np.max(np.abs(AF))
        P = np.abs(AF) ** 2
        SLL_dB = compute_SLL(P, theta_deg, D, lam)

        # Figure: [aperture | principal-plane | polar]
        fig = plt.figure(figsize=fig_size, layout='constrained')

        # (1) Aperture weights (tight squares via pcolormesh)
        ax1 = fig.add_subplot(1, 3, 1)
        xe = np.arange(g[0] - d / 2, g[-1] + d / 2 + 1e-12, d)
        ye = np.arange(g[0] - d / 2, g[-1] + d / 2 + 1e-12, d)
        W_plot = W_full.copy()
        W_plot[~mask] = np.nan

        pcm = ax1.pcolormesh(
            xe / lam, ye / lam, W_plot, cmap=cmap_w2b, vmin=0.0, vmax=1.0, shading='flat'
        )
        circ = plt.Circle((0, 0), R / lam, fc='none', ec=(0.6, 0, 0), lw=1.0)
        ax1.add_patch(circ)
        ax1.set_aspect('equal', 'box')
        ax1.grid(True, ls=':', alpha=0.35)
        ax1.set_xlim([(g[0] - d / 2) / lam, (g[-1] + d / 2) / lam])
        ax1.set_ylim([(g[0] - d / 2) / lam, (g[-1] + d / 2) / lam])
        cb = fig.colorbar(pcm, ax=ax1, fraction=0.046, pad=0.04, ticks=[1.0, 0.5, 0.25])
        cb.set_label('Relative weighting')
        cb.ax.set_yticklabels(['1.00', '0.50', '0.25'])
        ax1.set_title(f'D = {Dlam}λ — Kaiser β = {beta_opt:.2f}')

        # (2) Principal-plane ||AF||^2 (dB)
        ax2 = fig.add_subplot(1, 3, 2)
        ax2.plot(theta_deg, 10 * np.log10(np.maximum(P, np.finfo(float).tiny)), lw=1.5)
        ax2.grid(True, ls=':')
        ax2.set_xlim([-90, 90])
        ax2.set_ylim([-60, 0])
        ax2.axhline(targetSLL_dB, ls='--', color='0.5')
        ax2.text(-88, targetSLL_dB + 0.8, f'Target = {targetSLL_dB} dB',
                 ha='left', va='bottom', color='0.35')
        ax2.set_xlabel(r'$\theta$ (deg)')
        ax2.set_ylabel(r'$||AF||^2$ (dB)')
        ax2.set_title(f'Principal plane (φ = {phi_deg}°) — SLL = {SLL_dB:.1f} dB')

        # (3) Polar (dB), -90..90, 0° at North, clockwise
        ax3 = fig.add_subplot(1, 3, 3, projection='polar')
        dB_floor = -60.0
        PdB = 10 * np.log10(np.maximum(P, np.finfo(float).tiny))
        PdB = np.maximum(PdB, dB_floor)

        ax3.plot(np.deg2rad(theta_deg), PdB, lw=1.2)
        ax3.set_theta_zero_location('N')
        ax3.set_theta_direction(-1)
        ax3.set_thetalim(np.deg2rad([-90, 90]))
        ax3.set_thetagrids(np.arange(-90, 91, 30))
        ax3.set_rlim(dB_floor, 0)
        ax3.set_rticks([-60, -50, -40, -30, -20, -10, 0])
        ax3.set_rlabel_position(135)
        ax3.set_title(f'Polar — D = {Dlam}λ (dB)', pad=14)

        print(f'D = {Dlam:2d} λ:  β ≈ {beta_opt:.2f}  →  SLL = {SLL_dB:.2f} dB')

    plt.show()


if __name__ == '__main__':
    main()
