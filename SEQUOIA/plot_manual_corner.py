from matplotlib import pyplot as plt
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.lines as mlines
from scipy.stats import gaussian_kde
import numpy as np
import bilby
import h5py


def plot_manual_corner(folder, zenodo_file, event_dir, outdir):
        bilby_result = os.path.join(outdir, "NNSur_result.json")
        #result = bilby.result.read_in_result(filename=bilby_result)
        OUTPUT_PLOT = os.path.join(outdir, "corner_comparison_NNSur_vs_Zenodo_" + folder + '.png' )
        zenodo_file = os.path.join(event_dir,zenodo_file)

        plt.rcParams.update({
            "font.size": 12,
            "axes.labelsize": 18,
            "axes.titlesize": 18,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 15,
            "mathtext.fontset": "stix",
            "font.family": "serif",
        })

        def latex_to_plain(lbl):
            mapping = {
                "$q$": "mass_ratio",
                "$\\mathcal{M}$": "chirp_mass",
                "$d_L$": "luminosity_distance",
                "$\\mathrm{DEC}$": "dec",
                "$\\mathrm{RA}$": "ra",
                "$\\theta_{JN}$": "theta_jn",
                "$\\psi$": "psi",
                "$\\chi_1$": "chi_1",
                "$\\chi_2$": "chi_2",
                "$\\chi_1$": "spin_1z",
                "$\\chi_2$": "spin_2z",
            }
            return mapping.get(lbl, lbl.replace("$", ""))

        def plain_symbol(lbl):
            mapping = {
                "mass_ratio": r"q",
                "chirp_mass": r"\mathcal{M} \,[M_\odot]",
                "luminosity_distance": r"d_L \, [Mpc]",
                "dec": r"\delta \,[rad]",
                "ra": r"\alpha \,[rad]",
                "theta_jn": r"\theta_{JN} \,[rad]",
                "psi": r"\psi \,[rad]",
                "spin_1z": r"\chi_1",
                "spin_2z": r"\chi_2",
            }
            return mapping.get(lbl, lbl)

        def median_and_ci(x):
            q16, q50, q84 = np.percentile(np.asarray(x), [16, 50, 84])
            return q50, q84 - q50, q50 - q16

        def kde2d(x, y, xgrid, ygrid):
            kde = gaussian_kde(np.vstack([x, y]))
            X, Y = np.meshgrid(xgrid, ygrid)
            Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
            return X, Y, Z

        def hpd_levels(Z, probs=(0.68, 0.90)):
            z = np.ravel(Z)
            dz = np.sort(z)[::-1]
            cdf = np.cumsum(dz)
            cdf /= cdf[-1]
            levels = []
            for p in probs:
                idx = np.searchsorted(cdf, p)
                levels.append(dz[min(idx, len(dz) - 1)])
            return np.array(levels)

        def force_spines(ax):
            ax.set_frame_on(True)
            ax.patch.set_visible(True)
            ax.patch.set_alpha(0)
            for side in ("top", "right", "bottom", "left"):
                sp = ax.spines[side]
                sp.set_visible(True)
                sp.set_linewidth(1.0)
                sp.set_color("0.55")
                sp.set_zorder(20)
            ax.xaxis.set_ticks_position("bottom")
            ax.yaxis.set_ticks_position("left")
            ax.tick_params(axis="both", which="both", labelsize=10, width=0.6, color="0.5")
            ax.set_axisbelow(True)
            ax.minorticks_on()
            ax.grid(True, which="major", color="0.75", linestyle="-", linewidth=0.8, alpha=0.8)
            ax.grid(True, which="minor", color="0.88", linestyle=":", linewidth=0.5, alpha=0.7)
        def fmt_unc(label, m, up, lo):
            return rf"{plain_symbol(label)} = {m:.2f}^{{+{up:.2f}}}_{{-{lo:.2f}}}"

        def x_formatter(x, pos):
            if abs(x) < 1:
                return f"{x:.1f}"
            return f"{x:.0f}"

        def plot_manual_corner(samples_a, samples_b, labels, output_path):
            ndim = len(labels)
            fig, axes = plt.subplots(ndim, ndim, figsize=(24, 12), dpi=300)

            
            def fix_levels(levels):
                levels = np.sort(levels)
                levels = np.unique(levels)

                if len(levels) < 2:
                    return None

                # asegurar estrictamente creciente
                if np.any(np.diff(levels) <= 0):
                    levels = levels + np.linspace(0, 1e-12, len(levels))

                return levels

            for i in range(ndim):
                for j in range(ndim):
                    ax = axes[i, j]
                    force_spines(ax)

                    if i < j:
                        ax.axis("off")
                        continue

                    xlab = labels[j]

                    if i == j:
                        xa = samples_a[:, j]
                        xb = samples_b[:, j]

                        ax.hist(xb, bins=32, density=True, histtype="step", color="k", lw=1.1)
                        ax.hist(xa, bins=32, density=True, histtype="step", color="#1f77b4", lw=1.1)

                        mb, upb, lob = median_and_ci(xb)
                        ma, upa, loa = median_and_ci(xa)
                        #Cambiar tamaño numeros
                        ax.text(
                            0.03, 1.30,
                            f"${fmt_unc(xlab, mb, upb, lob)}$",
                            transform=ax.transAxes,
                            ha="left", va="bottom",
                            fontsize=15, color="k", clip_on=False
                        )

                        ax.text(
                            0.03, 1.05,
                            f"${fmt_unc(xlab, ma, upa, loa)}$",
                            transform=ax.transAxes,
                            ha="left", va="bottom",
                            fontsize=15, color="#1f77b4", clip_on=False
        )

                    else:
                        xb = samples_b[:, j]
                        yb = samples_b[:, i]
                        xa = samples_a[:, j]
                        ya = samples_a[:, i]

                        x_min = min(xa.min(), xb.min())
                        x_max = max(xa.max(), xb.max())
                        y_min = min(ya.min(), yb.min())
                        y_max = max(ya.max(), yb.max())

                        padx = 0.05 * (x_max - x_min)
                        pady = 0.05 * (y_max - y_min)

                        xgrid = np.linspace(x_min - padx, x_max + padx, 20)
                        ygrid = np.linspace(y_min - pady, y_max + pady, 20)

                        X, Y, Zb = kde2d(xb, yb, xgrid, ygrid)
                        _, _, Za = kde2d(xa, ya, xgrid, ygrid)

                        #  evitar KDE degenerado
                        if np.allclose(Zb, Zb[0, 0]) or np.allclose(Za, Za[0, 0]):
                            continue

                        #  aplicar fix a ambos
                        levels_b = fix_levels(hpd_levels(Zb, probs=(0.68, 0.90)))
                        levels_a = fix_levels(hpd_levels(Za, probs=(0.68, 0.90)))

                        if levels_b is not None:
                            ax.contour(X, Y, Zb, levels=levels_b, colors="k", linewidths=1.05)

                        if levels_a is not None:
                            ax.contour(X, Y, Za, levels=levels_a, colors="#1f77b4", linewidths=1.05)

                        # medianas
                        xm_b = np.median(xb)
                        ym_b = np.median(yb)
                        xm_a = np.median(xa)
                        ym_a = np.median(ya)

                        ax.axvline(xm_b, color="k", lw=0.8, alpha=0.9)
                        ax.axhline(ym_b, color="k", lw=0.8, alpha=0.9)

                        ax.axvline(xm_a, color="#1f77b4", lw=0.8, alpha=0.9)
                        ax.axhline(ym_a, color="#1f77b4", lw=0.8, alpha=0.9)

                        ax.plot(xm_b, ym_b, "o", color="k", ms=2.8)
                        ax.plot(xm_a, ym_a, "o", color="#1f77b4", ms=2.8)

                    if i < ndim - 1:
                        ax.set_xticklabels([])
                    else:
                        ax.set_xlabel(rf"${plain_symbol(xlab)}$", fontsize=20)

                    if j > 0:
                        ax.set_yticklabels([])
                    else:
                        ax.set_ylabel(rf"${plain_symbol(labels[i])}$", fontsize=20)

                    # asegurar spines visibles
                    for side in ("top", "right", "bottom", "left"):
                        ax.spines[side].set_visible(True)
                        ax.spines[side].set_linewidth(1.0)
                        ax.spines[side].set_color("0.55")
                        ax.spines[side].set_zorder(20)

            handles = [
                mlines.Line2D([], [], color="k", lw=8, label="LVK samples"),
                mlines.Line2D([], [], color="#1f77b4", lw=8, label="DANSur"),
            ]
            
            fig.legend(handles=handles, loc="upper right", frameon=False, fontsize=20, bbox_to_anchor=(0.8, 0.94))
            fig.suptitle(str(folder), fontsize=24, x=0.4, y=0.99, ha="center")

            plt.subplots_adjust(left=0.065, right=0.82, bottom=0.065, top=0.94, wspace=0.26, hspace=0.1)
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            plt.close(fig)

            print(f"Guardado en: {output_path}")
        def main():
            result = bilby.result.read_in_result(bilby_result)
            samples_bilby = result.samples
            labels_bilby_plain = [latex_to_plain(lbl) for lbl in result.parameter_labels]

            with h5py.File(zenodo_file, "r") as f:
                try:
                    structured = f["C01:Mixed/posterior_samples"][:]
                except KeyError:
                    group_name = None
                    for key in f.keys():
                        if key.startswith("C00:Mixed"):
                            group_name = key
                            break
                    structured = f[f"{group_name}/posterior_samples"][:]
            labels_zenodo = list(structured.dtype.names)

            requested = [
                "mass_ratio",
                "chirp_mass",  #Evitar el source (para no tener en cuenta la expansión del universo)
                "luminosity_distance",
                "dec",
                "ra",
                "theta_jn",
                "psi",
                "spin_1z",
                "spin_2z",
            ]

            zenodo_map = {
                "chi_1": "spin_1z",
                "chi_2": "spin_2z",
            }

            common = []
            for p in requested:
                if p in labels_bilby_plain:
                    if p in zenodo_map:
                        if zenodo_map[p] in labels_zenodo:
                            common.append(p)
                    elif p in labels_zenodo:
                        common.append(p)

            if len(common) < 2:
                raise ValueError(
                    f"No hay suficientes parámetros comunes. "
                    f"Bilby={labels_bilby_plain}, Zenodo={labels_zenodo}"
                )

            idx_b = []
            idx_z = []
            for p in common:
                idx_b.append(labels_bilby_plain.index(p))
                if p in zenodo_map:
                    idx_z.append(list(labels_zenodo).index(zenodo_map[p]))
                else:
                    idx_z.append(list(labels_zenodo).index(p))

            samples_nnsur = samples_bilby[:, idx_b]
            samples_zenodo = np.column_stack([structured[labels_zenodo[i]] for i in idx_z])
            plot_manual_corner(samples_nnsur, samples_zenodo, common, OUTPUT_PLOT)

        main()
