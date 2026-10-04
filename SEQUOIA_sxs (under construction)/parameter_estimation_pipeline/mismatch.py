import numpy as np
import matplotlib.pyplot as plt
def inner_product(h, hhat, fs):
        """
        <h | hhat> = 4 Re integral h(f) hhat*(f) df

        No PSD weighting.
        """

        # print("ANTES DE ASARRAY:")
        # print("h    =", h)
        # print("hhat =", hhat)
        # print("type(h) =", type(h))

        h = np.asarray(h)
        hhat = np.asarray(hhat)

        # print("DESPUÉS DE ASARRAY:")
        # print("h =", h)
        # print("h.dtype =", h.dtype)
        # print("h.shape =", h.shape)
        if len(h) != len(hhat):
            print(len(h))
            print(len(hhat))
            raise ValueError("h and hhat must have the same length")

        N = len(h)

        # Fourier transforms
        hf = np.fft.rfft(h)
        hhat_f = np.fft.rfft(hhat)

        # Frequency resolution
        df = fs / N

        # Inner product
        integrand = hf * np.conj(hhat_f)

        return 4.0 * np.real(np.sum(integrand) * df)


def overlap(h, hhat, fs):
        """
        O(h, hhat) =
            <h|hhat> /
            sqrt(<h|h><hhat|hhat>). Uses the expression from: 	arXiv:2412.06946
        """
        zero_h= np.flatnonzero(h)
        zero_hat= np.flatnonzero(hhat)
        print(zero_h[0] , zero_hat[0])
        if zero_h[0]>zero_hat[0]:
            h = h[zero_h[0]:]
            hhat = hhat[zero_h[0]:]
        else:
            h = h[zero_hat[0]:]
            hhat = hhat[zero_hat[0]:]
              
        hh = inner_product(h, h, fs)
        hhat_hhat = inner_product(hhat, hhat, fs)
        h_hhat = inner_product(h, hhat, fs)

        denominator = np.sqrt(hh * hhat_hhat)
        h_array = np.asarray(
            h.value if hasattr(h, "value") else h
        )

        hhat_array = np.asarray(
            hhat.value if hasattr(hhat, "value") else hhat
        )

        if denominator == 0:
            raise ValueError("One waveform has zero norm")

        O = h_hhat / denominator

        if hasattr(h, "times"):
            t = h.times.value
        else:
            t = np.arange(len(h_array)) / fs
        plt.figure(figsize=(10, 5))

        plt.plot(t, h_array, label="h sxs")
        plt.plot(t, hhat_array, label=r"$\hat{h}$ (dansur)",
                 linestyle="--")

        plt.xlabel("Time [s]")
        plt.ylabel("Strain")
        plt.title(f"Waveforms — Overlap = {O:.6f}")

        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        plt.savefig('overlap 22 mode')

        return h_hhat / denominator


def L2(h, hhat, fs):
        """
        L2(h, hhat) = 1 - O(h, hhat)
        """

        return 1.0 - overlap(h, hhat, fs)
