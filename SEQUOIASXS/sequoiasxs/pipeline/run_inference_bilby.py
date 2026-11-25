import os

# ---------------------------------------------------------
# Limit every thread pool to 1 thread.
# ---------------------------------------------------------
os.environ["NUMBA_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# Optional: only if it still hangs with the default layer
# os.environ["NUMBA_THREADING_LAYER"] = "omp"

try:
    import threadpoolctl
    threadpoolctl.threadpool_limits(limits=1, user_api="blas")
    threadpoolctl.threadpool_limits(limits=1, user_api="openmp")
    print("-> Threads limited with threadpoolctl.")
except ImportError:
    pass


import time
import json
import numpy as np
import bilby
import matplotlib.pyplot as plt

from bilby.gw.likelihood import GravitationalWaveTransient
from gwpy.timeseries import TimeSeries
from sequoiasxs.waveforms.dansur.waveform_generator import waveform_generator
from sequoiasxs.waveforms.nrybsur.nrh_waveform_generator import nrh_waveform_generator
from sequoiasxs.waveforms.nrsur7dq4.nrsur_waveform_generator import nrsur_waveform_generator
from sequoiasxs.output.dif import diference


# ---------------------------------------------------------
# 8. Run inference
# ---------------------------------------------------------

def run_inference(ifos, dicc,targ_keys, priors, outdir,folder):
    """
    Run Bayesian parameter estimation for a gravitational-wave event.

    The function selects the waveform model specified in ``dicc``, creates
    the corresponding Bilby waveform generator, constructs a
    ``GravitationalWaveTransient`` likelihood, and runs parameter estimation
    with the Nessai sampler. It then saves the execution time, produces a
    corner plot, reconstructs the maximum-likelihood waveform, saves the
    reconstructed waveform plot, and computes the corresponding waveform
    difference diagnostics.

    INPUTS:
    ------
    
    -ifos : bilby.gw.detector.InterferometerList
        List of interferometers containing the detector data, power spectral
        densities, sampling information, and detector strain data for the
        event.

    -dicc : The created dicctionay in config.py 

    -targ_keys : iterable
        Names of the source parameters used by the waveform generator.
        These parameters define the subset of parameters passed from the
        Bilby parameter dictionary to the waveform model.

    -priors : bilby.core.prior.PriorDict
        Prior distributions for the parameters sampled during the Bayesian
        inference.

    -outdir :
        Directory in which Bilby output files, sampler results, execution
        time information, and diagnostic plots are saved.

    -folder : str
        Identifier of the event or simulation folder. It is used to label
        the reconstructed waveform plot.

    Returns
    -------
    -Result:

        Bilby result object containing the posterior samples, likelihood
        evaluations, evidence estimates, and sampler metadata.

    - A JSON file containing the posterior medians and the lower and upper bounds of their 68% credible intervals, 
    together with the injection parameters of the SXS event.   
    """

    # ---------------------------------------------------------
    # 8.1 Checks if DANSur or an approximant will be used
    # ---------------------------------------------------------

    print('=' * 50)
    print('Loading waveform generator...')
    print(dicc["surrogate_model"].lower())
    if dicc["surrogate_model"].lower() in 'dansur':
        generated_waveform = waveform_generator(dicc,targ_keys)
        label = "DANSur"

    elif dicc["surrogate_model"].lower() in 'NRHybSur3dq8'.lower() :
        generated_waveform = nrh_waveform_generator(dicc,targ_keys)
        label = "NRHybSur3dq8"

    elif dicc["surrogate_model"].lower() in 'NRSur7dq4'.lower() :        
        generated_waveform = nrsur_waveform_generator(dicc,targ_keys)
        label = "NRSur"

    dicc["label"] =label
    print('Waveform generator loaded')
    print('=' * 50)

    # ---------------------------------------------------------
    # 8.2 Runs PE
    # ---------------------------------------------------------
    
    # Force initialization of frequency-domain data and frequency mask

    # for ifo in ifos:
    #         ifo.frequency_domain_strain.shape
    #         ifo.frequency_mask.shape
    #         ifo.frequency_array.shape
    class NRSURLikelihood(GravitationalWaveTransient):

        def log_likelihood_ratio(self, parameters):

            try:
                return super().log_likelihood_ratio(parameters)

            except Exception as e:

                if "omega_ref" in str(e):
                    return -np.inf

                raise
    t0 = time.perf_counter()

    likelihood = NRSURLikelihood(
        interferometers=ifos,
        waveform_generator=generated_waveform,
        priors=priors,
        distance_marginalization=False,
        phase_marginalization=False,
        time_marginalization=False,
        jitter_time=False,
    )
#     likelihood.parameters.update(dicc["inyection"])
# # --- BLOQUE DE DIAGNÓSTICO A ELIMINAR ---
#     likelihood.parameters.update(dicc["inyection"])
#     print("Log-likelihood ratio en parámetros exactos:", likelihood.log_likelihood_ratio(dicc["inyection"].copy()))

#     test_params = dicc["inyection"].copy()
#     phases = np.linspace(0, 2 * np.pi, 36, endpoint=False)
#     llrs = []
#     for ph in phases:
#         test_params["phase"] = float(ph)
#         llrs.append(likelihood.log_likelihood_ratio(test_params.copy()))

#     best_idx = int(np.argmax(llrs))
#     print(f"LLR con phase=0: {llrs[0]:.2f}")
#     print(f"Mejor LLR variando solo phase: {llrs[best_idx]:.2f} (en phase = {phases[best_idx]:.4f} rad)")

#     ifo0 = ifos[0]
#     wf_fd = generated_waveform.frequency_domain_strain(test_params.copy())
#     h_det = ifo0.get_detector_response(wf_fd, test_params.copy())

#     d_td = ifo0.time_domain_strain
#     h_td = bilby.core.utils.infft(h_det, ifo0.sampling_frequency)
#     t_arr = ifo0.time_array

#     print(f"Max |d(t)| H1: {np.max(np.abs(d_td)):.3e} en t = {t_arr[np.argmax(np.abs(d_td))]:.6f}")
#     print(f"Max |h(t)| H1: {np.max(np.abs(h_td)):.3e} en t = {t_arr[np.argmax(np.abs(h_td))]:.6f}")
#     print(f"Diferencia temporal entre picos (ms): {(t_arr[np.argmax(np.abs(d_td))] - t_arr[np.argmax(np.abs(h_td))]) * 1000:.3f} ms")

#     mask = ifo0.frequency_mask
#     psd = ifo0.power_spectral_density_array[mask]
#     df = ifo0.frequency_array[1] - ifo0.frequency_array[0]

#     s_clean_td = dicc["signal"][0]
#     s_clean_fd, _ = bilby.core.utils.nfft(s_clean_td, ifo0.sampling_frequency)
#     mask = ifo0.frequency_mask
#     freqs_band = ifo0.frequency_array[mask]
#     psd_band = ifo0.power_spectral_density_array[mask]
#     df_band = freqs_band[1] - freqs_band[0]

#     def snr_banda(h_fd):
#         return np.sqrt(np.real(
#             4.0 * np.sum(np.conj(h_fd[mask]) * h_fd[mask] / psd_band)
#         ) * df_band)
#     snr_s = snr_banda(s_clean_fd)
#     snr_h = snr_banda(h_det)

#     corr_clean = np.fft.irfft(s_clean_fd * np.conj(h_det) / np.where(ifo0.power_spectral_density_array > 0, ifo0.power_spectral_density_array, np.inf))
#     match_clean = (np.max(np.abs(corr_clean)) * 4.0 * df * len(t_arr)) / (snr_s * snr_h)

#     wf_fd_flip = {"plus": wf_fd["plus"], "cross": -wf_fd["cross"]}
#     h_det_flip = ifo0.get_detector_response(wf_fd_flip, test_params)
#     corr_flip = np.fft.irfft(s_clean_fd * np.conj(h_det_flip) / np.where(ifo0.power_spectral_density_array > 0, ifo0.power_spectral_density_array, np.inf))
#     snr_h_flip = np.sqrt(np.real(ifo0.optimal_snr_squared(h_det_flip)))
#     match_flip = (np.max(np.abs(corr_flip)) * 4.0 * df * len(t_arr)) / (snr_s * snr_h_flip)

#     params_pi = test_params.copy()
#     params_pi["theta_jn"] = np.pi - test_params["theta_jn"]
#     wf_fd_pi = generated_waveform.frequency_domain_strain(params_pi)
#     h_det_pi = ifo0.get_detector_response(wf_fd_pi, params_pi)
#     corr_pi = np.fft.irfft(s_clean_fd * np.conj(h_det_pi) / np.where(ifo0.power_spectral_density_array > 0, ifo0.power_spectral_density_array, np.inf))
#     snr_h_pi = np.sqrt(np.real(ifo0.optimal_snr_squared(h_det_pi)))
#     match_pi = (np.max(np.abs(corr_pi)) * 4.0 * df * len(t_arr)) / (snr_s * snr_h_pi)

#     print(f"Match limpio (normal): {match_clean:.6f}")
#     print(f"Match limpio (cross invertido): {match_flip:.6f}")
#     print(f"Match limpio (theta_jn = pi - theta): {match_pi:.6f}")
#     print(f"Diferencia picos limpios (ms): {(t_arr[np.argmax(np.abs(s_clean_td))] - t_arr[np.argmax(np.abs(h_td))]) * 1000:.3f} ms")

#     dt_grid = np.linspace(-0.025, 0.025, 501)
#     freqs = ifo0.frequency_array
#     overlaps = []
#     for dt_try in dt_grid:
#         h_shifted = h_det * np.exp(-1j * 2.0 * np.pi * freqs * dt_try)
#         ip_clean = 4.0 * df * np.sum(s_clean_fd[mask] * np.conj(h_shifted[mask]) / psd)
#         overlaps.append(ip_clean)

#     overlaps = np.array(overlaps)
#     best_i = int(np.argmax(np.abs(overlaps)))
#     real_match = np.abs(overlaps[best_i]) / (snr_s * snr_h)
#     best_dt_ms = dt_grid[best_i] * 1000.0
#     best_dphi = np.angle(overlaps[best_i]) % (2 * np.pi)

#     print(f"SNR(SXS limpia): {snr_s:.2f} | SNR(DANSur): {snr_h:.2f}")
#     print(f"Match real normalizado [0, 1]: {real_match:.6f}")
#     print(f"Desfase temporal necesario en geocent_time: {best_dt_ms:.4f} ms")
#     print(f"Fase orbital (phase) necesaria: {best_dphi:.4f} rad")

#     test_params["geocent_time"] += dt_grid[best_i]
#     test_params["phase"] = float(best_dphi)
#     print(f"Nuevo LLR con dt y phase corregidos: {likelihood.log_likelihood_ratio(test_params):.2f}")

#     for candidate_phase in [best_dphi / 2.0, (best_dphi / 2.0) + np.pi, (2*np.pi - best_dphi) / 2.0]:
#         p_try = dicc["inyection"].copy()
#         p_try["geocent_time"] += dt_grid[best_i]
#         p_try["phase"] = float(candidate_phase % (2 * np.pi))
#         print(f"LLR con phase={p_try['phase']:.4f}: {likelihood.log_likelihood_ratio(p_try):.2f}")
# #     # --- FIN DEL BLOQUE A ELIMINAR ---
#     raise Exception("Stopping execution after log-likelihood evaluation at injection parameters.")
    dt = time.perf_counter() - t0
    # for i in range(3):
    #     parameters = priors.sample()
    #     start = time.perf_counter()
    #     try:
    #         log_likelihood = likelihood.log_likelihood(parameters)
    #         print(
    #             f"Test draw {i + 1}: log_likelihood={log_likelihood}, "
    #             f"elapsed={time.perf_counter() - start:.2f}s",
    #             flush=True,
    #         )
    #     except Exception as exc:
    #         print(f"Test draw {i + 1} failed: {exc!r}", flush=True)
    start_time = time.time()
        # phase_marginalization=dicc["marginalization"]["phase"],
        # time_marginalization=dicc["marginalization"]["time"],
        # jitter_time=dicc["marginalization"]["jitter"],

        # distance_marginalization=False,
        # phase_marginalization=True,
        # time_marginalization=True,
        # jitter_time=False,
    result=bilby.run_sampler(
                likelihood=likelihood,
                priors=priors,
                sampler="nessai",
                use_ratio=False,
                flow_class="gwflowproposal",
                npoints=dicc["inference"]["npoints"],
                resume=dicc["inference"]["resume"],
                injection_parameters=dicc["inyection"],
                clean=not dicc["inference"]["resume"],
                outdir=outdir,
                label=label,
                npool=1,
                nparallel=1,
                stopping=dicc["inference"]["stopping"],
                sampling_seed=dicc["inference"]["sampling-seed"]
            )
    plt.figure()
    plt.plot(ifos[0].time_array, ifos[0].whitened_time_domain_strain, label='H1 strain')
    #plot waveform generated with top likelihood parameters
    inp = dict(result.posterior.iloc[result.posterior.log_likelihood.idxmax()])
    print(inp)
    plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(inp)['plus'], '--')
    execution_time = time.time() - start_time
    plt.savefig(os.path.join(outdir, label + "_H1_strain.png"))
    timing = {
        "execution_time_seconds": execution_time,
        "execution_time_minutes": execution_time / 60,
        "execution_time_hours":  execution_time / 3600,
    }


    # ---------------------------------------------------------
    # 8.3 Save execution time
    # ---------------------------------------------------------


    title = label + "_execution_time.json"
    time_outdir = os.path.join(outdir, title)
    with open(time_outdir, "w") as f:
            json.dump(timing, f, indent=4)

    print(f"Execution time: {execution_time / 60:.2f} min")

    # ---------------------------------------------------------
    # 8.5 Plots reconstructed waveform and corner plot
    # ---------------------------------------------------------
    result.plot_corner()

    plt.figure(figsize=(10, 5))
    maxll_params = dict(result.posterior[result.posterior.log_likelihood == result.posterior.log_likelihood.max()].iloc[0])
    generated_wave= generated_waveform.time_domain_strain(maxll_params)['plus']
    dicc["generated_waveform"] = generated_wave
    if len(generated_wave) > len(generated_waveform.time_array):
        plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'][:-1], '--')
    else:
        plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'], '--')

    title = "Reconstructed waveform for event " + folder +' with ' + label
    
    plt.title(title, fontsize=11)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, title + '.jpg' ))
    diference(outdir,label,dicc)
    post = bilby.gw.conversion.generate_all_bbh_parameters(result.posterior)
    print("chi_eff:", np.median(post["chi_eff"]), np.percentile(post["chi_eff"], [16, 84]))

    return result


