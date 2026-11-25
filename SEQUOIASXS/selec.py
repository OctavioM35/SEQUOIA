import ast
import re
import sxs
import numpy as np
import pandas as pd

df = sxs.load("dataframe")

def parse_sxs_num(sxs_id):
    m = re.search(r"SXS:BBH:(\d+)", str(sxs_id))
    return int(m.group(1)) if m else 999999

def parse_ecc_strict(val):
    # Rechazar NaN, None o vacíos: exigir medición explícita de baja excentricidad
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return np.inf
    s = str(val).strip().lower()
    if s in ("nan", "none", ""):
        return np.inf
    if s.startswith("<"):
        return float(s[1:].strip())
    try:
        return float(s)
    except ValueError:
        return np.inf

def get_spin_vec(row, prefix, idx):
    col_z = f"{prefix}{idx}_z"
    if col_z in row and not pd.isna(row[col_z]):
        return np.array([
            float(row.get(f"{prefix}{idx}_x", 0.0)),
            float(row.get(f"{prefix}{idx}_y", 0.0)),
            float(row[col_z]),
        ])
    val = row.get(f"{prefix}{idx}", None)
    if isinstance(val, str):
        try:
            val = ast.literal_eval(val)
        except Exception:
            return None
    if isinstance(val, (list, tuple, np.ndarray)) and len(val) == 3:
        arr = np.asarray(val, dtype=float)
        return None if np.any(np.isnan(arr)) else arr
    return None

sxs_nums = [parse_sxs_num(idx) for idx in df.index]
ecc = df["reference_eccentricity"].apply(parse_ecc_strict)
q_ref = df["reference_mass_ratio"].astype(float)
q_norm = np.where(q_ref <= 1.0, q_ref, 1.0 / q_ref)

valid_rows = []
for i, (sxs_id, row) in enumerate(df.iterrows()):
    # 1. Solo catálogo clásico estándar (<= 2265) y excentricidad explícita <= 1e-3
    if sxs_nums[i] > 2265 or ecc.iloc[i] > 1e-3:
        continue

    # 2. Extraer espines tanto iniciales como de referencia
    s1_ref = get_spin_vec(row, "reference_dimensionless_spin", 1)
    s2_ref = get_spin_vec(row, "reference_dimensionless_spin", 2)
    s1_ini = get_spin_vec(row, "initial_dimensionless_spin", 1)
    s2_ini = get_spin_vec(row, "initial_dimensionless_spin", 2)

    if s1_ref is None or s2_ref is None:
        continue

    # 3. Cero precesión estricto (componentes en el plano orbital < 1e-5)
    in_plane_ref = max(np.hypot(s1_ref[0], s1_ref[1]), np.hypot(s2_ref[0], s2_ref[1]))
    in_plane_ini = 0.0
    if s1_ini is not None and s2_ini is not None:
        in_plane_ini = max(np.hypot(s1_ini[0], s1_ini[1]), np.hypot(s2_ini[0], s2_ini[1]))

    if max(in_plane_ref, in_plane_ini) >= 1e-5:
        continue

    # 4. Rango de validez de DANSur: q >= 0.125 y |chi_z| <= 0.8
    q = float(q_norm[i])
    c1, c2 = float(s1_ref[2]), float(s2_ref[2])
    if not (0.125 <= q <= 1.0 and abs(c1) <= 0.8 and abs(c2) <= 0.8):
        continue

    # 5. Regímenes óptimos para romper degeneración de espines:
    #    - Asimétrico con espines distintos: 0.25 <= q <= 0.65 y |chi1 - chi2| >= 0.35
    #    - O sin espines (|chi1|, |chi2| < 1e-3) hasta q <= 0.85
    is_asym_spin = (0.25 <= q <= 0.65) and (abs(c1 - c2) >= 0.35)
    is_non_spin = (0.25 <= q <= 0.85) and (abs(c1) < 1e-3 and abs(c2) < 1e-3)

    if is_asym_spin or is_non_spin:
        valid_rows.append({
            "sxs_id": sxs_id,
            "q": round(q, 4),
            "chi1_z": round(c1, 4),
            "chi2_z": round(c2, 4),
            "ecc": row["reference_eccentricity"],
            "regimen": "spinning_asym" if is_asym_spin else "non_spinning",
        })

candidatas = pd.DataFrame(valid_rows)
print(f"Simulaciones 100% circulares y sin precesión encontradas: {len(candidatas)}")
with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 1000):
    print(candidatas.to_string(index=False))