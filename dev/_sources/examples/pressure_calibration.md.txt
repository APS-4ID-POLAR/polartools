# Pressure calibration

The `polartools.pressure_calibration` module computes pressure from
diamond-anvil-cell (DAC) experiments using the equation-of-state of an
internal calibrant — Au, Ag, or Pt. You measure a Bragg peak from the
calibrant; this module fits its position, converts it to a unit-cell
volume, and inverts the calibrant's EOS to give pressure.

## End-to-end calibration in one call

```python
from polartools.pressure_calibration import xrd_calibrate_pressure

result = xrd_calibrate_pressure(
    scan_id=42,
    source="myfile.dat",
    folder="/data/2024-1",
    x="tth", y="apd", monitor="I0",
    calibrant="Au",          # "Au", "Ag", or "Pt"
    energy_keV=20.0,
    temperature=300,         # K
    hkl=(1, 1, 1),           # which Bragg reflection
    model="PseudoVoigt",
)

print(f"Pressure: {result['P_GPa']:.2f} GPa")
print(f"Volume:   {result['V_A3']:.3f} Å³")
print(f"2θ:       {result['tth_deg']:.4f} °")
```

This is the high-level helper most users want. Under the hood it:

1. Loads scan `42` from `myfile.dat`
2. Fits the (111) Bragg peak with the chosen lmfit model
3. Converts peak 2θ → d-spacing → unit-cell volume (using the calibrant's
   crystal symmetry)
4. Inverts the calibrant's third-order Birch–Murnaghan EOS (Au/Ag) or the
   Holzapfel EOS (Pt) at the given temperature
5. Returns a dict of derived quantities

## EOS parameters at a given temperature

```python
from polartools.pressure_calibration import (
    load_ag_params, load_au_params, load_pt_params,
)

params_au = load_au_params(temperature=300)
# {'V0': 67.847, 'K0': 167.0, 'Kp0': 5.5, ...}

params_ag = load_ag_params(temperature=295)
params_pt = load_pt_params()  # T-independent reference values
```

## Manual fit + pressure calculation

If you want more control over the fit, do the steps yourself:

```python
from polartools.pressure_calibration import (
    fit_bragg_peak, calculate_tth, calculate_pressure, load_au_params,
)

# 1. Fit the peak
peak_fit = fit_bragg_peak(
    scan_id=42,
    source="myfile.dat",
    folder="/data/2024-1",
    x="tth", y="apd", monitor="I0",
    model="Gaussian",
)
tth = peak_fit.params["center"].value

# 2. Convert to volume (cubic FCC: a from d-spacing of (111))
import numpy as np
d = wavelength_A / (2 * np.sin(np.deg2rad(tth) / 2))
a = d * np.sqrt(1**2 + 1**2 + 1**2)
V = a**3

# 3. Invert EOS for pressure
params = load_au_params(temperature=300)
P_GPa = calculate_pressure(V, calibrant="Au", **params)
```

## Going the other way (predicting 2θ at a known P)

If you know the pressure and want to predict where the calibrant peak should
appear at a given energy:

```python
from polartools.pressure_calibration import calculate_tth

tth = calculate_tth(
    pressure=10.0,
    calibrant="Au",
    energy_keV=20.0,
    temperature=300,
    hkl=(1, 1, 1),
)
print(f"Au (111) at 10 GPa, 20 keV, 300 K: 2θ = {tth:.3f} °")
```

Useful for choosing the detector position before a run.

## Calibrant choice

| Calibrant | EOS                      | Useful range          |
|-----------|--------------------------|------------------------|
| Au        | Birch–Murnaghan (3rd ord.)| 0–200 GPa, 0–2000 K   |
| Ag        | Birch–Murnaghan (3rd ord.)| 0–100 GPa             |
| Pt        | Holzapfel                 | 0–550 GPa             |

For room-temperature DAC work below ~50 GPa, Au is the most common choice.

## See also

- [Diffraction](diffraction.md) — peak fitting machinery used here
