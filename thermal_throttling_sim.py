"""

THERMAL THROTTLING SIMULATION — Si-C vs. Graphite Li-ion Battery Anodes

Model:      Newman–Tiedemann–Gu–Kim (NTGK) multiscale thermal framework
Numerics:   Explicit Finite Difference Method (EFDM), 1D transverse slab
Protocol:   100 W fast-charge (≈4C–5C) on a 5000 mAh dual-cell pouch battery
Geometry:   100 mm × 60 mm × 8 mm pouch cell, z-axis only (transverse)
Author:     Harshith Bhardwaz| NITW 

"""
# Importing libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# Domain and geometry parameters
L  = 0.008          # slab thickness [m] — 8 mm pouch-cell transverse dimension
N  = 20             # spatial node count — enforced per project spec (not adjustable)
dz = L / (N - 1)   # uniform node spacing [m]

V_cell = 100e-3 * 60e-3 * 8e-3   # bulk cell volume [m³] = 4.8 × 10⁻⁵ m³
                                   # used to derive volumetric heat generation from watt-level losses

T_amb  = 25.0       # ambient temperature [°C] — standard JEDEC/IEC 62133 room-temp test condition
t_max  = 600.0      # hard simulation ceiling [s] — extended to 600 s to allow stabilization observation


# Material properties

#  All values are effective transverse (z-direction) bulk properties for a
#  calendared composite electrode laminate (active material + PVDF binder +
#  liquid electrolyte in pores).  In-plane values are ~10× higher but are
#  irrelevant for the 1D transverse model.

# --- Silicon-Carbon (Si-C) composite anode ---
rho_sic = 2100.0    # density [kg/m³]  — intermediate between pure Si (2330) and graphite (2200)
                    #                   lower than graphite due to carbon matrix dilution
Cp_sic  =  900.0    # specific heat [J/kg·K] — higher Cₚ → greater thermal inertia → slower ΔT rise
k_sic   =    4.5    # eff. transverse conductivity [W/m·K] — attenuated by PVDF binder & liquid electrolyte

# --- Graphite anode ---
rho_gr  = 2200.0    # density [kg/m³]  — calendared graphite is denser than Si-C composite
Cp_gr   =  700.0    # specific heat [J/kg·K] — LOWER than Si-C → heats up faster per joule deposited
k_gr    =    5.0    # eff. transverse conductivity [W/m·K] — slightly higher than Si-C (more crystalline)


# Heat generation parameters

# Total volumetric heat rate aggregates three NTGK terms:
#   q_ohm  = I²·R_int  — Joule heating through tabs & current collectors
#   q_irr  = I(U−V)    — irreversible polarization losses (SEI, pore diffusion)
#   q_rev  = I·T·∂U/∂T — entropic reversible heating from lithiation reactions
#
# Derivation for q_nominal:
#   I = 10 A (4C–5C on 2500 mAh/cell × 2 cells = 5000 mAh pack)
#   R_int ≈ 0.05 Ω (internal impedance at elevated SOC, measured via EIS)
#   q_ohm = I²·R = 100 × 0.05 = 5 W
#   Polarization + entropic scaling at 4C–5C adds ~40% → ≈ 7 W thermal load
#   q = 7 W / 4.8×10⁻⁵ m³ ≈ 145,833 W/m³ → normalized to 100,000 W/m³
#   (aligns with established CFD literature baselines for 4C–5C charge rates)

q_nominal   = 100_000.0   # volumetric heat rate at full 100 W charge [W/m³]
q_throttled =  27_000.0   # firmware step-down to ~27 W — heat gen falls BELOW passive dissipation rate
                           # at h = 2.5 W/m²·K, net dT/dt becomes negative → cooling begins


# Convective heat transfer coefficients (h) for boundary conditions — 3 representative cooling scenarios

h_cases = [
    (2.5,  "Low  — near-sealed smartphone chassis, restricted natural convection"),
    (6.0,  "Mid  — open air, partial spreading via internal graphene sheets"),
    (12.0, "High — light forced airflow or vapor chamber-assisted passive spreading"),
]


# Throttling parameters - Sch,itt trigger thresholds for BMS control logic

THROTTLE_ON  = 45.0   # [°C] upper trip point — Android THERMAL_STATUS_EMERGENCY / SEI degradation limit
THROTTLE_OFF = 43.0   # [°C] lower trip point — mandatory 2°C dead band prevents chatter oscillation

core = N // 2         # index of the mid-plane node (z = L/2 = 4 mm)
                      # T[core] = T_core = point of maximum thermal stagnation. This is the hottest point because heat generated everywhere in the slab has to escape to both surfaces.
                      # BMS infers T_core from surface telemetry + current integration;
                      # surface thermistors alone CANNOT measure it directly


# Criteria for numerical stability of the explicit finite difference method (EFDM):

alpha_sic      = k_sic / (rho_sic * Cp_sic)    # thermal diffusivity Si-C  [m²/s]
alpha_graphite = k_gr  / (rho_gr  * Cp_gr )    # thermal diffusivity Graphite [m²/s]
alpha_max      = max(alpha_sic, alpha_graphite) # use the FASTEST diffuser to set dt
                                                 # the fastest diffuser demands the smallest dt;
                                                 # all slower materials are automatically stable.
                                                 # (Using alpha_min would under-constrain dt and
                                                 #  let the faster material's Fo exceed 0.5.)

Fo = 0.4           # Fourier number — 20% safety buffer below the 0.5 stability limit
dt = Fo * dz**2 / alpha_max  # [s] — both materials share this dt for time-aligned x-axes 
                             # von Neumann stability criterion  

# Verify stability for both materials (assertion will raise if violated)
Fo_sic = alpha_sic      * dt / dz**2
Fo_gr  = alpha_graphite * dt / dz**2
assert Fo_sic <= 0.5, f"Stability violated for Si-C: Fo = {Fo_sic:.4f}"
assert Fo_gr  <= 0.5, f"Stability violated for Graphite: Fo = {Fo_gr:.4f}"

print(f"{'='*70}")
print(f"EFDM THERMAL THROTTLING SIMULATION — NTGK Framework")
print(f"{'='*70}")
print(f"  Spatial nodes     N  = {N}")
print(f"  Node spacing      dz = {dz*1000:.4f} mm")
print(f"  Timestep          dt = {dt:.6f} s")
print(f"  Fo (Si-C)            = {Fo_sic:.4f}  ✓ stable (≤ 0.5)")
print(f"  Fo (Graphite)        = {Fo_gr:.4f}  ✓ stable (≤ 0.5)")
print(f"  alpha_Si-C           = {alpha_sic:.6e} m²/s")
print(f"  alpha_Graphite       = {alpha_graphite:.6e} m²/s")
print(f"  q_nominal            = {q_nominal:.0f} W/m³")
print(f"  q_throttled          = {q_throttled:.0f} W/m³")
print(f"  t_max                = {t_max:.0f} s")
print(f"{'='*70}\n")

# Sampling interval for plotting — retain ~1 data point per simulated second
# prevents Matplotlib from rendering 300,000+ points (dt is on the order of ms)
plot_interval = max(1, int(1.0 / dt))


# Simulation loops

# Structure: results[h_label][material] = dict with arrays and event logs

materials = {
    "Si-C":     {"rho": rho_sic, "Cp": Cp_sic, "k": k_sic,
                 "color": "red",  "label": "Si-C"},
    "Graphite": {"rho": rho_gr,  "Cp": Cp_gr,  "k": k_gr,
                 "color": "blue", "label": "Graphite"},
}

results = {}   # top-level storage for all runs

for h_val, h_desc in h_cases:
    h_key = f"h={h_val:.1f}"
    results[h_key] = {}

    for mat_name, mat in materials.items():
        rho  = mat["rho"]
        Cp   = mat["Cp"]
        k    = mat["k"]

        # Pre-compute EFDM coefficients for this material
        alpha = k / (rho * Cp)
        Fo_m  = alpha * dt / dz**2                # Fourier number for this material
        coeff_conv = 2.0 * h_val * dt / (rho * Cp * dz)  # boundary convection coefficient
        coeff_q    = q_nominal * dt / (rho * Cp)           # volumetric heating term (changes with throttle)

        # Initial conditions
        T = np.full(N, T_amb, dtype=np.float64)   # uniform 25°C everywhere at t = 0
        t = 0.0                                    # simulation clock [s]
        q = q_nominal                              # current heat generation rate [W/m³]
        throttled = False                          # BMS state flag

        #  Plot data containers 
        t_plot     = []
        T_core_plot = []

        # Event log containers 
        throttle_on_events  = []   # list of (t, T_core) tuples
        throttle_off_events = []   # list of (t, T_core) tuples
        perm_throttle       = False
        throttle_never      = True  # assume never triggered until proven otherwise

        # Stabilization tracking
        T_core_snapshot = T[core]  # T_core value 10 simulated seconds ago (rolling update)
        snapshot_timer  = 0.0      # accumulates simulated time since last snapshot

        n_steps = int(t_max / dt)

        for step in range(n_steps):

            T_core_now = T[core]   # mid-plane temperature this timestep

            
            # THROTTLE LOGIC — Schmitt Trigger (evaluated every timestep)
            if not throttled and T_core_now >= THROTTLE_ON:
                throttled = True
                throttle_never = False
                q = q_throttled
                throttle_on_events.append((t, T_core_now))
                print(f"  [{mat_name:<8} | h={h_val:4.1f}]  Throttle ON   at t={t:7.2f}s  |  T_core = {T_core_now:.2f}°C")

            elif throttled and T_core_now <= THROTTLE_OFF:
                throttled = False
                q = q_nominal
                throttle_off_events.append((t, T_core_now))
                print(f"  [{mat_name:<8} | h={h_val:4.1f}]  Throttle OFF  at t={t:7.2f}s  |  T_core = {T_core_now:.2f}°C")

            
            # STABILIZATION CHECK — every 10 simulated seconds
            snapshot_timer += dt
            if snapshot_timer >= 10.0:
                delta_T = abs(T_core_now - T_core_snapshot)

                if throttled and delta_T < 0.05:
                    # T_core is stuck between 43°C and 45°C while throttled
                    # Recovery condition (≤ 43°C) is physically unreachable under current cooling
                    perm_throttle = True
                    print(f"  [{mat_name:<8} | h={h_val:4.1f}]  PERMANENT THROTTLE — T_core stabilized in 43–45°C dead band at t={t:.2f}s")
                    t_plot.append(t)
                    T_core_plot.append(T_core_now)
                    break

                elif not throttled and delta_T < 0.05 and step > int(30.0 / dt):
                    # Thermal stabilization in open-loop (unthrottled) state
                    # Guard: skip the first 30 s to avoid false-positive on flat initial transient
                    print(f"  [{mat_name:<8} | h={h_val:4.1f}]  Thermally stabilized at t={t:.2f}s  |  T_core = {T_core_now:.2f}°C")
                    t_plot.append(t)
                    T_core_plot.append(T_core_now)
                    break

                T_core_snapshot = T_core_now
                snapshot_timer  = 0.0

        
            # EFDM TEMPERATURE UPDATE — vectorized NumPy for speed
            T_new = T.copy()
            coeff_q_now = q * dt / (rho * Cp)   # recomputing because q may have changed

            # Interior nodes  (i = 1 … N-2):  standard second-order central difference
            T_new[1:-1] = (T[1:-1]
                           + Fo_m * (T[2:] - 2.0 * T[1:-1] + T[:-2])
                           + coeff_q_now)

            # Left boundary node (z = 0):  convective BC with ghost-node elimination
            T_new[0] = (T[0]
                        + 2.0 * Fo_m * (T[1] - T[0])
                        - coeff_conv * (T[0] - T_amb)
                        + coeff_q_now)

            # Right boundary node (z = L):  mirror of left BC using T[N-2] as interior neighbor
            T_new[-1] = (T[-1]
                         + 2.0 * Fo_m * (T[-2] - T[-1])
                         - coeff_conv * (T[-1] - T_amb)
                         + coeff_q_now)

            T = T_new
            t += dt

            
            # SAMPLING — retain ~1 point per simulated second
            if step % plot_interval == 0:
                t_plot.append(t)
                T_core_plot.append(T[core])

        else:
            # Loop completed without break — t_max reached
            if throttle_never:
                print(f"  [{mat_name:<8} | h={h_val:4.1f}]  Throttle NEVER triggered      |  T_core_max = {max(T_core_plot):.2f}°C")
            else:
                print(f"  [{mat_name:<8} | h={h_val:4.1f}]  Simulation reached t_max={t_max:.0f}s  |  T_core_final = {T[core]:.2f}°C")

        # Store all run data for plotting
        results[h_key][mat_name] = {
            "t":                  np.array(t_plot),
            "T_core":             np.array(T_core_plot),
            "throttle_on":        throttle_on_events,
            "throttle_off":       throttle_off_events,
            "perm_throttle":      perm_throttle,
            "throttle_never":     throttle_never,
            "color":              mat["color"],
            "label":              mat["label"],
        }

    print()   # blank line between h-value groups


# 7.  PLOTTING
fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=150)
fig.suptitle(
    "Thermal Throttling — Si-C vs. Graphite  |  100 W Fast Charge  |  EFDM Simulation\n"
    "NTGK Framework · 1D Transverse Slab · N=20 nodes · Fo=0.4",
    fontsize=12, fontweight="bold", y=1.02
)

for ax, (h_val, _) in zip(axes, h_cases):
    h_key = f"h={h_val:.1f}"
    run   = results[h_key]

    for mat_name, data in run.items():
        color = data["color"]
        t_arr = data["t"]
        T_arr = data["T_core"]

        # Main temperature trace
        ax.plot(t_arr, T_arr,
                color=color, linewidth=1.8,
                label=data["label"])

        # Throttle ON events — dashed vertical line
        for (t_ev, T_ev) in data["throttle_on"]:
            ax.axvline(t_ev, color=color, linestyle="--", alpha=0.6, linewidth=1.0)
            ax.annotate(f"ON\n{t_ev:.0f}s",
                        xy=(t_ev, T_ev),
                        xytext=(t_ev + 3, T_ev + 0.4),
                        fontsize=6.5, color=color, alpha=0.85)

        # Throttle OFF events — dotted vertical line
        for (t_ev, T_ev) in data["throttle_off"]:
            ax.axvline(t_ev, color=color, linestyle=":", alpha=0.4, linewidth=1.0)
            ax.annotate(f"OFF\n{t_ev:.0f}s",
                        xy=(t_ev, T_ev),
                        xytext=(t_ev + 3, T_ev - 1.2),
                        fontsize=6.5, color=color, alpha=0.75)

    # Reference threshold lines
    ax.axhline(THROTTLE_ON,  color="black", linestyle="--", linewidth=0.8,
               label="Throttle threshold (45°C)")
    ax.axhline(THROTTLE_OFF, color="gray",  linestyle=":",  linewidth=0.8,
               label="Recovery threshold (43°C)")

    # Shaded dead band between 43°C and 45°C
    ax.axhspan(THROTTLE_OFF, THROTTLE_ON, color="lightyellow",
               alpha=0.45, label="2°C dead band")

    ax.set_xlabel("Time (s)", fontsize=10)
    ax.set_ylabel("Core Temperature (°C)", fontsize=10)
    ax.set_title(f"h = {h_val:.1f} W/m²·K", fontsize=11, fontweight="bold")
    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.85)
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.6)
    ax.set_xlim(left=0)

plt.tight_layout()
plt.savefig("thermal_throttle_simulation.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nPlot saved → thermal_throttle_simulation.png")


# 8.  SUMMARY TABLE
print(f"\n{'='*70}")
print(f"  SIMULATION SUMMARY")
print(f"{'='*70}")
print(f"  {'Material':<10} {'h (W/m²K)':<12} {'Throttle ON (s)':<18} {'T_core_max (°C)':<18} {'Status'}")
print(f"  {'-'*68}")
for h_val, _ in h_cases:
    h_key = f"h={h_val:.1f}"
    for mat_name, data in results[h_key].items():
        t_on   = data["throttle_on"][0][0]  if data["throttle_on"]  else "—"
        t_on_s = f"{t_on:.1f}" if isinstance(t_on, float) else t_on
        T_max  = float(np.max(data["T_core"]))
        status = ("Permanent throttle" if data["perm_throttle"] or
                                          (data["throttle_on"] and not data["throttle_off"] and not data["throttle_never"])
                  else "Never triggered" if data["throttle_never"]
                  else "Normal cycling")
        print(f"  {mat_name:<10} {h_val:<12.1f} {t_on_s:<18} {T_max:<18.2f} {status}")
print(f"{'='*70}\n")
