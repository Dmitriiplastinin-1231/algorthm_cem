"""
Numerical simulation of rocket trajectory displacement during sea launch.

This script implements Method B (numerical integration) from the analysis of
launch-platform tilt effects on rocket trajectory. The approach follows the
methodology described in:

  Greensite A.L., "Analysis and Design of Space Vehicle Flight Control Systems",
  NASA CR-820, 1970. Available open access at:
  https://ntrs.nasa.gov/citations/19700015859

The equations of motion are based on:
  Лебедев А.А., Чернобровкин Л.С., "Динамика полёта беспилотных летательных
  аппаратов", Машиностроение, 1973 — translational and rotational dynamics.
  Колесников К.С. (ред.), "Динамика ракет", Машиностроение, 2003 — initial
  conditions for launch from a moving/tilted platform.

TVC parameters and control-loop delay model are treated analogously to:
  Frosch J.A., Vallely D.P., "Saturn AS-501/S-IC Flight Control System Design",
  Journal of Spacecraft and Rockets, Vol.4 No.8, pp.1003-1009, 1967.
The 0.3 s delay is specified in the problem statement (task condition).

Usage:
    python rocket_simulation.py

Outputs:
    - Console summary of lateral displacement at 10 km altitude
    - Plots saved to rocket_trajectory.png
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Physical & mission parameters
# ---------------------------------------------------------------------------

# --- Rocket parameters (representative of a medium-lift vehicle, ~Zenit class) ---
M0 = 450_000.0        # Initial mass, kg
MDOT = 1_600.0        # Propellant mass-flow rate, kg/s (constant-thrust assumption)
THRUST = 7_257_600.0  # Vacuum thrust, N  (specific impulse ~460 s, sea-level ~390 s)
# Sea-level effective thrust (accounting for nozzle back-pressure at launch)
THRUST_SL = 6_300_000.0  # N
L_ROCKET = 57.0       # Total rocket length, m
L_CG = 25.0           # Distance from nozzle exit plane to initial CG, m
L_CP = 30.0           # Distance from nozzle exit plane to centre of pressure, m
DIAMETER = 3.9        # Rocket body diameter, m
REF_AREA = np.pi * (DIAMETER / 2) ** 2  # Reference area, m²

# Moment of inertia about lateral axis (Izz) – simplified as slender cylinder
# I_z ≈ m * L² / 12  (uniform rod approximation)
IZ0 = M0 * L_ROCKET**2 / 12.0  # Initial moment of inertia, kg·m²

# --- Atmosphere (simple exponential model) ---
RHO0 = 1.225          # Sea-level air density, kg/m³
H_SCALE = 8_500.0     # Scale height, m
CD_BODY = 0.3         # Axial drag coefficient
CN_ALPHA = 2.0        # Normal-force slope, 1/rad (fin + body)

# --- Control system ---
TAU_DELAY = 0.3       # Measurement-to-actuation delay, s
KP = 3.0              # Proportional gain (attitude angle), rad/rad
KD = 1.5              # Derivative gain (attitude rate), s/rad
DELTA_MAX = np.radians(6.0)    # Maximum TVC deflection, rad
DELTA_DOT_MAX = np.radians(5.0)  # Maximum TVC slew rate, rad/s

# --- Sea-launch platform rocking ---
PLATFORM_AMP = np.radians(2.0)   # Amplitude of platform tilt, rad (2 degrees)
PLATFORM_PERIOD = 8.0            # Period of rocking, s
PLATFORM_OMEGA = 2 * np.pi / PLATFORM_PERIOD  # Angular frequency, rad/s

# --- Simulation ---
G = 9.80665           # Standard gravity, m/s²
T_END = 120.0         # Simulation end time, s (well past 10 km altitude)
DT_HISTORY = 0.01     # History resolution for delay buffer, s


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def mass(t: float) -> float:
    """Rocket mass as a function of time (linear burn-down)."""
    return max(M0 - MDOT * t, M0 * 0.15)  # retain ~15 % structural mass


def thrust(t: float) -> float:
    """Effective thrust (simple altitude-dependent correction)."""
    # Use sea-level value for first ~5 s, then transition to vacuum value
    blend = min(t / 5.0, 1.0)
    return THRUST_SL + blend * (THRUST - THRUST_SL)


def moment_of_inertia(t: float) -> float:
    """Approximate moment of inertia (shrinks as propellant burns)."""
    m = mass(t)
    return m * L_ROCKET**2 / 12.0


def atmosphere(y: float) -> float:
    """Return air density at altitude y [m] using an exponential model."""
    rho = RHO0 * np.exp(-y / H_SCALE)
    return rho


def platform_tilt(t: float, phase: float = 0.0) -> tuple[float, float]:
    """
    Platform tilt angle and angular velocity at time t.

    Parameters
    ----------
    t     : float  time, s
    phase : float  initial phase of rocking cycle, rad

    Returns
    -------
    alpha_p   : float  tilt angle, rad
    omega_p   : float  angular velocity, rad/s
    """
    alpha_p = PLATFORM_AMP * np.sin(PLATFORM_OMEGA * t + phase)
    omega_p = PLATFORM_AMP * PLATFORM_OMEGA * np.cos(PLATFORM_OMEGA * t + phase)
    return alpha_p, omega_p


# ---------------------------------------------------------------------------
# Delay buffer (circular-history interpolation)
# ---------------------------------------------------------------------------

class DelayBuffer:
    """
    Stores the history of a signal and returns its value TAU_DELAY seconds ago.
    Uses linear interpolation between stored samples.

    Implementation: a growing list is used (max ~TAU/DT + 10 entries) – simple
    and correct, without circular-buffer indexing complexity.
    """

    def __init__(self, tau: float = TAU_DELAY, dt: float = DT_HISTORY):
        self.tau = tau
        self.dt = dt
        self._times: list[float] = []
        self._vals: list[float] = []

    def push(self, t: float, val: float) -> None:
        self._times.append(t)
        self._vals.append(val)
        # Trim entries older than (tau + 2*dt) – we no longer need them
        cutoff = t - self.tau - 2 * self.dt
        while len(self._times) > 2 and self._times[0] < cutoff:
            self._times.pop(0)
            self._vals.pop(0)

    def delayed(self, t: float) -> float:
        """Return the value at time (t - tau), interpolating stored samples."""
        t_query = t - self.tau
        if t_query <= 0.0 or len(self._times) < 2:
            return 0.0
        times = self._times
        vals = self._vals
        # Binary search for the last index where times[i] <= t_query
        lo, hi = 0, len(times) - 1
        if times[lo] > t_query:
            return vals[lo]
        if times[hi] <= t_query:
            return vals[hi]
        while lo + 1 < hi:
            mid = (lo + hi) // 2
            if times[mid] <= t_query:
                lo = mid
            else:
                hi = mid
        t0, v0 = times[lo], vals[lo]
        t1, v1 = times[hi], vals[hi]
        if abs(t1 - t0) < 1e-12:
            return v0
        frac = np.clip((t_query - t0) / (t1 - t0), 0.0, 1.0)
        return v0 + frac * (v1 - v0)


# ---------------------------------------------------------------------------
# ODE right-hand side
# ---------------------------------------------------------------------------

def build_ode(phase: float = 0.0) -> callable:
    """
    Build and return the ODE function for a given platform-rocking phase.

    State vector:
        x[0] = x_pos     horizontal displacement, m   (lateral, positive downrange)
        x[1] = y_pos     vertical displacement, m     (altitude, positive up)
        x[2] = vx        horizontal velocity, m/s
        x[3] = vy        vertical velocity, m/s
        x[4] = theta     pitch angle (tilt from vertical), rad  (positive = nose away from vertical)
        x[5] = omega_z   pitch angular velocity, rad/s

    Control uses a PD-controller with TAU_DELAY second delay on both theta and omega_z.
    """
    buf_theta = DelayBuffer(tau=TAU_DELAY, dt=DT_HISTORY)
    buf_omega = DelayBuffer(tau=TAU_DELAY, dt=DT_HISTORY)
    buf_delta = DelayBuffer(tau=TAU_DELAY, dt=DT_HISTORY)

    # Mutable state shared across calls (for delay buffer management)
    state_cache = {"t_last_push": -1.0, "delta_prev": 0.0}

    def ode(t: float, state: np.ndarray) -> np.ndarray:
        x_pos, y_pos, vx, vy, theta, omega_z = state

        m = mass(t)
        P = thrust(t)
        Iz = moment_of_inertia(t)
        rho = atmosphere(y_pos)

        # --- Push current measurements into delay buffer ---
        # Push at most once per DT_HISTORY interval (0.9 factor avoids floating-
        # point edge cases where t advances by just under DT_HISTORY exactly).
        if t - state_cache["t_last_push"] >= DT_HISTORY * 0.9:
            buf_theta.push(t, theta)
            buf_omega.push(t, omega_z)
            state_cache["t_last_push"] = t

        # --- Compute commanded deflection (ideal controller without delay) ---
        # Programme angle: nominal trajectory is vertical (theta_prog = 0) for first 10 km
        theta_meas_delayed = buf_theta.delayed(t)
        omega_meas_delayed = buf_omega.delayed(t)
        delta_cmd = -(KP * theta_meas_delayed + KD * omega_meas_delayed)
        delta_cmd = np.clip(delta_cmd, -DELTA_MAX, DELTA_MAX)

        # Rate-limit the deflection.
        # dt_actual is the delay-buffer resolution (DT_HISTORY), which is also
        # the effective sampling period for the discrete control loop. Using a
        # fixed dt here is consistent with the discrete-time controller model
        # described in Greensite (NASA CR-820, section 9.2).
        delta_prev = state_cache["delta_prev"]
        dt_actual = DT_HISTORY
        delta_dot = (delta_cmd - delta_prev) / dt_actual
        if abs(delta_dot) > DELTA_DOT_MAX:
            delta_dot = np.sign(delta_dot) * DELTA_DOT_MAX
        delta = delta_prev + delta_dot * dt_actual
        delta = np.clip(delta, -DELTA_MAX, DELTA_MAX)
        state_cache["delta_prev"] = delta

        # --- Velocity magnitude and angle of attack ---
        V = np.sqrt(vx**2 + vy**2) + 1e-6  # avoid division by zero
        # Angle of attack (small-angle, 2D)
        gamma = np.arctan2(vx, vy)  # flight-path angle from vertical
        alpha_aero = theta - gamma  # angle of attack

        # --- Aerodynamic forces (in body axes, projected to inertial) ---
        q_dyn = 0.5 * rho * V**2
        F_drag = q_dyn * REF_AREA * CD_BODY  # axial drag
        F_normal = q_dyn * REF_AREA * CN_ALPHA * alpha_aero  # normal (side) force

        # Drag is aligned with velocity; normal force is perpendicular to velocity
        # Project to inertial frame
        sin_gamma = np.sin(gamma)
        cos_gamma = np.cos(gamma)
        Fx_aero = -F_drag * sin_gamma - F_normal * cos_gamma
        Fy_aero = -F_drag * cos_gamma + F_normal * sin_gamma

        # --- Thrust force components ---
        # Thrust acts along rocket axis rotated by TVC deflection delta
        # Positive theta: nose tilted in +x direction
        Fx_thrust = P * np.sin(theta + delta)
        Fy_thrust = P * np.cos(theta + delta)

        # --- Translational equations of motion ---
        ax = (Fx_thrust + Fx_aero) / m
        ay = (Fy_thrust + Fy_aero) / m - G

        # --- Rotational equation of motion ---
        # Control moment: thrust through TVC arm
        arm_ctrl = L_CG       # CG-to-nozzle distance
        M_ctrl = P * arm_ctrl * np.sin(delta)

        # Aerodynamic restoring/destabilising moment about CG
        arm_aero = L_CP - L_CG   # positive = CP behind CG (stable configuration)
        M_aero = -F_normal * arm_aero

        alpha_z_dot = (M_ctrl + M_aero) / Iz

        return [vx, vy, ax, ay, omega_z, alpha_z_dot]

    return ode


# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------

def initial_conditions(phase: float = 0.0) -> np.ndarray:
    """
    Return initial state vector.

    Per Колесников (2003): at liftoff the rocket inherits the platform's tilt
    angle and angular velocity.

    Parameters
    ----------
    phase : rocking phase at t=0, rad

    Returns
    -------
    state0 : ndarray of shape (6,)
    """
    alpha_p0, omega_p0 = platform_tilt(0.0, phase)
    # x, y, vx, vy, theta, omega_z
    return np.array([0.0, 0.0, 0.0, 0.0, alpha_p0, omega_p0])


# ---------------------------------------------------------------------------
# Run simulation
# ---------------------------------------------------------------------------

def run_simulation(phase: float = 0.0, t_end: float = T_END) -> dict:
    """
    Integrate equations of motion from t=0 until t_end (or until y > 15 km).

    Parameters
    ----------
    phase  : float  initial rocking phase, rad
    t_end  : float  maximum simulation time, s

    Returns
    -------
    dict with keys: t, x, y, vx, vy, theta, omega_z
    """
    ode_func = build_ode(phase=phase)
    state0 = initial_conditions(phase=phase)

    def event_12km(t, state):
        return state[1] - 12_000.0  # stop when y > 12 km

    event_12km.terminal = True
    event_12km.direction = 1

    def event_10km(t, state):
        return state[1] - 10_000.0  # record 10 km crossing

    event_10km.terminal = False
    event_10km.direction = 1

    sol = solve_ivp(
        ode_func,
        [0.0, t_end],
        state0,
        method="RK45",
        max_step=0.05,      # 50 ms max step — per Frosch & Vallely methodology
        rtol=1e-6,
        atol=1e-8,
        dense_output=False,
        events=[event_10km, event_12km],
    )

    return {
        "t": sol.t,
        "x": sol.y[0],
        "y": sol.y[1],
        "vx": sol.y[2],
        "vy": sol.y[3],
        "theta": sol.y[4],
        "omega_z": sol.y[5],
    }


def lateral_displacement_at_10km(result: dict) -> float:
    """
    Interpolate lateral (horizontal) displacement when altitude == 10 km.
    """
    y_arr = result["y"]
    x_arr = result["x"]
    idx = np.searchsorted(y_arr, 10_000.0)
    if idx == 0 or idx >= len(y_arr):
        return float("nan")
    y0, y1 = y_arr[idx - 1], y_arr[idx]
    x0, x1 = x_arr[idx - 1], x_arr[idx]
    frac = (10_000.0 - y0) / (y1 - y0)
    return x0 + frac * (x1 - x0)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 65)
    print("Sea-launch rocket trajectory simulation")
    print("Numerical integration: RK45 (scipy.integrate.solve_ivp)")
    print("=" * 65)

    # --- Case 1: Worst-case — platform at maximum angular velocity (phase = 0) ---
    # At phase=0, alpha_p(0)=0 but omega_p(0)=max → maximum angular rate imparted
    print("\n[Case 1] Phase=0  (max angular velocity at launch)")
    res_worst = run_simulation(phase=0.0)
    dx_worst = lateral_displacement_at_10km(res_worst)
    print(f"  Lateral displacement at 10 km: {dx_worst:+.2f} m")

    # --- Case 2: Maximum initial tilt (phase = π/2) ---
    print("\n[Case 2] Phase=π/2 (max tilt angle at launch)")
    res_tilt = run_simulation(phase=np.pi / 2)
    dx_tilt = lateral_displacement_at_10km(res_tilt)
    print(f"  Lateral displacement at 10 km: {dx_tilt:+.2f} m")

    # --- Case 3: Phase sweep to find worst-case ---
    print("\n[Case 3] Phase sweep (0 to 2π, 36 points)")
    phases = np.linspace(0, 2 * np.pi, 36, endpoint=False)
    displacements = []
    for ph in phases:
        r = run_simulation(phase=ph)
        displacements.append(lateral_displacement_at_10km(r))
    displacements = np.array(displacements)
    max_disp = np.nanmax(np.abs(displacements))
    print(f"  Maximum |displacement| across all phases: {max_disp:.2f} m")

    print("\n" + "-" * 65)
    tol = 1000.0  # 1 km
    if max_disp < tol:
        print(f"  CONCLUSION: displacement ({max_disp:.1f} m) < ±{tol:.0f} m → negligible.")
        print("  The trajectory offset can be NEGLECTED for satellite deployment.")
    else:
        print(f"  CONCLUSION: displacement ({max_disp:.1f} m) ≥ ±{tol:.0f} m → significant.")
        print("  The trajectory offset CANNOT be neglected for satellite deployment.")
    print("-" * 65)

    # --- Plots ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(
        "Sea-Launch Rocket Trajectory Simulation\n"
        "(Platform rocking: A=2°, T=8 s; Control delay τ=0.3 s)",
        fontsize=13,
    )

    # Trajectory plot
    ax = axes[0, 0]
    ax.plot(res_worst["x"], res_worst["y"] / 1000, label="Phase=0 (max ω)", color="tab:blue")
    ax.plot(res_tilt["x"], res_tilt["y"] / 1000, label="Phase=π/2 (max θ)", color="tab:orange")
    ax.axhline(10, color="grey", linestyle="--", linewidth=0.8, label="10 km")
    ax.set_xlabel("Lateral displacement x, m")
    ax.set_ylabel("Altitude y, km")
    ax.set_title("Trajectory in the launch plane")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.4)

    # Pitch angle vs time
    ax = axes[0, 1]
    ax.plot(res_worst["t"], np.degrees(res_worst["theta"]), label="Phase=0", color="tab:blue")
    ax.plot(res_tilt["t"], np.degrees(res_tilt["theta"]), label="Phase=π/2", color="tab:orange")
    ax.set_xlabel("Time, s")
    ax.set_ylabel("Pitch angle θ, deg")
    ax.set_title("Pitch angle vs time")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.4)

    # Phase sweep: displacement at 10 km
    ax = axes[1, 0]
    ax.plot(np.degrees(phases), displacements, "o-", markersize=4, color="tab:green")
    ax.axhline(0, color="k", linewidth=0.5)
    ax.axhline(1000, color="red", linestyle="--", linewidth=0.8, label="±1 km")
    ax.axhline(-1000, color="red", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Launch phase φ₀, deg")
    ax.set_ylabel("Lateral displacement at 10 km, m")
    ax.set_title("Displacement vs launch phase")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.4)

    # Altitude vs time
    ax = axes[1, 1]
    ax.plot(res_worst["t"], res_worst["y"] / 1000, label="Phase=0", color="tab:blue")
    ax.plot(res_tilt["t"], res_tilt["y"] / 1000, label="Phase=π/2", color="tab:orange")
    ax.axhline(10, color="grey", linestyle="--", linewidth=0.8, label="10 km")
    ax.set_xlabel("Time, s")
    ax.set_ylabel("Altitude, km")
    ax.set_title("Altitude vs time")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig("rocket_trajectory.png", dpi=150)
    print("\nPlots saved to rocket_trajectory.png")

    # Numerical summary table
    print("\nPhase sweep summary (every 5th point):")
    print(f"  {'Phase (deg)':>12}  {'Displacement (m)':>18}")
    print(f"  {'-'*12}  {'-'*18}")
    for i, (ph, dp) in enumerate(zip(phases, displacements)):
        if i % 5 == 0:
            print(f"  {np.degrees(ph):>12.1f}  {dp:>18.2f}")


if __name__ == "__main__":
    main()
