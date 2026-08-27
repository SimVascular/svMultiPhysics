# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Independent TP06 EPI, ENDO, and M model implementation."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
from math import exp, isfinite, log, sqrt
from typing import Dict, Iterable, Literal, Tuple


STIMULUS_AMPLITUDE = -52.0
STIMULUS_START_MS = 10.0
STIMULUS_END_MS = 11.0
Phenotype = Literal["epi", "endo", "m"]


@dataclass(frozen=True)
class Parameters:
    """TP06 parameters in the numerical convention used by svMP."""

    R: float = 8314.472
    T: float = 310.0
    F: float = 96485.3415
    Cm: float = 0.185
    V_c: float = 16.404e-3
    V_SR: float = 1.094e-3
    V_ss: float = 5.468e-5

    K_o: float = 5.4
    Na_o: float = 140.0
    Ca_o: float = 2.0
    P_KNa: float = 0.03

    G_Na: float = 14.838
    G_K1: float = 5.405
    G_to: float = 0.294
    G_Kr: float = 0.153
    G_Ks: float = 0.392
    G_CaL: float = 3.98e-5
    G_bNa: float = 2.9e-4
    G_bCa: float = 5.92e-4
    G_pK: float = 1.46e-2
    G_pCa: float = 0.1238
    K_pCa: float = 5.0e-4

    K_NaCa: float = 1000.0
    gamma: float = 0.35
    K_mCa: float = 1.38
    K_mNai: float = 87.5
    K_sat: float = 0.1
    alpha: float = 2.5

    P_NaK: float = 2.724
    K_mK: float = 1.0
    K_mNa: float = 40.0

    Vmax_up: float = 6.375e-3
    K_up: float = 2.5e-4
    V_rel: float = 0.102
    k1_prime: float = 0.15
    k2_prime: float = 4.5e-2
    k3: float = 6.0e-2
    k4: float = 5.0e-3
    EC: float = 1.5
    max_SR: float = 2.5
    min_SR: float = 1.0
    V_leak: float = 3.6e-4
    V_xfer: float = 3.8e-3

    Buf_c: float = 0.2
    K_buf_c: float = 1.0e-3
    Buf_SR: float = 10.0
    K_buf_SR: float = 0.3
    Buf_ss: float = 0.4
    K_buf_ss: float = 2.5e-4

    V_rest: float = -85.23


@dataclass(frozen=True)
class MainState:
    """The seven TP06 states stored in svMP's main X vector."""

    V: float
    Ki: float
    Nai: float
    Cai: float
    Ca_ss: float
    Ca_SR: float
    R_prime: float


@dataclass(frozen=True)
class Gates:
    """The twelve TP06 states advanced by Rush--Larsen in svMP."""

    Xr1: float
    Xr2: float
    Xs: float
    m: float
    h: float
    j: float
    d: float
    f: float
    f2: float
    fCass: float
    s: float
    r: float


@dataclass(frozen=True)
class ReversalPotentials:
    EK: float
    ENa: float
    ECa: float
    EKs: float


@dataclass(frozen=True)
class Currents:
    IK1: float
    Ito: float
    IKr: float
    IKs: float
    ICaL: float
    INaK: float
    INa: float
    IbNa: float
    INaCa: float
    IbCa: float
    IpK: float
    IpCa: float
    I_rel: float
    I_up: float
    I_leak: float
    I_xfer: float


@dataclass(frozen=True)
class StepResult:
    state: MainState
    gates: Gates


def default_epi_state() -> MainState:
    return MainState(
        V=-85.23,
        Ki=136.89,
        Nai=8.6040,
        Cai=1.26e-4,
        Ca_ss=3.6e-4,
        Ca_SR=3.64,
        R_prime=0.9073,
    )


def default_epi_gates() -> Gates:
    return Gates(
        Xr1=6.21e-3,
        Xr2=0.4712,
        Xs=9.5e-3,
        m=1.72e-3,
        h=0.7444,
        j=0.7045,
        d=3.373e-5,
        f=0.7888,
        f2=0.9755,
        fCass=0.9953,
        s=0.999998,
        r=2.42e-8,
    )


def endo_parameters() -> Parameters:
    """Return the TP06 ENDO conductance configuration."""

    return replace(Parameters(), G_to=0.073)


def default_endo_state() -> MainState:
    return MainState(
        V=-86.709,
        Ki=138.4,
        Nai=10.355,
        Cai=1.3e-4,
        Ca_ss=3.6e-4,
        Ca_SR=3.715,
        R_prime=0.9068,
    )


def default_endo_gates() -> Gates:
    return Gates(
        Xr1=4.48e-3,
        Xr2=0.476,
        Xs=8.7e-3,
        m=1.55e-3,
        h=0.7573,
        j=0.7225,
        d=3.164e-5,
        f=0.8009,
        f2=0.9778,
        fCass=0.9953,
        s=0.3212,
        r=2.235e-8,
    )


def m_parameters() -> Parameters:
    """Return the TP06 M-cell conductance configuration."""

    return replace(Parameters(), G_Ks=0.098)


def default_m_state() -> MainState:
    return MainState(
        V=-85.423,
        Ki=138.52,
        Nai=10.132,
        Cai=1.53e-4,
        Ca_ss=4.2e-4,
        Ca_SR=4.272,
        R_prime=0.8978,
    )


def default_m_gates() -> Gates:
    return Gates(
        Xr1=1.65e-2,
        Xr2=0.473,
        Xs=1.74e-2,
        m=1.65e-3,
        h=0.749,
        j=0.6788,
        d=3.288e-5,
        f=0.7026,
        f2=0.9526,
        fCass=0.9942,
        s=0.999998,
        r=2.347e-8,
    )


def stimulus_current(time_ms: float) -> float:
    """Return the canonical half-open TP06 stimulus at old-state time."""

    if not isfinite(time_ms):
        raise ValueError("stimulus time must be finite")
    if STIMULUS_START_MS <= time_ms < STIMULUS_END_MS:
        return STIMULUS_AMPLITUDE
    return 0.0


def reversal_potentials(
    state: MainState, parameters: Parameters = Parameters()
) -> ReversalPotentials:
    voltage_factor = parameters.R * parameters.T / parameters.F
    return ReversalPotentials(
        EK=voltage_factor * log(parameters.K_o / state.Ki),
        ENa=voltage_factor * log(parameters.Na_o / state.Nai),
        ECa=0.5 * voltage_factor * log(parameters.Ca_o / state.Cai),
        EKs=voltage_factor
        * log(
            (parameters.K_o + parameters.P_KNa * parameters.Na_o)
            / (state.Ki + parameters.P_KNa * state.Nai)
        ),
    )


def gate_kinetics(
    state: MainState, phenotype: Phenotype = "epi"
) -> Dict[str, Tuple[float, float]]:
    """Return ``gate_name: (steady_state, time_constant_ms)``."""

    V = state.V

    Xr1_inf = 1.0 / (1.0 + exp(-(26.0 + V) / 7.0))
    Xr1_tau = (
        450.0 / (1.0 + exp(-(45.0 + V) / 10.0))
        * 6.0
        / (1.0 + exp((30.0 + V) / 11.5))
    )

    Xr2_inf = 1.0 / (1.0 + exp((88.0 + V) / 24.0))
    Xr2_tau = (
        3.0 / (1.0 + exp(-(60.0 + V) / 20.0))
        * 1.12
        / (1.0 + exp(-(60.0 - V) / 20.0))
    )

    Xs_inf = 1.0 / (1.0 + exp(-(5.0 + V) / 14.0))
    Xs_tau = (
        1400.0 / sqrt(1.0 + exp((5.0 - V) / 6.0))
        / (1.0 + exp((V - 35.0) / 15.0))
        + 80.0
    )

    m_inf = 1.0 / (1.0 + exp(-(56.86 + V) / 9.03)) ** 2
    m_tau = (
        1.0 / (1.0 + exp(-(60.0 + V) / 5.0))
        * (
            0.1 / (1.0 + exp((35.0 + V) / 5.0))
            + 0.1 / (1.0 + exp((V - 50.0) / 200.0))
        )
    )

    h_inf = 1.0 / (1.0 + exp((71.55 + V) / 7.43)) ** 2
    if V >= -40.0:
        alpha_h = 0.0
        beta_h = 0.77 / (0.13 * (1.0 + exp(-(10.66 + V) / 11.1)))
    else:
        alpha_h = 0.057 * exp(-(80.0 + V) / 6.8)
        beta_h = 2.7 * exp(0.079 * V) + 310000.0 * exp(0.3485 * V)
    h_tau = 1.0 / (alpha_h + beta_h)

    j_inf = h_inf
    if V >= -40.0:
        alpha_j = 0.0
        beta_j = 0.6 * exp(0.057 * V) / (1.0 + exp(-0.1 * (V + 32.0)))
    else:
        alpha_j = -(
            25428.0 * exp(0.2444 * V) + 6.948e-6 * exp(-0.04391 * V)
        ) * (V + 37.78) / (1.0 + exp(0.311 * (79.23 + V)))
        beta_j = 0.02424 * exp(-0.01052 * V) / (
            1.0 + exp(-0.1378 * (40.14 + V))
        )
    j_tau = 1.0 / (alpha_j + beta_j)

    d_inf = 1.0 / (1.0 + exp(-(8.0 + V) / 7.5))
    d_tau = (
        (1.4 / (1.0 + exp(-(35.0 + V) / 13.0)) + 0.25)
        * 1.4
        / (1.0 + exp((5.0 + V) / 5.0))
        + 1.0 / (1.0 + exp((50.0 - V) / 20.0))
    )

    f_inf = 1.0 / (1.0 + exp((20.0 + V) / 7.0))
    f_tau = (
        1102.5 * exp(-((V + 27.0) ** 2) / 225.0)
        + 200.0 / (1.0 + exp((13.0 - V) / 10.0))
        + 180.0 / (1.0 + exp((30.0 + V) / 10.0))
        + 20.0
    )

    f2_inf = 0.67 / (1.0 + exp((35.0 + V) / 7.0)) + 0.33
    f2_tau = (
        562.0 * exp(-((27.0 + V) ** 2) / 240.0)
        + 31.0 / (1.0 + exp((25.0 - V) / 10.0))
        + 80.0 / (1.0 + exp((30.0 + V) / 10.0))
    )

    calcium_factor = 1.0 / (1.0 + (state.Ca_ss / 0.05) ** 2)
    fCass_inf = 0.6 * calcium_factor + 0.4
    fCass_tau = 80.0 * calcium_factor + 2.0

    if phenotype in ("epi", "m"):
        s_inf = 1.0 / (1.0 + exp((20.0 + V) / 5.0))
        s_tau = (
            85.0 * exp(-((V + 45.0) ** 2) / 320.0)
            + 5.0 / (1.0 + exp((V - 20.0) / 5.0))
            + 3.0
        )
    elif phenotype == "endo":
        s_inf = 1.0 / (1.0 + exp((28.0 + V) / 5.0))
        s_tau = 1000.0 * exp(-((V + 67.0) ** 2) / 1000.0) + 8.0
    else:
        raise ValueError(f"unsupported TP06 phenotype: {phenotype}")

    r_inf = 1.0 / (1.0 + exp((20.0 - V) / 6.0))
    r_tau = 9.5 * exp(-((V + 40.0) ** 2) / 1800.0) + 0.8

    return {
        "Xr1": (Xr1_inf, Xr1_tau),
        "Xr2": (Xr2_inf, Xr2_tau),
        "Xs": (Xs_inf, Xs_tau),
        "m": (m_inf, m_tau),
        "h": (h_inf, h_tau),
        "j": (j_inf, j_tau),
        "d": (d_inf, d_tau),
        "f": (f_inf, f_tau),
        "f2": (f2_inf, f2_tau),
        "fCass": (fCass_inf, fCass_tau),
        "s": (s_inf, s_tau),
        "r": (r_inf, r_tau),
    }


def ionic_and_calcium_currents(
    state: MainState,
    gates: Gates,
    parameters: Parameters = Parameters(),
) -> Currents:
    p = parameters
    reversal = reversal_potentials(state, p)
    V = state.V
    voltage_factor = p.R * p.T / p.F
    potassium_scale = sqrt(p.K_o / 5.4)

    alpha_K1 = 0.1 / (1.0 + exp(0.06 * (V - reversal.EK - 200.0)))
    beta_K1 = (
        3.0 * exp(0.0002 * (V - reversal.EK + 100.0))
        + exp(0.1 * (V - reversal.EK - 10.0))
    ) / (1.0 + exp(-0.5 * (V - reversal.EK)))
    xK1_inf = alpha_K1 / (alpha_K1 + beta_K1)

    IK1 = p.G_K1 * potassium_scale * xK1_inf * (V - reversal.EK)
    Ito = p.G_to * gates.r * gates.s * (V - reversal.EK)
    IKr = p.G_Kr * potassium_scale * gates.Xr1 * gates.Xr2 * (V - reversal.EK)
    IKs = p.G_Ks * gates.Xs**2 * (V - reversal.EKs)
    INa = p.G_Na * gates.m**3 * gates.h * gates.j * (V - reversal.ENa)
    IbNa = p.G_bNa * (V - reversal.ENa)
    IbCa = p.G_bCa * (V - reversal.ECa)
    IpK = p.G_pK * (V - reversal.EK) / (1.0 + exp((25.0 - V) / 5.98))
    IpCa = p.G_pCa * state.Cai / (p.K_pCa + state.Cai)

    calcium_exponent = 2.0 * (V - 15.0) / voltage_factor
    # This deliberately mirrors the published/svMP expression.  It has a
    # removable singularity at exactly V=15 mV; no alternate regularization is
    # introduced in this oracle.
    calcium_flux_factor = (
        2.0
        * calcium_exponent
        * p.F
        * (0.25 * state.Ca_ss * exp(calcium_exponent) - p.Ca_o)
        / (exp(calcium_exponent) - 1.0)
    )
    ICaL = (
        p.G_CaL
        * gates.d
        * gates.f
        * gates.f2
        * gates.fCass
        * calcium_flux_factor
    )

    exchanger_forward = exp(p.gamma * V / voltage_factor)
    exchanger_reverse = exp((p.gamma - 1.0) * V / voltage_factor)
    INaCa = p.K_NaCa * (
        exchanger_forward * state.Nai**3 * p.Ca_o
        - exchanger_reverse * p.Na_o**3 * state.Cai * p.alpha
    ) / (
        (p.K_mNai**3 + p.Na_o**3)
        * (p.K_mCa + p.Ca_o)
        * (1.0 + p.K_sat * exchanger_reverse)
    )

    INaK = p.P_NaK * p.K_o * state.Nai / (
        (p.K_o + p.K_mK)
        * (state.Nai + p.K_mNa)
        * (
            1.0
            + 0.1245 * exp(-0.1 * V / voltage_factor)
            + 0.0353 * exp(-V / voltage_factor)
        )
    )

    kcasr = p.max_SR - (p.max_SR - p.min_SR) / (
        1.0 + (p.EC / state.Ca_SR) ** 2
    )
    k1 = p.k1_prime / kcasr
    open_probability = (
        k1 * state.R_prime * state.Ca_ss**2
        / (p.k3 + k1 * state.Ca_ss**2)
    )

    I_rel = p.V_rel * open_probability * (state.Ca_SR - state.Ca_ss)
    I_up = p.Vmax_up / (1.0 + (p.K_up / state.Cai) ** 2)
    I_leak = p.V_leak * (state.Ca_SR - state.Cai)
    I_xfer = p.V_xfer * (state.Ca_ss - state.Cai)

    return Currents(
        IK1=IK1,
        Ito=Ito,
        IKr=IKr,
        IKs=IKs,
        ICaL=ICaL,
        INaK=INaK,
        INa=INa,
        IbNa=IbNa,
        INaCa=INaCa,
        IbCa=IbCa,
        IpK=IpK,
        IpCa=IpCa,
        I_rel=I_rel,
        I_up=I_up,
        I_leak=I_leak,
        I_xfer=I_xfer,
    )


def main_state_rhs(
    state: MainState,
    gates: Gates,
    parameters: Parameters = Parameters(),
    stimulus_current: float = 0.0,
    sac_current: float = 0.0,
) -> MainState:
    """Evaluate the seven main-state derivatives at one old-state point."""

    p = parameters
    current = ionic_and_calcium_currents(state, gates, p)

    dV = -(
        current.INa
        + current.Ito
        + current.IK1
        + current.IKr
        + current.IKs
        + current.ICaL
        + current.INaCa
        + current.INaK
        + current.IpCa
        + current.IpK
        + current.IbCa
        + current.IbNa
        + stimulus_current
    ) + sac_current

    concentration_scale = p.Cm / (p.V_c * p.F)
    dKi = -concentration_scale * (
        current.IK1
        + current.Ito
        + current.IKr
        + current.IKs
        + current.IpK
        - 2.0 * current.INaK
        + stimulus_current
    )
    dNai = -concentration_scale * (
        current.INa + current.IbNa + 3.0 * (current.INaK + current.INaCa)
    )

    cytosolic_total = (
        (current.I_leak - current.I_up) * p.V_SR / p.V_c
        + current.I_xfer
        - concentration_scale
        * (current.IbCa + current.IpCa - 2.0 * current.INaCa)
        / 2.0
    )
    cytosolic_buffer = 1.0 + p.K_buf_c * p.Buf_c / (
        state.Cai + p.K_buf_c
    ) ** 2
    dCai = cytosolic_total / cytosolic_buffer

    subspace_total = (
        -current.ICaL * p.Cm / (2.0 * p.F)
        + current.I_rel * p.V_SR
        - p.V_c * current.I_xfer
    ) / p.V_ss
    subspace_buffer = 1.0 + p.K_buf_ss * p.Buf_ss / (
        state.Ca_ss + p.K_buf_ss
    ) ** 2
    dCa_ss = subspace_total / subspace_buffer

    sr_total = current.I_up - current.I_leak - current.I_rel
    sr_buffer = 1.0 + p.K_buf_SR * p.Buf_SR / (
        state.Ca_SR + p.K_buf_SR
    ) ** 2
    dCa_SR = sr_total / sr_buffer

    kcasr = p.max_SR - (p.max_SR - p.min_SR) / (
        1.0 + (p.EC / state.Ca_SR) ** 2
    )
    k2 = p.k2_prime * kcasr
    dR_prime = -k2 * state.Ca_ss * state.R_prime + p.k4 * (
        1.0 - state.R_prime
    )

    return MainState(dV, dKi, dNai, dCai, dCa_ss, dCa_SR, dR_prime)


def _rush_larsen(old_value: float, steady_state: float, tau: float, dt: float) -> float:
    return steady_state - (steady_state - old_value) * exp(-dt / tau)


def hybrid_forward_euler_step(
    state: MainState,
    gates: Gates,
    dt: float,
    parameters: Parameters = Parameters(),
    phenotype: Phenotype = "epi",
    stimulus_current: float = 0.0,
    sac_coefficient: float = 0.0,
) -> StepResult:
    """Advance one svMP-style FE/Rush--Larsen step.

    The main-state RHS and every gate coefficient use ``state`` and ``gates``
    from the beginning of the step.  The stretch current follows the public
    svMP convention ``Ksac * (V_rest - V_old)``.
    """

    if not isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be finite and positive")

    sac_current = sac_coefficient * (parameters.V_rest - state.V)
    derivative = main_state_rhs(
        state,
        gates,
        parameters,
        stimulus_current=stimulus_current,
        sac_current=sac_current,
    )
    kinetics = gate_kinetics(state, phenotype)

    new_state = MainState(
        **{
            item.name: getattr(state, item.name) + dt * getattr(derivative, item.name)
            for item in fields(MainState)
        }
    )
    new_gates = Gates(
        **{
            item.name: _rush_larsen(
                getattr(gates, item.name),
                kinetics[item.name][0],
                kinetics[item.name][1],
                dt,
            )
            for item in fields(Gates)
        }
    )
    return StepResult(new_state, new_gates)


def named_values(state: MainState, gates: Gates) -> Iterable[Tuple[str, float]]:
    """Yield the canonical seven-main-state then twelve-gate ordering."""

    for item in fields(MainState):
        yield item.name, getattr(state, item.name)
    for item in fields(Gates):
        yield item.name, getattr(gates, item.name)
