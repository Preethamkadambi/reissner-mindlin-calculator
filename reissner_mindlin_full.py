import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# HEADER & CONFIG
# ==========================================
st.set_page_config(layout="wide", page_title="Full Reissner-Mindlin Plate Theory")

st.title("Reissner-Mindlin Plate Theory: Complete Implementation")
st.markdown("""
This application digitizes the entire derivation from **Aircraft Structures - Plates - Reissner-Mindlin Theory** (University of Liège).
It follows the document structure strictly from **Kinematics** to **Equilibrium**.
""")

# ==========================================
# SIDEBAR: MATERIAL & GEOMETRY (Pages 2-3)
# ==========================================
st.sidebar.header("Global Inputs (Pages 2-3)")
st.sidebar.info("Reference: *StructAeroPlates.pdf*")

# Page 2: Linear elastic, homogeneous, isotropic material
E = st.sidebar.number_input("Young's Modulus (E) [Pa]", value=70e9, format="%.2e")
nu = st.sidebar.number_input("Poisson's Ratio (ν)", value=0.3, min_value=0.0, max_value=0.499)
rho = st.sidebar.number_input("Density (ρ) [kg/m³]", value=2700.0)

# Page 3: Geometry
h0 = st.sidebar.number_input("Thickness (h₀) [m]", value=0.01, format="%.4f")
k_shear = st.sidebar.number_input("Shear Correction Factor (A'/A) [Page 38/44]", value=5/6)

st.sidebar.markdown("---")

# ==========================================
# UTILS
# ==========================================
def latex_log(title, formula, result, page):
    """Logs a calculation step with LaTeX and page reference."""
    with st.expander(f"{title} (Page {page})", expanded=True):
        st.latex(formula)
        st.write(f"**Calculated Value:** {result}")

# ==========================================
# MAIN TABS
# ==========================================
tab1, tab2, tab3, tab4 = st.tabs([
    "1. Kinematics (Deformations)", 
    "2. Material Law (Stresses)", 
    "3. Resultants (Forces/Moments)",
    "4. Equilibrium (Balance)"
])

# ==========================================
# TAB 1: KINEMATICS (Pages 3 - 35)
# ==========================================
with tab1:
    st.header("1. Kinematics: Displacement to Strain")
    st.markdown("Calculates strains $\epsilon, \kappa, \gamma$ from displacement gradients.")
    
    col_in, col_out = st.columns(2)
    
    with col_in:
        st.subheader("Inputs: Displacement Gradients")
        st.markdown("**Membrane Gradients (Page 31):** $u_{\\alpha,\\beta}$")
        u1_1 = st.number_input("∂u₁/∂x (u1,1)", value=0.001, format="%.6f")
        u2_2 = st.number_input("∂u₂/∂y (u2,2)", value=0.0005, format="%.6f")
        u1_2 = st.number_input("∂u₁/∂y (u1,2)", value=0.000, format="%.6f")
        u2_1 = st.number_input("∂u₂/∂x (u2,1)", value=0.000, format="%.6f")
        
        st.markdown("**Rotation Gradients (Page 32):** $\Delta t_{\\alpha,\\beta}$")
        dt1_1 = st.number_input("∂Δt₁/∂x (Δt1,1)", value=0.05, format="%.6f")
        dt2_2 = st.number_input("∂Δt₂/∂y (Δt2,2)", value=0.01, format="%.6f")
        dt1_2 = st.number_input("∂Δt₁/∂y (Δt1,2)", value=0.00, format="%.6f")
        dt2_1 = st.number_input("∂Δt₂/∂x (Δt2,1)", value=0.00, format="%.6f")
        
        st.markdown("**Shear Gradients (Page 33):** $u_{3,\\alpha} + \Delta t_{\\alpha}$")
        u3_1 = st.number_input("∂u₃/∂x (u3,1)", value=0.00, format="%.6f")
        u3_2 = st.number_input("∂u₃/∂y (u3,2)", value=0.00, format="%.6f")
        dt1 = st.number_input("Rotation Δt₁", value=0.0001, format="%.6f")
        dt2 = st.number_input("Rotation Δt₂", value=0.0000, format="%.6f")

    with col_out:
        st.subheader("Calculated Deformation Modes")
        
        # --- Membrane Mode (Page 31) ---
        eps_11 = u1_1
        eps_22 = u2_2
        eps_12 = 0.5 * (u1_2 + u2_1)
        
        latex_log("Membrane Strain", 
                  r"\epsilon_{\alpha\beta} = \frac{u_{\alpha,\beta} + u_{\beta,\alpha}}{2} \quad [Page 31]",
                  f"ε₁₁={eps_11:.6f}, ε₂₂={eps_22:.6f}, ε₁₂={eps_12:.6f}", 31)
        
        # --- Bending Mode (Page 32) ---
        kappa_11 = dt1_1
        kappa_22 = dt2_2
        kappa_12 = 0.5 * (dt1_2 + dt2_1)
        
        latex_log("Curvature", 
                  r"\kappa_{\alpha\beta} = \frac{\Delta t_{\alpha,\beta} + \Delta t_{\beta,\alpha}}{2} \quad [Page 32]",
                  f"κ₁₁={kappa_11:.6f}, κ₂₂={kappa_22:.6f}, κ₁₂={kappa_12:.6f}", 32)

        # --- Shear Mode (Page 33) ---
        # Note: Page 33 defines gamma = 2 * delta. Usually gamma is engineering shear.
        # Formula: 2*delta_alpha = gamma_alpha = u3,alpha + Delta_t_alpha
        gamma_1 = u3_1 + dt1
        gamma_2 = u3_2 + dt2
        
        latex_log("Transverse Shear Strain", 
                  r"\gamma_{\alpha} = u_{3,\alpha} + \Delta t_{\alpha} \quad [Page 33]",
                  f"γ₁={gamma_1:.6f}, γ₂={gamma_2:.6f}", 33)
        
        # Store in session state for next tabs
        st.session_state['strain'] = {'e11': eps_11, 'e22': eps_22, 'e12': eps_12}
        st.session_state['curv'] = {'k11': kappa_11, 'k22': kappa_22, 'k12': kappa_12}
        st.session_state['shear'] = {'g1': gamma_1, 'g2': gamma_2}

# ==========================================
# TAB 2: MATERIAL LAW & STRESSES (Pages 36 - 41)
# ==========================================
with tab2:
    st.header("2. Material Law: Stresses & Through-Thickness Elongation")
    
    if 'strain' not in st.session_state:
        st.warning("Please calculate deformations in Tab 1 first.")
    else:
        s = st.session_state['strain']
        k = st.session_state['curv']
        g = st.session_state['shear']
        
        # --- Coefficients ---
        # Page 2/36: Hooke's Law Coefficients
        C1 = (E * nu) / ((1 + nu) * (1 - 2 * nu))
        C2 = E / (1 + nu)
        
        # --- Through-Thickness Elongation (Page 39/40) ---
        st.subheader("Through-Thickness Elongation ($\lambda_h$)")
        # Formula Page 40: lambda_h = (nu / (nu - 1)) * (u_gamma,gamma + xi^3 * Delta_t_gamma,gamma)
        # Note: u_gamma,gamma is trace of strain (div u)
        div_u = s['e11'] + s['e22']
        div_dt = k['k11'] + k['k22']
        
        xi3 = st.slider("Select thickness coordinate (ξ³)", -h0/2, h0/2, 0.0, format="%.4f")
        
        lambda_h = (nu / (nu - 1)) * (div_u + xi3 * div_dt)
        
        latex_log("Through-Thickness Elongation", 
                  r"\lambda_h(\xi^3) = \frac{\nu}{\nu - 1} (u_{\gamma,\gamma} + \xi^3 \Delta t_{\gamma,\gamma}) \quad [Page 40]",
                  f"λ_h at z={xi3:.4f} is {lambda_h:.6e}", 40)
        
        # --- Stress Calculation (Page 41) ---
        st.subheader("Stress Recovery ($\sigma^{\alpha\\beta}$)")
        st.markdown(r"Total Stress $\sigma^{\alpha\beta} = \sigma^{\alpha\beta}_{\tilde{n}} + \xi^3 \sigma^{\alpha\beta}_{\tilde{m}}$")
        
        # Membrane Stress (Page 41)
        # sigma_n_ab = E*nu/(1-nu^2) * div_u * delta_ab + E/(1+nu) * eps_ab
        coeff_n1 = (E * nu) / (1 - nu**2)
        coeff_n2 = E / (1 + nu)
        
        sig_n_11 = coeff_n1 * div_u + coeff_n2 * s['e11']
        sig_n_22 = coeff_n1 * div_u + coeff_n2 * s['e22']
        sig_n_12 = coeff_n2 * s['e12']
        
        # Bending Stress (Page 41)
        # sigma_m_ab = same formula but with kappa and div_dt
        sig_m_11 = coeff_n1 * div_dt + coeff_n2 * k['k11']
        sig_m_22 = coeff_n1 * div_dt + coeff_n2 * k['k22']
        sig_m_12 = coeff_n2 * k['k12']
        
        # Total Stress at xi3
        total_sig_11 = sig_n_11 + xi3 * sig_m_11
        
        latex_log("Total Stress Formula", 
                  r"\sigma^{\alpha\beta} = \left[ \frac{E\nu}{1-\nu^2} u_{\gamma,\gamma} \delta^{\alpha\beta} + \frac{E}{1+\nu} \epsilon_{\alpha\beta} \right] + \xi^3 \left[ \dots \kappa \dots \right] \quad [Page 41]",
                  f"σ₁₁ at z={xi3:.4f} is {total_sig_11/1e6:.2f} MPa", 41)

        # Plotting Stress Distribution
        st.markdown("### Stress Distribution Through Thickness")
        z_vals = np.linspace(-h0/2, h0/2, 100)
        sig_11_vals = sig_n_11 + z_vals * sig_m_11
        
        fig, ax = plt.subplots()
        ax.plot(sig_11_vals/1e6, z_vals*1000)
        ax.set_xlabel("Stress σ₁₁ [MPa]")
        ax.set_ylabel("Thickness Coordinate ξ³ [mm]")
        ax.set_title("Stress σ₁₁ vs Thickness")
        ax.grid(True)
        ax.axvline(0, color='k', linestyle='--')
        ax.axhline(0, color='k', linestyle='--')
        st.pyplot(fig)

# ==========================================
# TAB 3: RESULTANTS (Pages 42 - 49)
# ==========================================
with tab3:
    st.header("3. Resultant Forces & Moments")
    
    if 'strain' not in st.session_state:
        st.warning("Calculate Kinematics first.")
    else:
        # --- Membrane Resultants (Page 43) ---
        st.subheader("Membrane Forces (Page 43)")
        # H_n tensor
        factor_n = (h0 * E) / (1 - nu**2)
        
        N11 = factor_n * (s['e11'] + nu * s['e22'])
        N22 = factor_n * (s['e22'] + nu * s['e11'])
        N12 = factor_n * (1 - nu) * s['e12']
        
        latex_log("Membrane Resultant", 
                  r"\tilde{n}^{\alpha\beta} = \int \sigma^{\alpha\beta} d\xi^3 = \mathcal{H}_n \epsilon \quad [Page 43]",
                  f"N₁₁ = {N11:.2f} N/m, N₂₂ = {N22:.2f} N/m, N₁₂ = {N12:.2f} N/m", 43)

        # --- Shear Resultants (Page 44) ---
        st.subheader("Shear Forces (Page 44)")
        # Note: Page 44 uses A'/A shear correction
        factor_q = (E * h0) / (1 + nu) * k_shear
        
        Q1 = 0.5 * factor_q * g['g1']
        Q2 = 0.5 * factor_q * g['g2']
        
        latex_log("Shear Resultant", 
                  r"\tilde{q}^{\alpha} = \frac{E h_0}{1+\nu} \frac{A'}{A} \frac{\gamma_{\alpha}}{2} \quad [Page 44]",
                  f"Q₁ = {Q1:.2f} N/m, Q₂ = {Q2:.2f} N/m", 44)
        
        # --- Bending Resultants (Page 47) ---
        st.subheader("Bending Moments (Page 47)")
        D = (h0**3 * E) / (12 * (1 - nu**2))
        
        M11 = D * (k['k11'] + nu * k['k22'])
        M22 = D * (k['k22'] + nu * k['k11'])
        M12 = D * (1 - nu) * k['k12']
        
        latex_log("Bending Resultant", 
                  r"\tilde{m}^{\alpha\beta} = \int \sigma^{\alpha\beta} \xi^3 d\xi^3 = \mathcal{H}_m \kappa \quad [Page 47]",
                  f"M₁₁ = {M11:.2f} N, M₂₂ = {M22:.2f} N, M₁₂ = {M12:.2f} N", 47)
        
        st.session_state['resultants'] = {'Q1': Q1, 'Q2': Q2}

# ==========================================
# TAB 4: EQUILIBRIUM (Pages 50 - 52)
# ==========================================
with tab4:
    st.header("4. Equilibrium Equations")
    st.markdown("Verifies the balance between internal forces gradients and external loads.")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Inputs: Gradients & Loads")
        st.markdown("**Gradients of Moments (Page 51):**")
        # Need div(m) to check shear balance
        dm11_dx = st.number_input("∂M₁₁/∂x", value=100.0)
        dm12_dy = st.number_input("∂M₁₂/∂y", value=50.0)
        
        dm22_dy = st.number_input("∂M₂₂/∂y", value=100.0)
        dm21_dx = st.number_input("∂M₂₁/∂x", value=50.0)
        
        st.markdown("**External Loads (Page 51):**")
        m_bar_1 = st.number_input("Ext Moment m_bar_1", value=0.0)
        m_bar_2 = st.number_input("Ext Moment m_bar_2", value=0.0)
        
    with col2:
        st.subheader("Equilibrium Check")
        
        if 'resultants' not in st.session_state:
            st.warning("Calculate resultants first.")
        else:
            Q = st.session_state['resultants']
            
            # Equation for rotation (Page 51, assuming static: I_p_dot = 0)
            # 0 = m_bar - Q + div(m)
            # Component 1: m_bar_1 - Q1 + (M11,1 + M12,2) = 0?
            
            # Note on Page 51: The equation is: 
            # -Ip * d(dt)/dt = m_bar - Q + div(m)
            # For static: Q_alpha = m_bar_alpha + m_alpha_beta,beta
            
            # LHS (Internal Shear)
            LHS_1 = Q['Q1']
            LHS_2 = Q['Q2']
            
            # RHS (Gradients + Ext)
            RHS_1 = m_bar_1 + (dm11_dx + dm12_dy)
            RHS_2 = m_bar_2 + (dm21_dx + dm22_dy)
            
            latex_log("Equilibrium Eq (Rotation)", 
                      r"\tilde{q}^{\alpha} = \tilde{\tilde{m}}_{\alpha} + \tilde{m}^{\alpha\beta}_{,\beta} \quad [Page 51]",
                      f"X-Dir: Q₁ ({LHS_1:.2f}) vs RHS ({RHS_1:.2f}) \n Y-Dir: Q₂ ({LHS_2:.2f}) vs RHS ({RHS_2:.2f})", 51)
            
            delta_1 = LHS_1 - RHS_1
            if abs(delta_1) < 1e-5:
                st.success("X-Rotation Equilibrium Satisfied")
            else:
                st.error(f"X-Rotation Unbalanced (Residual: {delta_1:.4f})")

st.markdown("---")
[cite_start]st.caption("Implementation based on 'Aircraft Structures - Plates - Reissner-Mindlin Theory', University of Liège [cite: 1-1286].")