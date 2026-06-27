include("ImportParameters.jl")

function infant_gutv5!(du, u, p, t)
    # Extract 9 state variables 
    B, P, D, C, H, O, F, G, S = max.(1e-12, u) 
    
    # Unpack parameters 
    (; washout, 
       μ_fit_B, μ_fit_P, μ_fit_D, μ_fit_C,
       μ_basal_B, μ_basal_P, μ_basal_D, μ_basal_C, 
       C_max, K_I, 
       # Interaction Matrix (α_ij = effect of j ON i)
       α_BP, α_BD, α_BC,  # Effects on B
       α_PB, α_PD, α_PC,  # Effects on P
       α_DB, α_DP, α_DC,  # Effects on D
       α_CB, α_CP, α_CD,  # Effects on C
       # Substrate affinities and consumption rates
       K_HB, K_HD, K_FB, K_FD, K_FC, K_O, K_HP, 
       κ_H, κ_F, κ_O, κ_G, κ_S, 
       K_GB, K_GD, K_SB, K_SD, 
       # Conversion efficiencies
       a_HB, a_HD, a_FB, a_FD, a_FC, a_O,    
       a_GB, a_GD, a_SB, a_SD,
       ϕ_H, ϕ_F, ϕ_G, ϕ_S) = p
    
    # Monod Saturation Terms
    sat_H = H / (κ_H + H)
    sat_F = F / (κ_F + F)
    sat_O = O / (κ_O + O)
    sat_G = G / (κ_G + G)
    sat_S = S / (κ_S + S)
    

    # --- GROWTH EQUATIONS ---
    
    # 1. Bifidobacteriaceae (B) on Basal, H, F, G, S
    tB_basal = 1.0 - (μ_basal_B / μ_fit_B)
    tB_H     = 1.0 - ((a_HB * K_HB * sat_H) / μ_fit_B)
    tB_F     = 1.0 - ((a_FB * K_FB * sat_F) / μ_fit_B)
    tB_G     = 1.0 - ((a_GB * K_GB * sat_G) / μ_fit_B)
    tB_S     = 1.0 - ((a_SB * K_SB * sat_S) / μ_fit_B)
    μ_dyn_B  = μ_fit_B * (1.0 - (tB_basal * tB_H * tB_F * tB_G * tB_S)) 

    # 2. Enterobacteriaceae (P) on Basal, O
    tP_basal = 1.0 - (μ_basal_P / μ_fit_P)
    tP_O     = 1.0 - ((a_O * K_O * sat_O) / μ_fit_P)
    μ_dyn_P  = μ_fit_P * (1.0 - (tP_basal * tP_O)) 

    # 3. Bacteroidaceae (D) on Basal, H, F, G, S
    tD_basal = 1.0 - (μ_basal_D / μ_fit_D)
    tD_H     = 1.0 - ((a_HD * K_HD * sat_H) / μ_fit_D)
    tD_F     = 1.0 - ((a_FD * K_FD * sat_F) / μ_fit_D)
    tD_G     = 1.0 - ((a_GD * K_GD * sat_G) / μ_fit_D)
    tD_S     = 1.0 - ((a_SD * K_SD * sat_S) / μ_fit_D)
    μ_dyn_D  = μ_fit_D * (1.0 - (tD_basal * tD_H * tD_F * tD_G * tD_S)) 

    # 4. Clostridiales (C) on Basal and F ONLY
    tC_basal = 1.0 - (μ_basal_C / μ_fit_C)
    tC_F     = 1.0 - ((a_FC * K_FC * sat_F) / μ_fit_C)
    μ_dyn_C  = μ_fit_C * (1.0 - (tC_basal * tC_F)) 
    
    # --- CARRYING CAPACITY LIMITS ---
    
    limit_B = 1.0 - (B + (α_BP * P) + (α_BD * D) + (α_BC * C)) / C_max
    limit_P = 1.0 - (P + (α_PB * B) + (α_PD * D) + (α_PC * C)) / C_max 
    limit_D = 1.0 - (D + (α_DB * B) + (α_DP * P) + (α_DC * C)) / C_max 
    limit_C = 1.0 - (C + (α_CB * B) + (α_CP * P) + (α_CD * D)) / C_max 

    # --- DIFFERENTIAL EQUATIONS ---
    
    # 1. Bacterial Populations
    du[1] = μ_dyn_B * B * limit_B - (washout * B)
    du[2] = μ_dyn_P * P * limit_P - (washout * P)
    du[3] = μ_dyn_D * D * limit_D - (washout * D)
    du[4] = μ_dyn_C * C * limit_C - (washout * C)
    
    # 2. Substrate/Resource Flux 
    
    du[5] = -(K_HB * B * sat_H) - (K_HD * D * sat_H) - (ϕ_H * H)                   # HMO (H)
    du[6] = - (K_O * P * sat_O)                                                    # Oxygen (O)
    du[7] = -(K_FB * B * sat_F) - (K_FD * D * sat_F) - (K_FC * C * sat_F) - (ϕ_F * F) # Fiber (F) - Now includes C
    du[8] = -(K_GB * B * sat_G) - (K_GD * D * sat_G) - (ϕ_G * G)                   # scGOS (G)
    du[9] = -(K_SB * B * sat_S) - (K_SD * D * sat_S) - (ϕ_S * S)                   # lcFOS (S)
end

# -------------------------------------------------------
# 1. Weaning Function W(t)
# -------------------------------------------------------
# Standard logistic curve representing the introduction of solid food.
# Transitions smoothly from 0.0 (pure milk) to 1.0 (pure solid food) around Day 180.
function W(t)
    return 1.0 / (1.0 + exp(-0.05 * (t - 180.0)))
end
