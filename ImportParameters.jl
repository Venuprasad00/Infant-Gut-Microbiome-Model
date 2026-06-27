using CSV
using DataFrames

# =======================================================
# 1. LOAD ESTIMATED PARAMETERS 
# =======================================================
param_file = joinpath("Parameters", "Estimated_Parameters.csv")
df_params = CSV.read(param_file, DataFrame)
est = Vector(df_params[:, :Estimate])

# Growth rates
μ_fit_B = est[1]; μ_fit_P = est[2]; μ_fit_D = est[3]; μ_fit_C = est[4]

# Interaction terms
α_BP = est[5];  α_BD = est[6];  α_BC = est[7]
α_PB = est[8];  α_PD = est[9];  α_PC = est[10]
α_DB = est[11]; α_DP = est[12]; α_DC = est[13]
α_CB = est[14]; α_CP = est[15]; α_CD = est[16]

# Calculate Basal Rates strictly as 50% of maximums
μ_basal_B = 0.50 * μ_fit_B
μ_basal_P = 0.50 * μ_fit_P
μ_basal_D = 0.50 * μ_fit_D
μ_basal_C = 0.50 * μ_fit_C

# =======================================================
# 2. LOAD SYSTEM PARAMETERS & CONSTANTS 
# =======================================================
sys_file = joinpath("Parameters", "Model_Parameters.csv")
df_sys = CSV.read(sys_file, DataFrame)

# Convert to Dictionary for robust name-based extraction
sys = Dict(df_sys.Parameter .=> df_sys.Value)

# =======================================================
# 3. PARAMETER DEFINITION (NamedTuple)
# =======================================================
params = (
    # A. INITIAL CONDITIONS
    b0_v = sys["b0_v"], p0_v = sys["p0_v"], d0_v = sys["d0_v"], c0_v = sys["c0_v"],
    b0_c = sys["b0_c"], p0_c = sys["p0_c"], d0_c = sys["d0_c"], c0_c = sys["c0_c"],
    
    # B. GLOBAL LIMITS
    washout = sys["washout"],
    C_max   = sys["C_max"],
    
    # C. GROWTH RATES
    μ_fit_B = μ_fit_B, μ_fit_P = μ_fit_P, μ_fit_D = μ_fit_D, μ_fit_C = μ_fit_C,
    μ_basal_B = μ_basal_B, μ_basal_P = μ_basal_P, μ_basal_D = μ_basal_D, μ_basal_C = μ_basal_C,

    # D. COMPETITIVE INTERACTION MATRIX
    α_BP = α_BP, α_BD = α_BD, α_BC = α_BC,
    α_PB = α_PB, α_PD = α_PD, α_PC = α_PC,
    α_DB = α_DB, α_DP = α_DP, α_DC = α_DC,
    α_CB = α_CB, α_CP = α_CP, α_CD = α_CD,

    # E. SUBSTRATE CONSUMPTION RATES (K) 
    K_HB = sys["K_HB"], K_HD = sys["K_HD"], K_HP = sys["K_HP"], 
    K_FB = sys["K_FB"], K_FD = sys["K_FD"], K_FC = sys["K_FC"],             
    K_O  = sys["K_O"],                              
    K_GB = sys["K_GB"], K_GD = sys["K_GD"],                 
    K_SB = sys["K_SB"], K_SD = sys["K_SD"],                

    # F. MONOD HALF-SATURATION CONSTANTS (κ)
    κ_H = sys["kappa_H"], κ_F = sys["kappa_F"], κ_O = sys["kappa_O"], 
    κ_G = sys["kappa_G"], κ_S = sys["kappa_S"], K_I = sys["K_I"],

    # G. BIOMASS YIELD COEFFICIENTS (a)
    a_HB = (1.00 * (μ_fit_B - μ_basal_B)) / sys["K_HB"],
    a_FB = (0.20 * (μ_fit_B - μ_basal_B)) / sys["K_FB"],
    a_GB = (0.80 * (μ_fit_B - μ_basal_B)) / sys["K_GB"],
    a_SB = (0.80 * (μ_fit_B - μ_basal_B)) / sys["K_SB"],

    a_HD = (0.20 * (μ_fit_D - μ_basal_D)) / sys["K_HD"],
    a_FD = (1.00 * (μ_fit_D - μ_basal_D)) / sys["K_FD"],
    a_GD = (1.00 * (μ_fit_D - μ_basal_D)) / sys["K_GD"],
    a_SD = (1.00 * (μ_fit_D - μ_basal_D)) / sys["K_SD"],

    a_FC = (1.00 * (μ_fit_C - μ_basal_C)) / sys["K_FC"],
    a_O  = (1.00 * (μ_fit_P - μ_basal_P)) / sys["K_O"],

    # H. SUBSTRATE WASHOUT / DECAY RATES (ϕ)
    ϕ_H = sys["phi_H"], 
    ϕ_F = sys["phi_F"], 
    ϕ_G = sys["phi_G"], 
    ϕ_S = sys["phi_S"]
)