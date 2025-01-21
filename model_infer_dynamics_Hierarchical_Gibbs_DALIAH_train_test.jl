"""
model_infer_dynamics_Hierarchical_Gibbs_DALIAH_train_test(VAF,df_p,effect,p,master_curve)

    Model function for infering VAF dynamics in the model using Turing. This version uses a 
    hierarchical approach designed for Gibbs sampling. This version is made specifically for the DALIAH cohort
    with train test splitting.

    Input:
    VAF: Vector with VAF data for all patients

    df: Data frame for all patients containing patientID, days, IFN, RUX, and JAK

    effect: String containing the name of the effect of the treatment

    p: Default parameter tuple

    master_curve: Named tuple holding t, x0, x1, x2, y0, y1, y2, a, s, and VAF from the master curve

  Output:

    ab: Vector containing the VAF at the time points in pDays
"""
@model function model_infer_dynamics_Hierarchical_Gibbs_DALIAH_train_test(VAF,df,effect,p,master_curve)
    
    # Unique patients
    unique_patients = unique(df.patientID) 

    # Extract number of patients
    P = length(unique_patients)

    # Extract number of observations
    N = size(df,1)

    # Hyper parameters
    rhobar ~ arraydist([Uniform(0,0.5),Uniform(0,0.5)])
    rhotau ~ arraydist([truncated(Cauchy(0,1), lower=0.0), truncated(Cauchy(0,1), lower=0.0)])

    Theta = []
    theta_2 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_2)
    theta_3 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_3)
    theta_4 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_4)
    # theta_6 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_6)
    theta_11 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_11)
    theta_14 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_14)
    # theta_16 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_16)
    theta_20 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_20)
    theta_22 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_22)
    # theta_25 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_25)
    # theta_27 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_27)
    theta_34 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_34)
    theta_35 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_35)
    theta_38 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_38)
    theta_42 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_42)
    # theta_44 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_44)
    # theta_45 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_45)
    theta_48 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_48)
    theta_61 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_61)
    theta_65 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_65)
    theta_66 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_66)
    theta_69 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_69)
    theta_73 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_73)
    theta_76 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_76)
    theta_77 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_77)
    theta_80 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_80)
    # theta_82 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_82)
    theta_84 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_84)
    theta_85 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_85)
    theta_87 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_87)
    # theta_102 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_102)
    theta_103 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_103)
    # theta_108 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_108)
    theta_110 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_110)
    theta_111 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_111)
    theta_112 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_112)
    theta_114 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_114)
    theta_119 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_119)
    theta_121 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_121)
    # theta_122 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_122)
    # theta_126 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_126)
    theta_128 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_128)
    theta_132 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_132)
    theta_136 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_136)
    theta_137 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_137)
    theta_139 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_139)
    theta_141 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_141)
    theta_142 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_142)
    theta_145 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_145)
    theta_147 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_147)
    theta_149 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_149)
    theta_153 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_153)
    theta_155 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_155)
    theta_159 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_159)
    # theta_160 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_160)
    theta_161 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_161)
    theta_163 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_163)
    theta_164 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_164)
    # theta_179 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_179)
    theta_181 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_181)
    theta_184 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_184)
    theta_187 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_187)
    theta_193 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_193)
    theta_194 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_194)
    theta_195 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_195)
    theta_196 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_196)
    theta_198 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_198)
    theta_199 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_199)
    # theta_200 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    # truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    # Uniform(0,0.99999), # InitJAK
    # truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    # push!(Theta, theta_200)
    theta_202 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_202)
   

   # Set counter to 0
   counter = 0
    
    # Calculate VAFs for patients in loop
    for i=1:P
        # Extract information for one patient
        pID = pID = unique_patients[i]
        ID = df.patientID .== pID
        df_p = df[df.patientID .== pID, :]

        # Extract treatment days
        pDays = df_p.days[.!isnan.(df_p.JAK)]

        # Extract number of data points
        n = length(pDays)

        # Extract rho values
        temp = Theta[i]
        rho_p = temp[1:2]
        initJAK = temp[3]
        sigma = temp[4]

        # Solve model with treatment
        model_VAF, sol = model_calc_VAF(rho_p,df_p,effect,pDays,p,master_curve,initJAK)

        # Statistical model for errors
        for j in 1:length(model_VAF)
            # Update counter
            counter = counter+1

            # Save VAF
            VAF[counter] ~ truncated(Normal(model_VAF[j], sigma), lower=0.0, upper=1.0)
        end
    end

end