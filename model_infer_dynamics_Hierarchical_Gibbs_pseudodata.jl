"""
model_infer_dynamics_Hierarchical_Gibbs_pseudodata(VAF,df_p,effect,p,master_curve)

    Model function for infering VAF dynamics in the model using Turing. This version uses a 
    hierarchical approach designed for Gibbs sampling. This version is made specifically for pseudodata.

    Input:
    VAF: Vector with VAF data for all patients

    df: Data frame for all patients containing patientID, days, IFN, RUX, and JAK

    effect: String containing the name of the effect of the treatment

    p: Default parameter tuple

    master_curve: Named tuple holding t, x0, x1, x2, y0, y1, y2, a, s, and VAF from the master curve

  Output:

    ab: Vector containing the VAF at the time points in pDays
"""
@model function model_infer_dynamics_Hierarchical_Gibbs_pseudodata(VAF,df,effect,p,master_curve)
    
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
    theta_1 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_1)
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
    theta_5 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_5)
    theta_6 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_6)
    theta_7 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_7)
    theta_8 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_8)
    theta_9 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_9)
    theta_10 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_10)
    theta_11 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_11)
    theta_12 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_12)
    theta_13 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_13)
    theta_14 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_14)
    theta_15 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_15)
    theta_16 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_16)
    theta_17 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_17)
    theta_18 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_18)
    theta_19 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_19)
    theta_20 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_20)
    theta_21 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_21)
    theta_22 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_22)
    theta_23 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_23)
    theta_24 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_24)
    theta_25 ~ arraydist([truncated(Normal(rhobar[1],rhotau[1]), lower=0), # py0
    truncated(Normal(rhobar[2],rhotau[2]), lower=0), # dy1
    Uniform(0,0.99999), # InitJAK
    truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)]) # sigma
    push!(Theta, theta_25)
   

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