"""
model_infer_dynamics_post_Hierarchical(VAF,df_p,effect,p,master_curve,estInitJAK=true)

    Model function for infering VAF dynamics in the model post hierarchical using Turing

    Input:
    VAF: Vector with VAF data for given patient

    df_p: Data frame for the given patient containing patientID, days, IFN, RUX, and JAK

    effect: String containing the name of the effect of the treatment

    pDays: Vector with days at which we want to calculate the VAF - should fit with VAF

    p: Default parameter tuple

    master_curve: Named tuple holding t, x0, x1, x2, y0, y1, y2, a, s, and VAF from the master curve

    estInitJAK: Boolean telling wether or not the initial condition should be estimated. Default is true.

  Output:

    ab: Vector containing the VAF at the time points in pDays
"""
@model function model_infer_dynamics_post_Hierarchical(VAF,df_p,effect,p,master_curve,estInitJAK=true)
    # Parameters - prior distributions
    if effect == "sy0dy1"
        rho1 ~ Uniform(0,4)
        rho2 ~ Uniform(0,2)
    elseif effect == "dy1"
        rho1 ~ Uniform(0,1)
    elseif effect == "py0py1"
        rho1 ~ Uniform(0,4)
        rho2 ~ Uniform(0,4)
    elseif effect == "py0"
        rho1 ~ Uniform(0,4)
    elseif effect == "py1"
        rho1 ~ Uniform(0,4)
    elseif effect == "dy0dy1"
        rho1 ~ Uniform(0,1)
        rho2 ~ Uniform(0,1) 
    elseif effect == "dy0"
        rho1 ~ Uniform(0,1)
    elseif effect == "dy1IFN"
        rho1 ~ Uniform(0,1)
    elseif effect == "py0dy0"
        rho1 ~ Uniform(0,4)
        rho2 ~ Uniform(0,1)
    elseif effect == "py0dy1"
        rho1 ~ truncated(Normal(0.0402, 0.283), lower=0.0)
        rho2 ~ truncated(Normal(0.0684, 0.312), lower=0.0)
    elseif effect == "dy0py1"
        rho1 ~ Uniform(0,1)
        rho2 ~ Uniform(0,4)
    elseif effect == "py1dy1"
        rho1 ~ Uniform(0,4)
        rho2 ~ Uniform(0,1)
    elseif effect == "COMBI_py0dy1"
        rho1 ~ Uniform(0,1)
        rho2 ~ Uniform(0,1)
        rho3 ~ Uniform(0,1)
    end

    # Collect in one vector
    if effect == "sy0dy1" || effect == "py0py1" || effect == "dy0dy1" || effect == "py0dy0" || effect == "py0dy1" || effect == "dy0py1" || effect == "py1dy1"
        rho = [rho1, rho2]
    elseif effect == "dy1" || effect == "py0" || effect == "py1" || effect == "dy0" || effect == "dy1IFN"
        rho = rho1
    elseif effect == "COMBI_py0dy1"
        rho = [rho1, rho2, rho3]
    end

    # Distribution for init_JAK if wanted
    if estInitJAK
        # initJAK ~ truncated(Normal(df_p.JAK[1],0.02),0,1)
        initJAK ~ Uniform(0,0.99999)
    else
        initJAK = df_p.JAK[1]
    end
    
    # Extract treatment days - remove NaN observations
    pDays = df_p.days
    pDays = pDays[.!isnan.(df_p.JAK)]

    # Solve model with treatment
    model_VAF, sol = model_calc_VAF(rho,df_p,effect,pDays,p,master_curve,initJAK)

    # Error standard deviation
    sigma ~ truncated(Cauchy(0,0.3), lower=0.0, upper=1.0)

    # Statistical model for errors
    for i in eachindex(VAF)
        VAF[i] ~ truncated(Normal(model_VAF[i], sigma), lower=0.0, upper=1.0)
    end

end