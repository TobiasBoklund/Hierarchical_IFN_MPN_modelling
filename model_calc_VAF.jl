     """
   model_calc_VAF(theta,df_p,effect,pDays,p,master_curve,initJAK)

Function for calculating the VAF at the given times according to the model

Input:
    theta: Parameter vector (to be estimated in optimisation/inference)

    df_p: Data frame for the given patient containing patientID, days, IFN, RUX, and JAK

    effect: String containing the effect of the treatment

    pDays: Vector with days at which we want to calculate the VAF

    p: Default parameter tuple

    master_curve: Named tuple holding t, x0, x1, x2, y0, y1, y2, a, s, and VAF from the master curve

    initJAK: Initial condition for the JAK used to calculate the initial condition for the variables

  Output:

    VAF: Vector containing the VAF at the time points in pDays
"""
function model_calc_VAF(theta,df_p,effect,pDays,p,master_curve,initJAK)
    # Extract dosing
    pDosing = df_p[:,[:days,:IFN,:RUX]]

    # Extract fitting parameters
    if effect == "sy0dy1" || effect == "py0py1" || effect == "dy0dy1" || effect == "py0dy0" || effect == "py0dy1" || effect == "dy0py1" || effect == "py1dy1"
        rho1 = theta[1]
        rho2 = theta[2]
    elseif effect == "dy1" || effect == "py0" || effect == "py1" || effect == "dy0" || effect == "dy1IFN"
        rho1 = theta[1]
    elseif effect == "COMBI_py0dy1"
        rho1 = theta[1]
        rho2 = theta[2]
        rho3 = theta[3]
    end

    # Update parameter tuple
    if effect == "sy0dy1" || effect == "py0py1" || effect == "dy0dy1" || effect == "py0dy0" || effect == "py0dy1" || effect == "dy0py1" || effect == "py1dy1"
        change = (dosing=pDosing, rho1=rho1, rho2=rho2, effect=effect)
    elseif effect == "dy1" || effect == "py0" || effect == "py1" || effect == "dy0" || effect == "dy1IFN"
        change = (dosing=pDosing, rho1=rho1, effect=effect)
    elseif effect == "COMBI_py0dy1"
        change = (dosing=pDosing, rho1=rho1, rho2=rho2, rho3=rho3, effect=effect)
    end
    p_treatment = merge(p, change)

    # Find master curve initial conditions using linear interpolation
    # Calculate difference between master curve and initJAK
    diffVAF = master_curve.VAF.-initJAK

    # Identify sign change to find closest two points
    if maximum(diffVAF)>0
        ID = findfirst(x -> x>0, diffVAF)

        # Extract VAFs
        V1 = master_curve.VAF[ID-1]
        V2 = master_curve.VAF[ID]

        # Find best times
        t1 = master_curve.t[ID-1]
        t2 = master_curve.t[ID]

        # Calculate interpolation time
        tstar = (t2-t1)/(V2-V1)*(initJAK-V1)+t1

        # ID = argmin(abs.(master_curve.VAF.-initJAK))

        # Extract initial conditions
        # x00 = master_curve.x0[ID]
        # x10 = master_curve.x1[ID]
        # x20 = master_curve.x2[ID]
        # y00 = master_curve.y0[ID]
        # y10 = master_curve.y1[ID]
        # y20 = master_curve.y2[ID]
        # a0 = master_curve.a[ID]
        # s0 = master_curve.s[ID]

        # Extract initial conditions using interpolation time
        x00 = (master_curve.x0[ID]-master_curve.x0[ID-1])/(t2-t1)*(tstar-t1)+master_curve.x0[ID-1]
        x10 = (master_curve.x1[ID]-master_curve.x1[ID-1])/(t2-t1)*(tstar-t1)+master_curve.x1[ID-1]
        x20 = (master_curve.x2[ID]-master_curve.x2[ID-1])/(t2-t1)*(tstar-t1)+master_curve.x2[ID-1]
        y00 = (master_curve.y0[ID]-master_curve.y0[ID-1])/(t2-t1)*(tstar-t1)+master_curve.y0[ID-1]
        y10 = (master_curve.y1[ID]-master_curve.y1[ID-1])/(t2-t1)*(tstar-t1)+master_curve.y1[ID-1]
        y20 = (master_curve.y2[ID]-master_curve.y2[ID-1])/(t2-t1)*(tstar-t1)+master_curve.y2[ID-1]
        a0 = (master_curve.a[ID]-master_curve.a[ID-1])/(t2-t1)*(tstar-t1)+master_curve.a[ID-1]
        s0 = (master_curve.s[ID]-master_curve.s[ID-1])/(t2-t1)*(tstar-t1)+master_curve.s[ID-1]
    end

    # Collect in one vector
    u0 = [x00,x10,x20,y00,y10,y20,a0,s0]

    # Find maximum time for which you need to solve
    tMin = 0
    tMax = maximum(pDays)

    # Setup and solve ODEproblem - save solution only at relevant times
    tspan = (tMin,tMax)
    prob = ODEProblem(model_rhs_treatment, u0, tspan, p_treatment)
    sol = solve(prob, TRBDF2(), reltol = 1e-3, abstol = 1e-3,  saveat = pDays)

    # Extract solutions
    t = sol.t
    x0 = sol[1,:]
    x1 = sol[2,:]
    x2 = sol[3,:]
    y0 = sol[4,:]
    y1 = sol[5,:]
    y2 = sol[6,:]
    a = sol[7,:]
    s = sol[8,:]
    VAF = y2./(x2+y2)

    # Return VAF
    return VAF, sol
end