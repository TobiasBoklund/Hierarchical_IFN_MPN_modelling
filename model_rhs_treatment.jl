"""
    model_rhs_treatment(du, u, p, t)

    Function for calculating rhs of the model for use in solving ODE problem during treatment

        Input: 

            du: In-place holder for rhs of ODE system

            u: Vector with values of variables of ODE system

            p: Named-tuple with parameter values for ODE system

            t: Time

        Output:

            du: rhs of ODE system
"""
function model_rhs_treatment(du, u, p, t)
    # Function for calculating rhs of model for use in solving ODE problem during treatment
    # Input: 
    #   du: In-place holder for rhs of ODE system
    #   u: Vector with values of variables of ODE system
    #   p: Named-tuple with parameter values for ODE system
    #   t: Time
    # Output:
    #   du: rhs of ODE system

    # Extract variables
    x0 = u[1]
    x1 = u[2]
    x2 = u[3]
    y0 = u[4]
    y1 = u[5]
    y2 = u[6]
    a = u[7]
    s = u[8]

    # Extract parameters
    alphax0 = p.alphax0
    px0 = p.px0
    dx0 = p.dx0
    rm = p.rm
    alphay0 = p.alphay0
    py0 = p.py0
    dy0 = p.dy0
    Ax0 = p.Ax0
    Ay0 = p.Ay0
    alphax1 = p.alphax1
    px1 = p.px1
    dx1 = p.dx1
    alphay1 = p.alphay1
    py1 = p.py1
    dy1 = p.dy1
    Ax1 = p.Ax1
    Ay1 = p.Ay1
    dx2 = p.dx2
    dy2 = p.dy2
    ea = p.ea
    rs = p.rs
    es = p.es
    I = p.I
    cxx = p.cxx
    cxy = p.cxy
    cyx = p.cyx
    cyy = p.cyy
    sx0 = p.sx0
    sy0 = p.sy0
    dosing = p.dosing
    ct = dosing.days
    cI = dosing.IFN
    cR = dosing.RUX
    effect = p.effect
    if effect == "sy0dy1" || effect == "py0py1" || effect == "dy0dy1" || effect == "py0dy0" || effect == "py0dy1" || 
       effect == "dy0py1" || effect == "py1dy1"
        rho1 = p.rho1
        rho2 = p.rho2
    elseif effect == "dy1" || effect == "py0" || effect == "py1" || effect == "dy0" || effect == "dy1IFN"
        rho1 = p.rho1
    elseif effect == "COMBI_py0dy1"
        rho1 = p.rho1
        rho2 = p.rho2
        rho3 = p.rho3
    end

    # Extract dose at current time point
    tdiff = ct.-t
    tID = tdiff.<=0
    if sum(tID) == 0
        cnowI = 0
        cnowR = 0
    else
        ID = findlast(tID)
        cnowI = cI[ID]
        cnowR = cR[ID]
    end

    # Update parameters
    if effect == "sy0dy1"
        sy0 = (1+rho1*cnowR)*sy0
        dy1 = (1+rho2*cnowR)*dy1
    elseif effect == "dy1"
        dy1 = (1+rho1*cnowR)*dy1
    elseif effect == "py0py1"
        py0 = (1-2*(1/(1+exp(-rho1*cnowI))-1/2))*py0
        py1 = (1-2*(1/(1+exp(-rho2*cnowI))-1/2))*py1
    elseif effect == "py0"
        py0 = (1-2*(1/(1+exp(-rho1*cnowI))-1/2))*py0
    elseif effect == "py1"
        py1 = (1-2*(1/(1+exp(-rho1*cnowI))-1/2))*py1
    elseif effect == "dy0dy1"
        dy0 = rho1*cnowI*1/500
        dy1 = (1+rho2*cnowI)*dy1
    elseif effect == "dy0"
        dy0 = rho1*cnowI*1/500
    elseif effect == "dy1IFN"
        dy1 = (1+rho1*cnowI)*dy1
    elseif effect == "py0dy0"
        py0 = (1-2*(1/(1+exp(-rho1*cnowI))-1/2))*py0
        dy0 = rho2*cnowI*1/500
    elseif effect == "py0dy1"
        py0 = (1-2*(1/(1+exp(-rho1*cnowI))-1/2))*py0
        dy1 = (1+rho2*cnowI)*dy1
    elseif effect == "dy0py1"
        dy0 = rho1*cnowI*1/500
        py1 = (1-2*(1/(1+exp(-rho2*cnowI))-1/2))*py1
    elseif effect == "py1dy1"
        py1 = (1-2*(1/(1+exp(-rho1*cnowI))-1/2))*py1
        dy1 = (1+rho2*cnowI)*dy1
    elseif effect == "COMBI_py0dy1"
        py0 = (1-2*(1/(1+exp(-rho1*cnowI))-1/2))*py0
        dy1 = (1+rho2*cnowI+rho3*cnowR)*dy1
    end

    # Calculate crowding functions phi_x and phi_y
    phix = 1/(1+cxx*x0+cxy*y0)
    phiy = 1/(1+cyx*x0+cyy*y0)
    
    # Rhs of ODEs
    x0dot = alphax0*(2*px0*phix*s/(sx0+s)-1)*x0-dx0*x0-rm*s*x0
    x1dot = alphax1*(2*px1-1)*x1+2*Ax0*alphax0*(1-px0*phix*s/(sx0+s))*x0-dx1*x1
    x2dot = 2*Ax1*alphax1*(1-px1)*x1-dx2*x2
    y0dot = alphay0*(2*py0*phiy*s/(sy0+s)-1)*y0-dy0*y0+rm*s*y0
    y1dot = alphay1*(2*py1-1)*y1+2*Ay0*alphay0*(1-py0*phiy*s/(sy0+s))*y0-dy1*y1 
    y2dot = 2*Ay1*alphay1*(1-py1)*y1-dy2*y2
    adot = dx0*x0+dy0*y0+dx1*x1+dy1*y1+dx2*x2+dy2*y2-ea*a*s
    sdot = rs*a-es*s+I
    
    # Collects rhses 
    du[1] = x0dot
    du[2] = x1dot
    du[3] = x2dot
    du[4] = y0dot
    du[5] = y1dot
    du[6] = y2dot
    du[7] = adot
    du[8] = sdot
end