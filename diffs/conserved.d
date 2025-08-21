// Alternative formulation derived in reformulate_jac.max
// This has been tested and works, but doesn't have heat or friction in it

double[3] increment_conserved(double x, double A, double dx, double dA, double Hdot, double f, ref Primitives P0, ref GasModel gm, ref GasState gs){

    double dAdx = dA/dx;
    double Ast = dAdx/(A/dx + dAdx); // FIXME: Note that the A here should be A1...

    double Fmass1 = P0.rho*P0.v;
    double Fmom1 = P0.rho*P0.v*P0.v + P0.p;
    double Fnrg1 = P0.v*(P0.rho*(P0.u + 0.5*P0.v*P0.v) + P0.p);

    double Fmass1_Ast = Fmass1*Ast;
    //double Fmom1_Ast = Fmom1*Ast;
    //double Fmom1_Ast = Fmom1*Ast - 0.5*P0.p*dAdx;
    //double Fmom1_Ast = (Fmom1 - 0.5*P0.p)*Ast;
    double Fmom1_Ast = (Fmom1 - P0.p)*Ast;
    double Fnrg1_Ast = Fnrg1*Ast;

    // These increments calculated using reformulate_jac.max. See 18/08 derivation.
    double gamma = gm.gamma(gs).re;
    double u = P0.u;
    double v = P0.v;

    // Double check the signs here, on the RHS specifically
    double dmass = -(2*Fmass1_Ast*u*gamma^^2-Fmass1_Ast*v^^2*gamma
                                       +2*Fmom1_Ast*v*gamma
                                       -2*Fmass1_Ast*u*gamma-2*Fnrg1_Ast*gamma
                                       -3*Fmass1_Ast*v^^2+2*Fnrg1_Ast)
                     /(2*v*(u*gamma^^2-u*gamma-v^^2));
 
    double dmom = -Fmass1_Ast; // FR?

    double dnrg = -(2*Fmass1_Ast*u*v*gamma^^2+Fmass1_Ast*v^^3*gamma
                                        -2*Fmom1_Ast*v^^2*gamma
                                        -6*Fmass1_Ast*u*v*gamma
                                        +2*Fnrg1_Ast*v*gamma
                                        +4*Fmom1_Ast*u*gamma-3*Fmass1_Ast*v^^3
                                        +4*Fmom1_Ast*v^^2-6*Fnrg1_Ast*v)
            /(4*(u*gamma^^2-u*gamma-v^^2));


    double mass0 = P0.rho;
    double mom0  = P0.rho*P0.v;
    double nrg0  = P0.rho*(P0.u + 0.5*P0.v*P0.v);
    writefln("dmass: %e dmom: %e dnrg: %e ", dmass, dmom, dnrg);

    double[3] U1 = [mass0 + dmass, mom0+dmom, nrg0 + dnrg];
    return U1;
}

// We need to check the original matrix by taking the increments from the
// primitive version and multiply them by J. We should get the right had side.

Primitives decode_conserved(double[3] U, ref GasModel gm, ref GasState gs){
    double rho = U[0];
    double v = U[1]/rho;
    double u = U[2]/rho - 0.5*v*v;
    gs.u = u;
    gs.rho = rho;
    gm.update_thermo_from_rhou(gs);

    double T = gs.T.re;
    double p = gs.p.re;
    return Primitives(rho = rho, p = p, v = v, u = u);
}

void check_flux(Primitives P0, Primitives P1, double A0, double dA, double dx, double gamma){

    double A1 = A0 + dA;
    double dAdx = dA/dx;

    double mass0 = P0.rho;
    double mom0  = P0.rho*P0.v;
    double nrg0  = P0.rho*(P0.u + 0.5*P0.v*P0.v);

    double mass1 = P1.rho;
    double mom1  = P1.rho*P1.v;
    double nrg1  = P1.rho*(P1.u + 0.5*P1.v*P1.v);

    double Fmass0 = P0.rho*P0.v;
    double Fmom0 = P0.rho*P0.v*P0.v + P0.p;
    double Fnrg0 = P0.v*(P0.rho*(P0.u + 0.5*P0.v*P0.v) + P0.p);

    double Fmass1 = P1.rho*P1.v;
    double Fmom1 = P1.rho*P1.v*P1.v + P1.p;
    double Fnrg1 = P1.v*(P1.rho*(P1.u + 0.5*P1.v*P1.v) + P1.p);


    double dFAdx_mass = (A1*Fmass1 - A0*Fmass0)/dx;
    double dFAdx_mom = (A1*Fmom1 - A0*Fmom0)/dx - P0.p*dAdx;
    double dFAdx_nrg = (A1*Fnrg1 - A0*Fnrg0)/dx;
    writefln("dFAdx_mass = %e",  dFAdx_mass);
    writefln("dFAdx_mom  = %e",  dFAdx_mom );
    writefln("dFAdx_nrg  =  %e", dFAdx_nrg );

    double v = P0.v;
    double u = P0.u;
    double[3] dFdU0 = [0,1,0];
    double[3] dFdU1 = [(v^^2*(gamma-3))/2,-v*(gamma-3),gamma-1];
    double[3] dFdU2 = [(v*(v^^2*gamma-2*u*gamma-2*v^^2))/2,-(2*v^^2*gamma-2*u*gamma-3*v^^2)/2,v*gamma];

    double[3] dU = [mass1-mass0, mom1-mom0, nrg1-nrg0];

    double Ast = dAdx/(A1/dx + dAdx);
    double[3] rhs = [-Fmass0*Ast,
                     //-Fmom0*Ast,
                     //-Fmom0*Ast + 0.5*P0.p*dAdx,
                     //-(Fmom0-0.5*P0.p)*Ast,
                     -(Fmom0-P0.p)*Ast,
                     -Fnrg0*Ast];

    double[3] lhs = [dFdU0[0]*dU[0] + dFdU0[1]*dU[1] + dFdU0[2]*dU[2],
                     dFdU1[0]*dU[0] + dFdU1[1]*dU[1] + dFdU1[2]*dU[2],
                     dFdU2[0]*dU[0] + dFdU2[1]*dU[1] + dFdU2[2]*dU[2]];

    writefln("lhs-> %15.15e, %15.15e, %15.15e", lhs[0], lhs[1], lhs[2]);
    writefln("rhs-> %15.15e, %15.15e, %15.15e", rhs[0], rhs[1], rhs[2]);
    writefln("diff> %15.15f, %15.15f, %15.15f\n",(lhs[0]-rhs[0])/rhs[0]*100.0,
                                                 (lhs[1]-rhs[1])/rhs[1]*100.0,
                                                 (lhs[2]-rhs[2])/rhs[2]*100.0);

}
