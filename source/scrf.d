/*
    scrf: A one-dimensional scramjet modelling code.

    References:

    @author: Nick Gibbons
*/

module scrf;

import std.stdio;
import std.math;
import std.mathspecial;
import std.format;
import std.string;
import std.conv;
import std.datetime.stopwatch : StopWatch;
import core.thread;
import core.sys.posix.signal;

// scrf specific modules
import io;

import gas;
import gas.physical_constants;
import kinetics;
import nm.bbla;
import nm.number;
import nm.complex;

void progress_bar(double x, double L){
    double percent = x/L*100.0;
    int filled = to!int(x/L*60);
    write("\r[");
    foreach(i; 0 .. filled) write('-');
    foreach(i; filled .. 60) write(' ');
    writef("] %3.0f%%\r", percent);
    stdout.flush();
}

void print_state(string name, double v, double M, double A, ref GasState gs, GasModel gm){
    writefln("%4s State: v=%12.12e (m/s) T=%12.12e (K) p=%12.12e Pa", name, v, gs.T.re, gs.p.re);
    writefln("            M=%12.12e rho=%12.12e (g/m3) A=%12.12e (m2) a=%12.12e (m/s) ", M, gs.rho.re*1000.0, A, gs.a.re);
    write("massf: [");
    foreach(isp; 0 .. gm.n_species){
        writef("%s:%3.3e", gm.species_name(isp), gs.massf[isp].re);
        if (isp!=gm.n_species-1) write(", ");
    }
    writeln("]");
}

double temp_from_u(GasState gs, GasModel gm){
	double cv = gm.dudT_const_v(gs).re;
    return gs.u.re/cv;
}

Primitives increment_primitives(double x, double A, double dx, double dA, double Hdot, double f, ref Primitives P0, ref GasModel gm, ref GasState gs){
    double rho = P0.rho;
    double p = P0.p;
    double v = P0.v;
    double u = P0.u;

    double du_chem = 0.0;
    double dp_chem = 0.0;

    double gamma = gm.gamma(gs).re;

	double R = gm.gas_constant(gs).re;
	double cv = gm.dudT_const_v(gs).re;
	double dfdr = R*gs.T.re;
	double dfdu = gs.rho.re*R/cv;

	// Friction factor
	double diameter = sqrt(4.0*A/PI);
	double tau = 1.0/8.0*f*rho*v*v;
	double taupiDdx = tau*PI*diameter*dx;

	// Rayleigh heat addition (I did this derivation at 2330)
	double Qdot = Hdot*A*dx;

	// Compute the accommodation increments using expressions from Maxima.
	// We get slightly different dp_chems to nenzf1d. I wonder why?
	double denom = A*rho^^2*v^^2 - A*dfdr*rho^^2 - A*dfdu*p;
	double drho = -(dA*rho^^3*v^^3+((-rho^^2)-dfdu*rho)*taupiDdx*v-Qdot*dfdu*rho)/(denom*v);
	double dv = -(((rho+dfdu)*taupiDdx-dA*dfdr*rho^^2-dA*dfdu*p)*v+Qdot*dfdu)/denom;
	double dp_gda = ((dfdu*rho*taupiDdx-dA*dfdr*rho^^3-A*dp_chem*rho^^2-dA*dfdu*p*rho)*v^^2
					 +Qdot*dfdu*rho*v+(dfdr*rho^^2+dfdu*p)*taupiDdx+A*dfdr*dp_chem*rho^^2
					 +A*dfdu*dp_chem*p)/denom;
	double du_gda = ((rho*taupiDdx-dA*p*rho)*v^^3+Qdot*rho*v^^2+(p-dfdr*rho)*taupiDdx*v
					 -Qdot*dfdr*rho)/(denom*v);

    return Primitives(rho = rho + drho,
                      p   = p   + dp_gda,
                      v   = v   + dv,
                      u   = u   + du_gda);
}

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

Primitives f_derivative(Primitives P0, Primitives P1, double f, double dA, double A0,
                                           double A1, double dx, GasState gs, GasModel gm){
    double rho = P1.rho;
    double p = P1.p;
    double v = P1.v;
    double u = P1.u;

    double rho0 = P0.rho;
    double p0 = P0.p;
    double v0 = P0.v;
    double u0 = P0.u;

    gs.u = u;
    gs.rho = rho;
    gs.p = p;
    gs.T = temp_from_u(gs, gm);

	double diameter = sqrt(4.0*A0/PI);
	double c = 1.0/8.0*PI*diameter*dx;

	double cv = gm.dudT_const_v(gs).re;
	double R = gm.gas_constant(gs).re;
    double A = A1;

    double mass=rho*v*A - rho0*v0*A0;
    double mom=rho*v*v*A + p*A - rho0*v0*v0*A0 - p0*A0 -(p+p0)/2*dA + c*f*rho0*v0*v0;
    double nrg=rho*v*A*(u+v*v/2) + p*A*v - rho0*v0*A0*(u0+v0*v0/2) - p0*A0*v0;
    double eos=p-rho*R*u/cv;

    writefln("mass %e", mass);
    writefln("mom %e", mom);
    writefln("nrg %e (%e vs. %e)", nrg, (rho*v*A*(u+v*v/2) + p*A*v), rho0*v0*A0*(u0+v0*v0/2) + p0*A0*v0);
    writefln("eos %e", eos);

    double D = (R*dA+2*A*cv)*rho*v^^2+(R*dA-2*A*R)*rho*u+(R*dA-2*A*R)*p;
    double drdf = ((2*c*cv+2*R*c)*rho*rho0*v0^^2)/D;
    double dvdf = -((2*c*cv+2*R*c)*rho0*v*v0^^2)/D;
    double dpdf = ((2*R*c*rho*rho0*v^^2+2*R*c*rho*rho0*u+2*R*c*p*rho0)*v0^^2)/D;
    double dudf = ((2*c*cv*rho*rho0*v^^2-2*R*c*rho*rho0*u+2*c*cv*p*rho0)*v0^^2)
                  /((R*dA+2*A*cv)*rho^^2*v^^2+(R*dA-2*A*R)*rho^^2*u+(R*dA-2*A*R)*p*rho);

    return Primitives(rho=drdf, p=dpdf, v=dvdf, u=dudf);
}


int main(string[] args)
{
    int exitFlag = 0;
    writefln("scrf: A Q1D Flow Analysis Tool");

    string config_file_name = "scrf.yaml";
    if (args.length>1) config_file_name = args[1];
    Config cfg = Config(config_file_name);

    GasModel gm = init_gas_model(cfg.gas_file_name);
    GasState gs = GasState(gm);
    GasState gs2= GasState(gm);

    gs.p = cfg.p0;
    gs.T = cfg.T0;
    foreach(sp,mf; cfg.Y0) gs.massf[gm.species_index(sp)] = mf;
    gm.update_thermo_from_pT(gs);
    gm.update_sound_speed(gs);

    double v = cfg.v0;
	double M = v/gs.a.re;

    double x = 0.0;
    double L = cfg.L;
    double rs = cfg.rs;
    double rf = cfg.rf;
    double dt = cfg.dt;
    double As = PI*rs*rs;
    double f = cfg.f;
    double Hdot = cfg.Hdot; // Volumetric heat addition rate W/m3

    print_state("Init", v, M, As, gs, gm);

    SimData[] simdata;
    double[] xs;
    size_t nreserve = to!size_t(L/(v*dt))*2;
    writefln("Timestep %e, Reserving space for %d simdatas", dt, nreserve);
    simdata.reserve(nreserve);
    xs.reserve(nreserve);

    double gamma = gm.gamma(gs).re;

    Primitives P0 = Primitives(gs.rho.re, gs.p.re, v, gs.u.re);
    double[3] U0 = [P0.rho,
                    P0.rho*v,
                    P0.rho*(P0.u + 0.5*v*v)];
    simdata ~= SimData(P0.p, gs.T.re, P0.rho, As, P0.v, M, gamma);
    xs ~= x;

    size_t iter = 0; 
    bool last_step = false;
    writeln("Running...");
    while (x<=L) {
        double r = (x-0.0)/L*(rf-rs) + rs;
        double A = PI*r*r;

        double dx = v*dt;
        if (x+dx>=L) {
            last_step = true;
            dx = L-x;
        }
        double x1 = x + dx;
        double r1 = (x1-0.0)/L*(rf-rs) + rs;
        double A1 = PI*r1*r1;
        double dA = A1-A;

        double eps = 1e-7;
        Primitives P1 = increment_primitives(x, A, dx, dA, Hdot, f,     P0, gm, gs);
        Primitives P1c= increment_primitives(x, A, dx, dA, Hdot, f+eps, P0, gm, gs);
        Primitives dPdf_f = (P1c-P1)/eps;
        Primitives dPdf_a = f_derivative(P0, P1, f, dA, A, A1, dx, gs2, gm);
        writefln("dPdf_f-> %s", dPdf_f);
        writefln("dPdf_a-> %s", dPdf_a);
        writefln("  diff-> %s", (dPdf_a-dPdf_f)/dPdf_f*100.0);
        if (M!=0.123) return 1;


        //double[3] U1  = increment_conserved( x, A, dx, dA, Hdot, f, P0, gm, gs);
        //Primitives P1c =  decode_conserved(U1, gm, gs);
        //Primitives P1 = P1p;

        // Add the increments
        x = x + dx;
        gs.rho = P1.rho;
        gs.u = P1.u;
        gs.p = P1.p;
        gs.T = temp_from_u(gs, gm);
        //gs.T = gs.T.re;
        gm.update_sound_speed(gs);
        gamma = gm.gamma(gs).re;

        M = P1.v/gs.a.re;
        simdata ~= SimData(P1.p, gs.T.re, P1.rho, A, P1.v, M, gamma);
        xs ~= x;
        P0 = P1;

        iter += 1;
        if ((iter%20==0)||(last_step)){
            progress_bar(x, L);
        }
        if (last_step) break;
    }
    writeln("");
    writefln("Done in %d iters", iter);
    print_state(" End", v, M, As, gs, gm);
    string output_file_name = format("%s.bin", config_file_name.chomp(".yaml"));
    writefln("Writing solution to file %s...", output_file_name);
    write_solution_to_file(xs, simdata, output_file_name);

    return exitFlag;
} // end main()
