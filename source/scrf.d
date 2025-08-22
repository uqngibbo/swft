/*
    scrf: A one-dimensional scramjet modelling code.

    References:

    Todo List:
     - Python bindings
     - Thermally Perfect Gas
     - Reactions
     - Multiple f and H values
     - Actually try the scramjet flow
     - Nonsymbolic derivatives that are a bit less messed up.

    @author: Nick Gibbons
*/

module scrf;

import std.stdio;
import std.math;
import std.format;
import std.string;
import std.conv;

// scrf specific modules
import io;
import derivatives;

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
    double cv = gm.dudT_const_v(gs).re;
    double R = gm.gas_constant(gs).re;

    double x = 0.0;
    double L = cfg.L;
    double rs = cfg.rs;
    double rf = cfg.rf;
    double dt = cfg.dt;
    double As = PI*rs*rs;
    double f = cfg.f;
    double Hdot = cfg.Hdot; // Volumetric heat addition rate W/m3

    print_state("Init", v, M, As, gs, gm);

    SimData[] fderivs;
    SimData[] Hderivs;
    SimData[] simdata;
    double[] xs;
    size_t nreserve = to!size_t(L/(v*dt))*2;
    writefln("Timestep %e, Reserving space for %d simdatas", dt, nreserve);
    simdata.reserve(nreserve);
    xs.reserve(nreserve);

    double gamma = gm.gamma(gs).re;

    Primitives P0 = Primitives(gs.rho.re, gs.p.re, v, gs.u.re);
    Primitives dPdf0 = Primitives(0.0, 0.0, 0.0, 0.0);
    Primitives dPdf;
    Primitives dPdH0 = Primitives(0.0, 0.0, 0.0, 0.0);
    Primitives dPdH;


    simdata ~= SimData(P0.p, gs.T.re, P0.rho, As, P0.v, M, gamma);
    xs ~= x;
    if (cfg.calc_derivatives) {
        fderivs~= SimData(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        Hderivs~= SimData(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

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

        Primitives P1 = increment_primitives(x, A, dx, dA, Hdot, f, P0, gm, gs);
        if (cfg.calc_derivatives) {
            dPdf = f_derivative(P0, P1, dPdf0, f, dA, A, A1, dx, gs2, gm);
            dPdH = H_derivative(P0, P1, dPdH0, f, Hdot, dA, A, A1, dx, gs2, gm);
        }

        // Add the increments
        x = x + dx;
        gs.rho = P1.rho;
        gs.u = P1.u;
        gs.p = P1.p;
        gs.T = temp_from_u(gs, gm);
        gm.update_sound_speed(gs);
        gamma = gm.gamma(gs).re;
        cv = gm.dudT_const_v(gs).re;
        R = gm.gas_constant(gs).re;

        M = P1.v/gs.a.re;
        simdata ~= SimData(P1.p, gs.T.re, P1.rho, A, P1.v, M, gamma);
        xs ~= x;
        P0 = P1;

        if (cfg.calc_derivatives){
            // Custom constructor that autodiffs the T and M maybe?
            double dTdf = dPdf.u/cv;
            double dadf = 0.5*sqrt(gamma*R/gs.T.re)*dTdf;
            double dMdf = (gs.a.re*dPdf.v - v*dadf)/gs.a.re/gs.a.re;
            fderivs ~= SimData(dPdf.p, dTdf, dPdf.rho, 0.0, dPdf.v, dMdf, 0.0);
            double dTdH = dPdH.u/cv;
            double dadH = 0.5*sqrt(gamma*R/gs.T.re)*dTdH;
            double dMdH = (gs.a.re*dPdH.v - v*dadH)/gs.a.re/gs.a.re;
            Hderivs ~= SimData(dPdH.p, dTdH, dPdH.rho, 0.0, dPdH.v, dMdH, 0.0);
            dPdf0 = dPdf;
            dPdH0 = dPdH;
        }

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

    if (cfg.calc_derivatives) {
        string derivs_file_name = format("fderivs-%s.bin", config_file_name.chomp(".yaml"));
        writefln("Writing derivs to file %s...", derivs_file_name);
        write_solution_to_file(xs, fderivs, derivs_file_name);

        derivs_file_name = format("Hderivs-%s.bin", config_file_name.chomp(".yaml"));
        writefln("Writing derivs to file %s...", derivs_file_name);
        write_solution_to_file(xs, Hderivs, derivs_file_name);
    }

    return exitFlag;
} // end main()
