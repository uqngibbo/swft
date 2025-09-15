/*
    swft: A one-dimensional scramjet modelling code.

    References:

    Todo List:
     - Python bindings
     - Reactions
     - Actually try the scramjet flow
     - Nonsymbolic derivatives that are a bit less messed up.

    @author: Nick Gibbons
*/

module swft;

import std.stdio;
import std.math;
import std.format;
import std.string;
import std.conv;
import std.typecons: Tuple;

// swft specific modules
import io;
import derivatives;

import gas;
import gas.physical_constants;
import kinetics;
import nm.bbla;
import nm.number;
import nm.complex;

void progress_bar(double x, double x0, double xf){
    double percent = (x-x0)/(xf-x0)*100.0;
    int filled = to!int((x-x0)/(xf-x0)*60);
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

double interpolate(double[] xi, double[] Ai, double x){
/*
    Non uniform linear interpolation using a binary search.

    @author: Nick Gibbons
*/
    size_t n = xi.length;

    size_t l = 0;
    size_t m = n/2;
    size_t u = n-1;
    double xl = xi[l]; double xm = xi[m]; double xu = xi[u];
    if (x<=xi[0]) return Ai[0];
    if (x>=xi[$-1]) return Ai[$-1];

    while (true) {
        if ((u-l)==1) break;
        if (x>xm) {
            l = m;
        } else if (x<=m) {
            u = m;
        } else {
            throw new Error("Bad binary search logic sorry.");
        }
        m = (l+u)/2;
        xl = xi[l]; xm = xi[m]; xu = xi[u];
    }

    double dAdx = (Ai[u]-Ai[l])/(xi[u]-xi[l]);
    double A = dAdx*(x-xi[l]) + Ai[l];
    return A;
}

Tuple!(size_t, size_t) get_bracketing_indices(double[] xi, double x){
    size_t n = xi.length;

    size_t l = 0;
    size_t m = n/2;
    size_t u = n-1;
    double xl = xi[l]; double xm = xi[m]; double xu = xi[u];
    if (x<xi[0]) throw new Error("x is more negative than xi[0]");
    if (x>xi[$-1]) throw new Error("x is more positive than xi[$-1]");

    while (true) {
        if ((u-l)==1) break;
        if (x>xm) {
            l = m;
        } else if (x<=m) {
            u = m;
        } else {
            throw new Error("Bad binary search logic sorry.");
        }
        m = (l+u)/2;
        xl = xi[l]; xm = xi[m]; xu = xi[u];
    }
    return Tuple!(size_t, size_t)(l, u);
}

double calc_dfdfi(double[] xf, double xj, size_t i){

    if (xf.length<2) throw new Error("xf has not enough entries in it");
    if (i>=xf.length) throw new Error(format("Bad node index %d given for xf of size %d",
                                             i, xf.length));

    Tuple!(size_t, size_t) indexes = get_bracketing_indices(xf, xj);
    size_t l = indexes[0];
    size_t u = indexes[1];

    if (i==l){
        // Is it always safe to +1 here? I think it is
        return 1.0 - (xj-xf[i])/(xf[i+1] - xf[i]);
    } else if (i==u){
        return (xj-xf[i-1])/(xf[i] - xf[i-1]);
    } else {
        return 0.0;
    }
}


Primitives increment_primitives(double x, double A, double dx, double dA, double Hdot, double CH, double uw, double f, ref Primitives P0, ref GasModel gm, ref GasState gs){
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
    double r = sqrt(0.71); // Recovery Factor

    // Friction factor
    double diameter = sqrt(4.0*A/PI);
    double tau = 1.0/8.0*f*rho*v*v;
    double taupiDdx = tau*PI*diameter*dx;

    // Rayleigh heat addition (I did this derivation at 2330)
    //double Qdot = Hdot*A*dx;
    double Qdot = Hdot*A*dx + rho*v*CH*(u + r*v*v/2.0 - uw)*dx;

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

    return Primitives(rho + drho,
                      p   + dp_gda,
                      v   + dv,
                      u   + du_gda);
}


int main(string[] args)
{
    int exitFlag = 0;
    writefln("swft: A Q1D Flow Analysis Tool");

    string config_file_name = "swft.yaml";
    if (args.length>1) config_file_name = args[1];
    Config cfg = Config(config_file_name);

    GasModel gm = init_gas_model(cfg.gas_file_name);
    GasState gs = GasState(gm);

    gs.p = cfg.p0;
    gs.T = cfg.T0;
    foreach(sp,mf; cfg.Y0) gs.massf[gm.species_index(sp)] = mf;
    gm.update_thermo_from_pT(gs);
    gm.update_sound_speed(gs);

    double v0 = cfg.v0;
    double M = v0/gs.a.re;
    double cv = gm.dudT_const_v(gs).re;
    double R = gm.gas_constant(gs).re;

    double x0 = cfg.xi[0];
    double xf = cfg.xi[$-1];
    double dt = cfg.dt;
    double As = interpolate(cfg.xi, cfg.Ai, x0);
    double Hdot = cfg.Hdot; // Volumetric heat addition rate W/m3
    double uw = cfg.Tw*cv;
    double CH = cfg.CH;

    print_state("Init", v0, M, As, gs, gm);

    SimData[][] fderivs;
    SimData[] Hderivs;

    double[] xs;
    SimData[] simdata;
    size_t nreserve = to!size_t((xf-x0)/(v0*dt))*2;
    writefln("Timestep %e, Reserving space for %d simdatas", dt, nreserve);
    simdata.reserve(nreserve);
    xs.reserve(nreserve);

    double gamma = gm.gamma(gs).re;

    Primitives P0 = Primitives(gs.rho.re, gs.p.re, v0, gs.u.re);
    
    Primitives[][] dPkdfj; // dUdf at each point (x)
    Primitives[] dPdfs;     // dUdf at each f schedule control point
    fderivs.length = cfg.f.length;
    dPdfs.length = cfg.f.length;

    Primitives dPdH0 = Primitives(0.0, 0.0, 0.0, 0.0);
    Primitives dPdH;

    simdata ~= SimData(P0.p, gs.T.re, P0.rho, As, P0.v, M, gamma);
    xs ~= x0;
    if (cfg.calc_derivatives) {
        foreach(i; 0 .. cfg.f.length){
            fderivs[i] ~= SimData(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        Hderivs~= SimData(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        dPkdfj.length=1;
    }

    size_t iter = 0; 
    bool last_step = false;
    double x = x0;
    writeln("Running...");
    while (x<=xf) {
        double A = interpolate(cfg.xi, cfg.Ai, x);
        double f = interpolate(cfg.xf, cfg.f, x);

        double dx = v0*dt; // This needs to be constant to get the right derivatives...
        if (x+dx>=xf) {
            last_step = true;
            dx = xf-x;
        }
        double x1 = x + dx;
        double A1 = interpolate(cfg.xi, cfg.Ai, x1);
        double dA = A1-A;

        Primitives P1 = increment_primitives(x, A, dx, dA, Hdot, CH, uw, f, P0, gm, gs);
        if (cfg.calc_derivatives) {

            Primitives[] dP1dfj;
            f_derivative_nonuniform(P0, P1, dPkdfj[iter],f, dA, dP1dfj, A, A1, dx, gs, gm);
            foreach(i; 0 .. cfg.f.length){
                dPdfs[i] = Primitives(0.0, 0.0, 0.0, 0.0);
                foreach(j, dP1df_x; dP1dfj){
                    double xj = xs[j];
                    double dfdfi =  calc_dfdfi(cfg.xf, xj, i);
                    dPdfs[i] += dP1df_x*dfdfi;
                }
            }
            dPkdfj ~= dP1dfj;
            dPdH = H_derivative(P0, P1, dPdH0, f, Hdot, dA, A, A1, dx, gs, gm);
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
            // TODO: Get this outta here.... I think the derivatives deserve their own
            // datastructure.
            foreach(i; 0 .. cfg.f.length){
                double dTdfi = dPdfs[i].u/cv;
                double dadfi = 0.5*sqrt(gamma*R/gs.T.re)*dTdfi;
                double dMdfi= (gs.a.re*dPdfs[i].v - P1.v*dadfi)/gs.a.re/gs.a.re;
                fderivs[i] ~= SimData(dPdfs[i].p, dTdfi, dPdfs[i].rho, 0.0, dPdfs[i].v, dMdfi, 0.0);
            }

            double dTdH = dPdH.u/cv;
            double dadH = 0.5*sqrt(gamma*R/gs.T.re)*dTdH;
            double dMdH = (gs.a.re*dPdH.v - P1.v*dadH)/gs.a.re/gs.a.re;
            Hderivs ~= SimData(dPdH.p, dTdH, dPdH.rho, 0.0, dPdH.v, dMdH, 0.0);
            dPdH0 = dPdH;
        }

        iter += 1;
        if ((iter%20==0)||(last_step)){
            progress_bar(x, x0, xf);
        }
        if (last_step) break;
    }
    writeln("");
    writefln("Done in %d iters", iter);
    print_state(" End", P0.v, M, As, gs, gm);

    string output_file_name = format("%s.bin", config_file_name.chomp(".yaml"));
    writefln("Writing solution to file %s...", output_file_name);
    write_solution_to_file(xs, simdata, output_file_name);

    if (cfg.calc_derivatives) {
        string derivs_file_name;
        foreach(i; 0 .. cfg.f.length){
            derivs_file_name = format("fderivs_%04d-%s.bin", i, config_file_name.chomp(".yaml"));
            writefln("Writing derivs to file %s...", derivs_file_name);
            write_solution_to_file(xs, fderivs[i], derivs_file_name);
        }

        derivs_file_name = format("Hderivs-%s.bin", config_file_name.chomp(".yaml"));
        writefln("Writing derivs to file %s...", derivs_file_name);
        write_solution_to_file(xs, Hderivs, derivs_file_name);
    }

    return exitFlag;
} // end main()
