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


int main(string[] args)
{
    int exitFlag = 0;

    writefln("Hello world!");
    GasModel gm = init_gas_model("gm.lua");
    GasState gs = GasState(gm);
    SimData[] simdata;
    double[] xs;

    gs.p = 968.0;
    gs.T = 361.0;
    gm.update_thermo_from_pT(gs);
    gm.update_sound_speed(gs);

    double v = 3623.0;
	double M = v/gs.a.re;
    writefln("M: %s", M);
    writefln("gs: %s", gs);

    double x = 0.0;
    double L = 1.0;
    double rs = 0.05;
    double rf = 0.10;
    double dt = 1e-6;
    double As = PI*rs*rs;

    size_t nreserve = to!size_t(L/(v*dt))*2;
    writefln("Reserving space for %d simdatas", nreserve);
    simdata.reserve(nreserve);
    xs.reserve(nreserve);

    double rho = gs.rho.re;
    double p = gs.p.re;
    double T = gs.T.re;
    double du_chem = 0.0;
    double dp_chem = 0.0;
    double gamma = gm.gamma(gs).re;

    simdata ~= SimData(p, T, rho, As, v, M, gamma);
    xs ~= x;

    size_t iter = 0; 
    bool last_step = false;
    write("Running"); stdout.flush();
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
        // THis is according to past nick who did these at 11pm and should be checked!
        double R = gm.gas_constant(gs).re;
        double cv = gm.dudT_const_v(gs).re;
        double dfdr = R*gs.T.re;
        double dfdu = gs.rho.re*R/cv;

        // Compute the accommodation increments using expressions from Maxima.
        double denom = A*(rho*rho*v*v - dfdr*rho*rho - dfdu*p);
        double drho = (A*(dp_chem - du_chem*dfdu)*rho*rho - dA*rho^^3*v*v) / denom;
        double dv = (A*(du_chem*dfdu - dp_chem)*rho +
                     dA*(dfdr*rho*rho + dfdu*p))*v / denom;
        double dp_gda = -((dA*dfdr*rho^^3 + A*dfdu*du_chem*rho*rho + dA*dfdu*p*rho)*v*v
                          - A*dfdr*dp_chem*rho*rho - A*dfdu*dp_chem*p) / denom;
        double du_gda = -(A*(du_chem*rho*rho*v*v - du_chem*dfdr*rho*rho - dp_chem*p)
                          + dA*p*rho*v*v) / denom;

        // Add the increments
        x = x + dx;
        v = v + dv;
        gs.rho += drho;
        gs.u += du_gda;
        gm.update_thermo_from_rhou(gs);
        gm.update_sound_speed(gs);

        // Reencode back to stack local variables
		rho = gs.rho.re;
		p = gs.p.re;
        T = gs.T.re;
        gamma = gm.gamma(gs).re;
		du_chem = 0.0;
		dp_chem = 0.0;
        M = v/gs.a.re;
        simdata ~= SimData(p, T, rho, A, v, M, gamma);
        xs ~= x;

        iter += 1;
		if (iter%1==0){ // Progress bar maybe
			write(".");
            stdout.flush();
		}
        if (last_step) break;
    }
    writeln("");
    writefln("Done in %d iters: x=%f v=%f M=%f", iter, x, v, M);
    writefln("Out gs: %s", gs);
    write_solution_to_file(xs, simdata, "solution.bin");

    return exitFlag;
} // end main()
