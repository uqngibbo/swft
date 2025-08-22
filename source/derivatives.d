/*
    Module for flowfield derivative calculations in scrf

    @author: Nick Gibbons
*/

module derivatives;

import std.stdio;
import std.math;
import std.format;
import std.string;
import std.conv;

import io;

import gas;
import nm.complex;
import nm.number;

double temp_from_u(GasState gs, GasModel gm){
/*
    This should copy the gasstate, so it won't affect the outer scope.
*/
    // Technically we use rho u in the increment calculation. Should that
    // be reflected here? Maybe it affects the error very slightly.
    gm.update_thermo_from_rhop(gs);
    return gs.T.re;
}

Primitives f_derivative(Primitives P1, Primitives P2, Primitives dP1df, double f, double dA,
                        double A1, double A2, double dx, GasState gs, GasModel gm){
/*
    Note that H does not appear in these expressions, even when included in R
*/
    double rho = P2.rho;
    double p = P2.p;
    double v = P2.v;
    double u = P2.u;

    double rho1 = P1.rho;
    double p1 = P1.p;
    double v1 = P1.v;
    double u1 = P1.u;

    gs.u = u;
    gs.rho = rho;
    gs.p = p;
    gs.T = temp_from_u(gs, gm);

    double diameter = sqrt(4.0*A1/PI);
    double c = 1.0/8.0*PI*diameter*dx;

    double cv = gm.dudT_const_v(gs).re;
    double R = gm.gas_constant(gs).re;
    double A = A2;

    double dr1df = dP1df.rho;
    double dv1df = dP1df.v;
    double dp1df = dP1df.p;
    double du1df = dP1df.u;

    double rhs0 = A1*dr1df*v1+A1*dv1df*rho1;
    double rhs1 = dr1df*(A1*v1^^2-c*f*v1^^2)-c*rho1*v1^^2
                                       +dv1df*(2*A1*rho1*v1-2*c*f*rho1*v1)
                                       +(dA/2+A1)*dp1df;
    double rhs2 = dv1df*(A1*rho1*v1^^2+A1*rho1*(v1^^2/2+u1)+A1*p1)
                +A1*dr1df*v1*(v1^^2/2+u1)+A1*du1df*rho1*v1+A1*dp1df*v1;
    double rhs3 = 0.0;


    double drdf = ((3*R*dA+8*A*cv+2*A*R)*rho*rhs0*v^^2
     +((4*A^^2*cv-2*A*cv*dA)*rho*rhs3+((-4*A*cv)-4*A*R)*rho*rhs1)*v
      +(2*R*dA-4*A*R)*rho*rhs0*u+(4*A*R-2*R*dA)*rho*rhs2+(2*R*dA-4*A*R)*p*rhs0)
     /((2*A*R*dA+4*A^^2*cv)*rho*v^^3+((2*A*R*dA-4*A^^2*R)*rho*u+(2*A*R*dA-4*A^^2*R)*p)
                                    *v);

    double dvdf = -((R*dA+4*A*cv+2*A*R)*rhs0*v^^2+((4*A^^2*cv-2*A*cv*dA)*rhs3
                                      +((-4*A*cv)-4*A*R)*rhs1)
                                      *v+(4*A*R-2*R*dA)*rhs2)
                 /((2*A*R*dA+4*A^^2*cv)*rho*v^^2+(2*A*R*dA-4*A^^2*R)*rho*u+(2*A*R*dA-4*A^^2*R)*p);

    double dpdf = (R*rho*rhs0*v^^3+(2*A*cv*rho*rhs3-2*R*rho*rhs1)*v^^2
                      +(2*R*rho*rhs0*u+2*R*rho*rhs2+2*R*p*rhs0)*v
                      -2*R*rho*rhs1*u-2*R*p*rhs1)
                   /((R*dA+2*A*cv)*rho*v^^2+(R*dA-2*A*R)*rho*u+(R*dA-2*A*R)*p);

    double dudf = (2*A*cv*rho*rhs0*v^^4+((-2*A*cv*dA*rho*rhs3)-4*A*cv*rho*rhs1)*v^^3
                           +(((-3*R*dA)-4*A*cv-2*A*R)*rho*rhs0*u
                            +4*A*cv*rho*rhs2+4*A*cv*p*rhs0)
                            *v^^2
                           +(4*A*R*rho*rhs1*u+(4*A^^2*cv-2*A*cv*dA)*p*rhs3
                                             -4*A*cv*p*rhs1)
                            *v+(4*A*R-2*R*dA)*rho*rhs0*u^^2
                           +((2*R*dA-4*A*R)*rho*rhs2+(4*A*R-2*R*dA)*p*rhs0)*u)
                    /((2*A*R*dA+4*A^^2*cv)*rho^^2*v^^3+((2*A*R*dA-4*A^^2*R)*rho^^2*u
                                 +(2*A*R*dA-4*A^^2*R)*p*rho)
                                 *v);

    return Primitives(rho=drdf, p=dpdf, v=dvdf, u=dudf);
}

Primitives H_derivative(Primitives P1, Primitives P2, Primitives dP1dH, double f,
         double Hdot, double dA,  double A1, double A2, double dx, GasState gs, GasModel gm){
    double rho = P2.rho;
    double p = P2.p;
    double v = P2.v;
    double u = P2.u;

    double rho1 = P1.rho;
    double p1 = P1.p;
    double v1 = P1.v;
    double u1 = P1.u;

    gs.u = u;
    gs.rho = rho;
    gs.p = p;
    gs.T = temp_from_u(gs, gm);

    double diameter = sqrt(4.0*A1/PI);
    double c = 1.0/8.0*PI*diameter*dx;

    double cv = gm.dudT_const_v(gs).re;
    double R = gm.gas_constant(gs).re;
    double A = A2;

    double dr1dH = dP1dH.rho;
    double dv1dH = dP1dH.v;
    double dp1dH = dP1dH.p;
    double du1dH = dP1dH.u;

    double rhs0 = A1*dr1dH*v1+A1*dv1dH*rho1;
    double rhs1 = dr1dH*(A1*v1^^2-c*f*v1^^2)+dv1dH*(2*A1*rho1*v1-2*c*f*rho1*v1)
                                       +(dA/2+A1)*dp1dH;
    double rhs2 = dv1dH*(A1*rho1*v1^^2+A1*rho1*(v1^^2/2+u1)+A1*p1)
                  +A1*dr1dH*v1*(v1^^2/2+u1)+A1*du1dH*rho1*v1+A1*dp1dH*v1+A1*dx;
    double rhs3 = 0.0;

    double drdH = ((3*R*dA+8*A*cv+2*A*R)*rho*rhs0*v^^2
            +((4*A^^2*cv-2*A*cv*dA)*rho*rhs3+((-4*A*cv)-4*A*R)*rho*rhs1)*v
            +(2*R*dA-4*A*R)*rho*rhs0*u+(4*A*R-2*R*dA)*rho*rhs2+(2*R*dA-4*A*R)*p*rhs0)
            /((2*A*R*dA+4*A^^2*cv)*rho*v^^3+((2*A*R*dA-4*A^^2*R)*rho*u+(2*A*R*dA-4*A^^2*R)*p)
                               *v);

    double dvdH = -((R*dA+4*A*cv+2*A*R)*rhs0*v^^2+((4*A^^2*cv-2*A*cv*dA)*rhs3
                                      +((-4*A*cv)-4*A*R)*rhs1)
                                      *v+(4*A*R-2*R*dA)*rhs2)
                /((2*A*R*dA+4*A^^2*cv)*rho*v^^2+(2*A*R*dA-4*A^^2*R)*rho*u
                +(2*A*R*dA-4*A^^2*R)*p);

    double dpdH = (R*rho*rhs0*v^^3+(2*A*cv*rho*rhs3-2*R*rho*rhs1)*v^^2
                      +(2*R*rho*rhs0*u+2*R*rho*rhs2+2*R*p*rhs0)*v
                      -2*R*rho*rhs1*u-2*R*p*rhs1)
                    /((R*dA+2*A*cv)*rho*v^^2+(R*dA-2*A*R)*rho*u+(R*dA-2*A*R)*p);

    double dudH = (2*A*cv*rho*rhs0*v^^4+((-2*A*cv*dA*rho*rhs3)-4*A*cv*rho*rhs1)*v^^3
                           +(((-3*R*dA)-4*A*cv-2*A*R)*rho*rhs0*u
                            +4*A*cv*rho*rhs2+4*A*cv*p*rhs0)
                            *v^^2
                           +(4*A*R*rho*rhs1*u+(4*A^^2*cv-2*A*cv*dA)*p*rhs3
                                             -4*A*cv*p*rhs1)
                            *v+(4*A*R-2*R*dA)*rho*rhs0*u^^2
                           +((2*R*dA-4*A*R)*rho*rhs2+(4*A*R-2*R*dA)*p*rhs0)*u)
                     /((2*A*R*dA+4*A^^2*cv)*rho^^2*v^^3+((2*A*R*dA-4*A^^2*R)*rho^^2*u
                                 +(2*A*R*dA-4*A^^2*R)*p*rho)
                                 *v);
    return Primitives(rho=drdH, p=dpdH, v=dvdH, u=dudH);
}

