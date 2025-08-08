// io.d: reading and writing and storage for edict
// @author: NNG

module io;

import std.stdio;
import std.format;
import std.string;
import std.conv;

import dyaml;
import nm.complex;
import nm.number;

struct SimData {
    double p,T,rho,A,v,M,gamma;
}

void write_solution_to_file(double[] xs, SimData[] simdata, string filename){

    File outfile = File(filename, "wb");
    size_t[1] ibuff; double[1] dbuff; // buffer arrays

    size_t neq = 7;
    size_t N = simdata.length;
    ibuff[0] = neq;  outfile.rawWrite(ibuff);
    ibuff[0] = N;   outfile.rawWrite(ibuff);

    foreach(i; 0 .. N)  { dbuff[0] = xs[i];  outfile.rawWrite(dbuff); }

    foreach(i; 0 .. N)  { dbuff[0] = simdata[i].p;     outfile.rawWrite(dbuff); }
    foreach(i; 0 .. N)  { dbuff[0] = simdata[i].T;     outfile.rawWrite(dbuff); }
    foreach(i; 0 .. N)  { dbuff[0] = simdata[i].rho;   outfile.rawWrite(dbuff); }
    foreach(i; 0 .. N)  { dbuff[0] = simdata[i].A;     outfile.rawWrite(dbuff); }
    foreach(i; 0 .. N)  { dbuff[0] = simdata[i].v;     outfile.rawWrite(dbuff); }
    foreach(i; 0 .. N)  { dbuff[0] = simdata[i].M;     outfile.rawWrite(dbuff); }
    foreach(i; 0 .. N)  { dbuff[0] = simdata[i].gamma; outfile.rawWrite(dbuff); }
    outfile.close();
    return;
}

