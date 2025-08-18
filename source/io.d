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

struct Config {
    string gas_file_name;
    string reaction_file_name;
    double L;
    double rs;
    double rf;
    double dt;
    double f;
    double Hdot;
    double T0;
    double p0;
    double v0;
    double[string] Y0;

    this(string filename){
        Node data = dyaml.Loader.fromFile(filename).load();

        gas_file_name      = data["gas_file_name"].as!string;
        reaction_file_name = data["reaction_file_name"].as!string;

        L    = to!double(data["L"].as!string);
        rs   = to!double(data["rs"].as!string);
        rf   = to!double(data["rf"].as!string);
        dt   = to!double(data["dt"].as!string);
        f    = to!double(data["f"].as!string);
        Hdot = to!double(data["Hdot"].as!string);
        T0   = to!double(data["T0"].as!string);
        p0   = to!double(data["p0"].as!string);
        v0   = to!double(data["v0"].as!string);

        foreach(Node nd; data["Y0"].mappingKeys) {
            string s = nd.as!string;
            double Ys= to!double(data["Y0"][s].as!string);
            Y0[s] = Ys;
        }

        // Optional parameters
        return;
    }
}

struct Primitives {
    double rho,p,v,u;

    Primitives opBinary(string op)(in Primitives rhs) if ( op == "+" || op == "-" || op == "*" || op == "/"  ) {
		return mixin("Primitives(rho"~op~"rhs.rho, p"~op~"rhs.p, v"~op~"rhs.v, u"~op~"rhs.u)");
    }

    Primitives opBinary(string op)(in double rhs) if ( op == "+" || op == "-" || op == "*" || op == "/"  ) {
		return mixin("Primitives(rho"~op~"rhs, p"~op~"rhs, v"~op~"rhs, u"~op~"rhs)");
    }

    string toString() const
    {
        char[] repr;
        repr ~= "Primitives(";
        repr ~= format("rho=%12.12e, ", rho);
        repr ~= format("p=%12.12e, ", p);
        repr ~= format("v=%12.12e, ", v);
        repr ~= format("u=%12.12e)", u);
        return to!string(repr);
    }

}

// Maybe there should be another one called Primitives or something??
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

