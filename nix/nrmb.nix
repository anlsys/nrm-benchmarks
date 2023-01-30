{ stdenv, autoreconfHook, pkgconfig, libnrm, zeromq }:
stdenv.mkDerivation {
  src = ../.;
  name = "nrm-benchmarks";
  nativeBuildInputs = [ autoreconfHook pkgconfig libnrm zeromq];
}
