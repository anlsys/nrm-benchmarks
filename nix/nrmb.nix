{ stdenv, autoreconfHook, pkgconfig }:
stdenv.mkDerivation {
  src = ../.;
  name = "nrm-benchmarks";
  nativeBuildInputs = [ autoreconfHook pkgconfig ];
}
