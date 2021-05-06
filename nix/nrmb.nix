{ stdenv, autoreconfHook, pkgconfig, libnrm }:
stdenv.mkDerivation {
  src = ../.;
  name = "nrm-benchmarks";
  nativeBuildInputs = [ autoreconfHook pkgconfig libnrm ];
}
