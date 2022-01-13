{ stdenv, autoreconfHook, pkgconfig, libnrm, zeromq, blas }:
stdenv.mkDerivation {
  src = ../.;
  name = "nrm-benchmarks";
  nativeBuildInputs = [ autoreconfHook pkgconfig libnrm zeromq blas];
  buildInputs = [ libnrm ];
}
