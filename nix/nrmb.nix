{ stdenv, autoreconfHook, pkgconfig, libnrm, zeromq, czmq, blas }:
stdenv.mkDerivation {
  src = ../.;
  name = "nrm-benchmarks";
  nativeBuildInputs = [ autoreconfHook pkgconfig libnrm czmq zeromq blas];
  buildInputs = [ libnrm ];
}
