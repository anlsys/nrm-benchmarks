{ stdenv, autoreconfHook, pkgconfig, libnrm, zeromq, czmq, blas }:
stdenv.mkDerivation {
  src = ../.;
  name = "nrm-benchmarks";
  nativeBuildInputs = [ autoreconfHook pkgconfig ];
  buildInputs = [ libnrm blas czmq zeromq ];
}
