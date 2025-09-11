{ stdenv, autoreconfHook, pkg-config, libnrm, zeromq, czmq, blas, mpich }:
stdenv.mkDerivation {
  src = ../.;
  name = "nrm-benchmarks";
  nativeBuildInputs = [ autoreconfHook pkg-config ];
  buildInputs = [ libnrm blas czmq zeromq mpich ];
}
