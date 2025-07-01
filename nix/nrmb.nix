{ stdenv, autoreconfHook, pkg-config, libnrm, zeromq, czmq, blas }:
stdenv.mkDerivation {
  src = ../.;
  name = "nrm-benchmarks";
  nativeBuildInputs = [ autoreconfHook pkg-config ];
  buildInputs = [ libnrm blas czmq zeromq ];
}
