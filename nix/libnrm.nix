{ stdenv, autoreconfHook, pkgconfig, zeromq, czmq, gfortran }:
stdenv.mkDerivation {
  src = fetchGit {
    url = "https://github.com/anlsys/libnrm.git";
    ref = "refs/heads/v0.7.x";
  };
  name = "libnrm";
  nativeBuildInputs = [ autoreconfHook pkgconfig ];
  buildInputs = [ zeromq czmq gfortran ];
}
