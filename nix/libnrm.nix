{ stdenv, autoreconfHook, pkgconfig, zeromq, czmq, gfortran }:
stdenv.mkDerivation {
  src = fetchGit {
    url = "https://github.com/anlsys/libnrm.git";
  };
  name = "libnrm";
  nativeBuildInputs = [ autoreconfHook pkgconfig ];
  buildInputs = [ zeromq czmq gfortran ];
}
