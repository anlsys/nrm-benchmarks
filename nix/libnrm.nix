{ stdenv, autoreconfHook, pkgconfig, zeromq, gfortran }:
stdenv.mkDerivation {
  src = fetchGit {
    url = "https://github.com/anlsys/libnrm.git";
  };
  name = "libnrm";
  nativeBuildInputs = [ autoreconfHook pkgconfig gfortran ];
  propagatedBuildInputs = [ zeromq ];
}
