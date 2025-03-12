{ stdenv
, autoreconfHook
, bats
, check
, czmq
, git
, hwloc
, jansson
, mpich
, openmp
, papi
, pkg-config
, protobufc
, zeromq
}:
stdenv.mkDerivation {
  src = fetchGit {
    url = "https://github.com/anlsys/libnrm.git";
    ref = "master";
  };
  name = "libnrm";
  prePatch = "echo 0.8.0 > .tarball-version";
  nativeBuildInputs = [ autoreconfHook pkg-config git ];
  buildInputs = [
    bats
    check
    hwloc
    jansson
    mpich
    openmp
    papi
    protobufc
  ];
  propagatedBuildInputs = [
    czmq
    zeromq
  ];
}
