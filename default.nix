{ pkgs ? import (builtins.fetchTarball "https://github.com/NixOS/nixpkgs/archive/24.05.tar.gz") {}
}:
pkgs // rec {
  libnrm = pkgs.callPackage ./nix/libnrm.nix { openmp = pkgs.llvmPackages_15.openmp; };
  nrmb = pkgs.callPackage ./nix/nrmb.nix { inherit libnrm; };
}
