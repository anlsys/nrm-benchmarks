{ pkgs ? import (builtins.fetchTarball "https://github.com/NixOS/nixpkgs/archive/25.05.tar.gz") {}
}:
pkgs // rec {
  libnrm = pkgs.callPackage ./nix/libnrm.nix { openmp = pkgs.llvmPackages_17.openmp; };
  nrmb = pkgs.callPackage ./nix/nrmb.nix { inherit libnrm; };
}
