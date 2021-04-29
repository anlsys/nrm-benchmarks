{ pkgs ? import (builtins.fetchTarball "https://github.com/NixOS/nixpkgs/archive/20.09.tar.gz") {}
}:
let
  callPackage = pkgs.lib.callPackageWith pkgs;
in pkgs // {
  nrmb = callPackage ./nix/nrmb.nix {};
}
