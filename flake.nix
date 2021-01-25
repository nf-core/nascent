{
  description = "A flake for Emulating virtualenv with nix-shell";
  # Provides abstraction to boiler-code when specifying multi-platform outputs.
  inputs.flake-utils.url = "github:numtide/flake-utils";
  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let pkgs = nixpkgs.legacyPackages.${system};
      in {
        devShell = pkgs.mkShell {
          name = "python";
          buildInputs = with pkgs; [
            python37
            python37Packages.pip
            python37Packages.setuptools
            python37Packages.virtualenvwrapper
          ];
          shellHook = ''
            export SOURCE_DATE_EPOCH=315532800
            alias pip="PIP_PREFIX='$HOME/.pip_packages' \pip"
            export PYTHONPATH="$HOME/.pip_packages/lib/python3.7/site-packages:$PYTHONPATH"
          '';
        };
      });
}
