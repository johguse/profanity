with import <nixpkgs> { };
# Make sure you have CUDA support enabled
# See: https://nixos.wiki/wiki/Nvidia
stdenv.mkDerivation {
  name = "dev";
  buildInputs = [
    gcc
    make
    opencl-headers
    cudatoolkit
  ];
}
