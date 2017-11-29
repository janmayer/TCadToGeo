# TCadToGeo
Convert your CAD files to ROOT geometries.

Warning: Highly experimental.

## Summary

Creating geometries from scratch in ROOT is slow and painful. In a lot of cases, the geometries already exist as CAD files, or can be created much easier with high-end CAD Software. This package extents ROOT with an optimized tetraheader shape `TGeoTet` (inspired by Geant4's `G4Tet`)and an import tool using the tetrahedral mesh builder `tetgen`, which provides a `TGeoVolumeAssembly`. 

This project is a ROOT adaptation of [C. Pooles `CADMesh` project](https://github.com/christopherpoole/CADMesh). Use `CADMesh` if you use Geant4 only and not CERN Virtual Monte Calo (strongly recommended). 


## Requirements

- CMake
- [`tetgen`](https://github.com/christopherpoole/tetgen)
- ROOT 6

and for the extended example:

- FairSoft
- FairRoot

for the usage with Virtual Monte Carlo.

CMake needs to find both dependencies. Use environment variables for tetgen, e.g.
```bash
export PATH=/path/to/tetgen/bin:$PATH
export LD_LIBRARY_PATH=/path/to/tetgen/lib:$LD_LIBRARY_PATH
```
and `. /path/to/root/bin/thisroot.sh` for ROOT (mind the space).


## Installation

- Clone the repo and move into it
- Create a build directory and move to it `mkdir build && cd build`
- `cmake .. -DCMAKE_INSTALL_PREFIX=/your/path/TCadToGeo`
- `make && make install`
- `export LD_LIBRARY_PATH=/your/path/TCadToGeo/lib:$LD_LIBRARY_PATH`

The ROOT CLI should then pick up everything it needs automatically,
test e.g. with the basic example `example/geometry/tet.C`
 
`root -l examples/geometry/tet.C`

## Troubleshooting

- `error: unknown type name 'TGeoTet'`
ROOT can't find the library or the associated helper files. Check environment variables.
- `no StreamerInfo found ... version mismatch` 
Must be compiled against the exact ROOT version you want to use it with. Chances are, somebody upgraded the system wide ROOT installation and now you need to recompile. We've all been there.
- Tetgen complains about the file beeing empty (or similar).
Either the path to the CAD file is wrong, or the file is binary. Tetgen can open only ASCII Files. Try converting your file to ASCII-ply.
- Tetgen complains about geometry errors.
Your geometry probably has some problems, e.g. olverlaps, inverted normals, etc. Use a tool to correct the problems. If you have created the mesh yourself, try other setting for the meshing process.
