# Efficient Fur and Hair Multiple Scattering Using Volumetric Approximation

This repository contains the Mitsuba-based offline implementation and release
scenes for:

**Efficient Fur and Hair Multiple Scattering Using Volumetric Approximation**

Ruike Hu, Junqiu Zhu, Minghao Lin, Ruian Zhang, Lu Wang, Jie Guo, Yanwen Guo,
and Lingqi Yan

ACM Transactions on Graphics, 45(4), Article 92, 2026

DOI: https://doi.org/10.1145/3811343

The main method is implemented as the `hairvolume` integrator. It approximates
multiple scattering in fur and hair using voxelized volumetric statistics, and
is intended to be used with the release scenes in `scenes/`.

## Repository Layout

- `src/integrators/hairvolume/`: the volumetric multiple-scattering integrator.
- `src/precomputer/`: hair voxelization plugins used to generate `.hairvol`
  and `.vol` files for `hairvolume`.
- `src/bsdfs/` and `src/phase/`: hair scattering models used by the examples.
- `scenes/`: cleaned release scenes, material presets, lighting presets, and
  usage instructions.

This codebase is based on a CMake/vcpkg fork of Mitsuba 0.6. See `LICENSE` for
the license terms inherited from Mitsuba.

## Requirements

The release has been tested on Windows 10/11 with Visual Studio 2022.

Required tools:

- CMake
- Visual Studio 2022 with C++ support
- vcpkg

## Build

Configure with your vcpkg toolchain file:

```powershell
cmake -S . -B cbuild -DCMAKE_TOOLCHAIN_FILE=<path-to-vcpkg>\scripts\buildsystems\vcpkg.cmake
```

Build the renderer and plugins:

```powershell
cmake --build cbuild --config RelWithDebInfo -j
```

The Mitsuba executable will be written to:

```text
cbuild/bin/RelWithDebInfo/mitsuba.exe
```

Plugins are written to:

```text
cbuild/bin/plugins/
```

## Scenes and Data

The release scenes are in `scenes/`. Hair geometry is not included in this
repository; see `scenes/README.md` for download links, expected filenames,
voxelization commands, and rendering examples.

Typical workflow:

1. Download the hair model files and place them in `scenes/`.
2. Generate the hair volume files with `method=hairvoxel`.
3. Render with `method=ours` and a matching `hairvolumePreset`.

## Citation

If you use this code, please cite:

```bibtex
@article{Hu2026EfficientFurHair,
  title={Efficient Fur and Hair Multiple Scattering Using Volumetric Approximation},
  author={Ruike Hu and Junqiu Zhu and Minghao Lin and Ruian Zhang and Lu Wang and Jie Guo and Yanwen Guo and Lingqi Yan},
  journal={ACM Transactions on Graphics},
  volume={45},
  number={4},
  articleno={92},
  year={2026},
  month={7},
  doi={10.1145/3811343},
  url={https://github.com/NJUCG/efficient-fur-hair-ms-vol-approx}
}
```

## Acknowledgement

This project builds on Mitsuba 0.6 and a CMake/vcpkg-based Mitsuba fork. We
thank the Mitsuba authors and contributors for making their renderer available.
