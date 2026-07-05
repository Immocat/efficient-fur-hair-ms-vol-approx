# Hair Multiple Scattering Scenes

This directory contains the release scenes used to test the `hairvolume`
integrator. The scenes are built around a single parameterized Mitsuba file,
`uber.xml`, which selects the hairstyle, material, lighting setup, and rendering
method through command-line `-D` parameters.

## Hair Data

Hair geometry files are not included in this repository. Place the downloaded
`.hair` files directly in this directory:

```text
scenes/
  wCurly.hair
  afro4.hair
```

Download sources:

- `wCurly.hair`: download the wCurly model from Cem Yuksel's hair model page:
  https://www.cemyuksel.com/research/hairmodels/
- `afro4.hair`: download the Afro hairstyle from:
  https://drive.google.com/file/d/1HJ7gz4qFhuJ5D02x_LayJKWMEux2YSpH/view?usp=sharing

The scene XML files expect these exact filenames. If your downloaded files have
different names, rename them before rendering.

## Scene Options

The main scene file is:

```text
scenes/uber.xml
```

Supported hairstyles:

- `wCurly`
- `afro4`

Supported materials:

- `dEon_blonde`
- `dEon_white`
- `fur2_red`

Supported rendering methods:

- `path`: reference path tracing with the hair BSDF.
- `ours`: the hair volume multiple-scattering approximation.
- `hairvoxel`: voxelizes the hair data using the release voxelizer.
- `hairvoxel_v2`: experimental path for GPU voxelization.

Supported lighting profiles:

- `const`
- `Directional`
- `Env`
- `SpotLight`
- `ThreePointLighting`
- `ThreePointLightingSmall`
- `ThreePointLighting_NoGeo`

## Generating Hair Volumes

The `ours` method requires voxelized hair volume files in this directory:

```text
<hairName>_128_128_128_atten_raw.hairvol
<hairName>_128_128_128_dir.vol
<hairName>_128_128_128_var.vol
```

Generate them with `method=hairvoxel`. From the repository root:

```powershell
.\cbuild\bin\RelWithDebInfo\mitsuba.exe `
  -o scenes\wCurly_voxel_preview.exr `
  -D hairName=wCurly `
  -D method=hairvoxel `
  -D shadingModel=dEon `
  -D material=dEon_blonde `
  -D lightingScene=const `
  -D spp=1 `
  -D filmRes=64 `
  scenes\uber.xml
```

For Afro, change `hairName=wCurly` to `hairName=afro4`.

The voxelizer also writes a small preview render. The preview image is not used
by `ours`; the required outputs are the `.hairvol` and `.vol` files listed
above.

The `hairvoxel_v2` method is an experimental GPU voxelization variant. For
normal offline use, use `method=hairvoxel`.

## Rendering

### Path Tracing Reference

```powershell
.\cbuild\bin\RelWithDebInfo\mitsuba.exe `
  -o scenes\wCurly_path.exr `
  -D hairName=wCurly `
  -D method=path `
  -D shadingModel=dEon `
  -D material=dEon_blonde `
  -D lightingScene=Env `
  -D spp=512 `
  -D filmRes=1024 `
  scenes\uber.xml
```

### Hair Volume Method

Run the voxelizer first, then render with `method=ours` and one of the provided
hair volume presets.

```powershell
.\cbuild\bin\RelWithDebInfo\mitsuba.exe `
  -o scenes\wCurly_blonde_ours.exr `
  -D hairName=wCurly `
  -D method=ours `
  -D hairvolumePreset=wCurly_blonde `
  -D shadingModel=dEon `
  -D shadingModelPhase=dEonPhase `
  -D material=dEon_blonde `
  -D lightingScene=Env `
  -D spp=512 `
  -D filmRes=1024 `
  scenes\uber.xml
```

Available `hairvolumePreset` values:

- `wCurly_blonde`
- `wCurly_fur2_red`
- `afro_blonde`
- `afro_white`

Use the matching material and phase settings:

```text
wCurly_blonde:   material=dEon_blonde, shadingModel=dEon, shadingModelPhase=dEonPhase
wCurly_fur2_red: material=fur2_red,    shadingModel=fur2, shadingModelPhase=fur2Phase
afro_blonde:     material=dEon_blonde, shadingModel=dEon, shadingModelPhase=dEonPhase
afro_white:      material=dEon_white,  shadingModel=dEon, shadingModelPhase=dEonPhase
```

## Useful Parameters

- `spp`: samples per pixel.
- `filmRes`: square image resolution.
- `hairName`: `wCurly` or `afro4`.
- `method`: `path`, `ours`, `hairvoxel`, or `hairvoxel_v2`.
- `lightingScene`: one of the lighting profiles listed above.
- `voxelCount`: voxel resolution used by the voxelization methods.
- `volResX`, `volResY`, `volResZ`: volume resolution read by `ours`.

If you generate volumes with a resolution other than `128`, pass the matching
`volResX`, `volResY`, and `volResZ` values when rendering with `method=ours`.
