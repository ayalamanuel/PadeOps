# dealias_field

Filters one full-domain, cell-centred PadeOps `.s3D` budget field using the
production FFT/dealias/inverse-FFT routines. Output retains the input grid,
Fortran array ordering, and raw binary format (no record markers or header).
Use a build with the same real precision and byte order as the original run.
File reads use the standard `decomp_2d_read_one` behavior without an exact
file-size check; the input namelist supplies the grid dimensions.
This is not a reader for cropped fields, profiles, or edge fields with nz+1 levels.

## Build and run

Load the same compiler, MPI, and dependency environment used to build PadeOps.
Reconfigure an existing build to discover the new source, then build the target:

```sh
cmake -S . -B build_opti/gcc
cmake --build build_opti/gcc --target dealias_field -j 4
mkdir -p filtered
mpiexec -n 4 build_opti/gcc/problems/postprocessing_igrid/dealias_field input_dealias.dat
```

Copy and edit the supplied `input_dealias.dat` before running. All dimensions and
lengths in it are examples; use the original simulation's values. Paths are
relative to the launch directory. The input is one explicitly named field;
there is no budget-number inference or automatic filtering of other terms.
The output directory must exist and the output file must not exist. Keeping
the original filename in a separate directory preserves compatibility with
filename-based postprocessing. The program reserves the output before FFT
initialization; a failed run may leave an empty or partial output file.

Use your own build directory if it differs from `build_opti/gcc`. Run the
synthetic integration checks (NumPy required) in the same MPI environment with:

```sh
python tests/test_dealias_field.py build_opti/gcc/problems/postprocessing_igrid/dealias_field
```

`prow=pcol=0` requests automatic decomposition. Otherwise both must be positive
and their product must equal the MPI rank count. As with other PadeOps tools,
the decomposition must fit the real and spectral grids.

## Filter settings

* `PeriodicInZ=.false.`: filter x/y at every z level, without vertical filtering.
  `dealiasType=0` uses the production strict 2/3 horizontal cutoff. The current
  production implementation ignores `dealiasFact` on this path.
* `PeriodicInZ=.true.`: also filter z using the production periodic path; nz
  must be even. For `dealiasType=0`, the current implementation uses
  `dealiasFact*pi/dx` as the cutoff in **all three directions**, even if grid
  spacings differ. This utility deliberately reproduces that behavior. The
  production vertical FFT setup uses exhaustive planning and can take time.
* `dealiasType=1`: production exponential smooth filter, which ignores
  `dealiasFact`. Repeated filtering further attenuates modes; apply once to
  the original file. The default binary cutoff is idempotent to roundoff.

Match the original run's settings. A simulation with dealiasing disabled does
not need this correction. The program does not change signs, normalization,
sample counts, or companion metadata.

Filtering commutes with averaging, so this can repair the missing linear
filter on an exported mean force. It cannot reconstruct correlations that
required filtering instantaneous force before multiplication: filtering
`<u*f>` does not generally recover `<u*G(f)>`. Nor does it repair missing work
components or the historical MKE/TKE split. See
`src/incompressible/pressure_turbine_budget_audit.md` for the audit.
