# Third-party software

TomoSTAR is distributed under the Apache License 2.0 (see `LICENSE` and `NOTICE`). It uses the
following libraries, each under its own permissive licence; the libraries are fetched from NuGet when
the solution is built. Apart from the translation listed under "Translated algorithms", no
third-party source code is copied into this repository.

## Run time

| Package | Version | Licence | Copyright | Used for |
|---|---|---|---|---|
| Silk.NET.OpenCL, Silk.NET.Core | 2.23.0 | MIT | Copyright (c) .NET Foundation and Contributors | OpenCL bindings of the eikonal and LSQR kernels |
| Microsoft.Data.Sqlite, Microsoft.Data.Sqlite.Core | 10.0.12 | MIT | Copyright (c) .NET Foundation and Contributors | Read-only access to the catalogue of QUIVER projects |
| SQLitePCLRaw (bundle_e_sqlite3, core, lib.e_sqlite3, provider.e_sqlite3) | 2.1.12 | Apache-2.0 | Copyright 2014-2024 SourceGear, LLC | SQLite native binding |
| SQLite | (bundled by SQLitePCLRaw) | Public domain | | The database engine |
| Microsoft.Extensions.DependencyModel, Microsoft.DotNet.PlatformAbstractions | 9.0.9, 3.1.6 | MIT | Copyright (c) .NET Foundation and Contributors | Dependencies of the above |

An OpenCL runtime, when present on the machine, is loaded at run time from the system; none is
distributed with TomoSTAR.

## Translated algorithms

| File | Origin | Licence | Used for |
|---|---|---|---|
| `src/TomoStar.Core/Source/OkadaDc3d.cs` | C# translation of the finite-rectangle part of the DC3D subroutine of Y. Okada (1991, revisions to 2002), implementing Okada (1992), Bull. Seismol. Soc. Am. 82(2), 1018-1040, doi:10.1785/BSSA0820021018; Fortran as distributed with okada_wrapper (https://github.com/cutde-org/okada_wrapper) | MIT (okada_wrapper) | Displacement gradients of rectangular dislocations in a half-space, for the Coulomb stress change |

## Tests only

| Package | Version | Licence |
|---|---|---|
| xunit, xunit.runner.visualstudio | 2.9.3, 3.1.4 | Apache-2.0 |
| Microsoft.NET.Test.Sdk | 17.14.1 | MIT |

## Published models

The 1-D velocity models of the built-in library (`tomostar model1d --list`) are numerical values
taken from the cited publications; each entry carries its full reference and DOI, which should be
cited when a model is used.
