# Validation of TomoSTAR

Every method of TomoSTAR is compared with an independent reference: a published, peer-reviewed code
run on the same input, or a published result of a peer-reviewed study. The input is real data: the
2016-2017 Amatrice-Visso-Norcia sequence of central Italy (the example of `examples/norcia2016`)
and published reference models. A comparison passes when its error, defined case by case below, is
below 5 %. Comparisons whose reference used other input (another catalogue's picks, another
criterion, another data set) are reported for information: their difference measures the input,
not the code, and each is shown next to the same comparison made for the reference code itself.

The scripts are in `validation/`; each writes its results to `validation/results/<case>.json`, and
`validation/report.py` copies them into the table below.

## Data

- Earthquakes and analyst picks: the INGV bulletin (ISIDe Working Group, 2007) for the M >= 2.5
  events of 24 August 2016 to 1 March 2017 inside the example's area (3569 events), downloaded by
  `examples/norcia2016/download.sh` from the INGV FDSN services.
- Stations and responses: StationXML of the networks IV (INGV, 2005), MN (MedNet Project Partner
  Institutions, 1990), 3A (INGV, CNR-IGAG and CNR-IDPA, 2018), XO (EMERSITO Working Group, 2018),
  8P (Marzorati et al., 2023), VM (INGV, 2023) and 5M (Wölbern et al., 2020), from the ORFEUS
  federator.
- Waveforms: the vertical and horizontal records of the IV stations for the 295 M >= 3.5 events.
- Published catalogues: CAT1 (absolute locations, NonLinLoc) and CAT2 (double-difference
  relocations, HypoDD; Michele et al., 2020) of Chiaraluce et al. (2022), from the British
  Geological Survey (Chiaraluce et al., 2022b).
- Published velocity models: the 1-D and 3-D models of Carannante et al. (2013) from the authors'
  data set (Carannante et al., 2025); ak135 (Kennett et al., 1995).

## Reference software

| Software | Version | Reference |
|---|---|---|
| NonLinLoc | 7.1.05, commit 7a62fa4 of github.com/alomax/NonLinLoc | Lomax et al. (2000) |
| HypoDD | 2.1b, commit fb40407 of github.com/fwaldhauser/HypoDD | Waldhauser and Ellsworth (2000); Waldhauser (2001) |
| AttenTIon | commit e34f1d8 of github.com/swei-seismo/AttenTIon | Stachnik et al. (2004); Wei and Wiens (2018) |
| PyKonal | 0.4.1 | White et al. (2020) |
| TauP (through ObsPy) | ObsPy 1.5.1 | Crotwell et al. (1999); Krischer et al. (2015) |
| ObsPy | 1.5.1 | Beyreuther et al. (2010); Krischer et al. (2015) |
| SciPy | 1.17.1 | Virtanen et al. (2020) |

NonLinLoc and HypoDD were compiled from source (gcc and gfortran); the Python packages ran on
Python 3.11 with NumPy 2.4.

## Results

<!-- results begin -->
<!-- results end -->

## Notes on the comparisons

<!-- notes begin -->
<!-- notes end -->

## Reproducing

```sh
# the example's data and results
cd examples/norcia2016 && bash download.sh && tomostar run norcia.tomo && cd ../..
# the reference data and codes (see the docstring of each script)
export TOMOSTAR_NORCIA_OUT=examples/norcia2016/out TOMOSTAR_NORCIA_DATA=examples/norcia2016/data
export TOMOSTAR_CHIARALUCE2022=/path/to/CAT_files CARANNANTE2025=/path/to/plane_files
export NLL_BIN=/path/to/NonLinLoc/bin HYPODD_BIN=/path/to/HypoDD/bin ATTENTION_DIR=/path/to/AttenTIon
cd validation
dotnet build TomoStar.Validation -c Release
for s in v0*.py; do python3 $s; done
python3 report.py
```

## References

- Beyreuther, M., Barsch, R., Krischer, L., Megies, T., Behr, Y., and Wassermann, J. (2010). ObsPy:
  A Python toolbox for seismology. Seismological Research Letters, 81(3), 530-533.
  https://doi.org/10.1785/gssrl.81.3.530
- Carannante, S., Monachesi, G., Cattaneo, M., Amato, A., and Chiarabba, C. (2013). Deep structure
  and tectonics of the northern-central Apennines as seen by regional-scale tomography and 3-D
  located earthquakes. Journal of Geophysical Research: Solid Earth, 118(10), 5391-5403.
  https://doi.org/10.1002/jgrb.50371
- Carannante, S., Cattaneo, M., and Monachesi, G. (2025). 1D and 3D velocity models of Umbria-Marche
  Region (Central Italy) [data set]. Zenodo. https://doi.org/10.5281/zenodo.16535187
- Chiaraluce, L., Michele, M., Waldhauser, F., Tan, Y. J., Herrmann, M., Spallarossa, D., et al.
  (2022). A comprehensive suite of earthquake catalogues for the 2016-2017 Central Italy seismic
  sequence. Scientific Data, 9, 710. https://doi.org/10.1038/s41597-022-01827-z
- Chiaraluce, L., Michele, M., Waldhauser, F., Tan, Y. J., et al. (2022b). A comprehensive suite of
  earthquake catalogues for the 2016-2017 Central Italy seismic sequence [data set]. NERC EDS
  National Geoscience Data Centre. https://doi.org/10.5285/5afccfe5-142e-4e93-a6cc-55216fa1db06
- Crotwell, H. P., Owens, T. J., and Ritsema, J. (1999). The TauP Toolkit: Flexible seismic
  travel-time and ray-path utilities. Seismological Research Letters, 70(2), 154-160.
  https://doi.org/10.1785/gssrl.70.2.154
- EMERSITO Working Group (2018). Rete sismica del gruppo EMERSITO, sequenza sismica del 2016 in
  Italia Centrale [data set, network XO]. Istituto Nazionale di Geofisica e Vulcanologia (INGV).
  https://doi.org/10.13127/SD/7TXEGDO5X8
- INGV (2005). Rete Sismica Nazionale (RSN) [data set, network IV]. Istituto Nazionale di Geofisica
  e Vulcanologia (INGV). https://doi.org/10.13127/SD/X0FXNH7QFY
- INGV (2023). Seismic Data acquired by Marche Seismic Network (MSN) [data set, network VM]. Istituto
  Nazionale di Geofisica e Vulcanologia (INGV). https://doi.org/10.13127/SD/Z7HOI9U3IX
- INGV, CNR-IGAG and CNR-IDPA (2018). Rete del Centro di Microzonazione Sismica (CentroMZ), sequenza
  sismica del 2016 in Italia Centrale [data set, network 3A]. Istituto Nazionale di Geofisica e
  Vulcanologia (INGV). https://doi.org/10.13127/SD/KU7XM12YY9
- ISIDe Working Group (2007). Italian Seismological Instrumental and Parametric Database (ISIDe)
  [data set]. Istituto Nazionale di Geofisica e Vulcanologia (INGV). https://doi.org/10.13127/ISIDE
- Kennett, B. L. N., Engdahl, E. R., and Buland, R. (1995). Constraints on seismic velocities in the
  Earth from traveltimes. Geophysical Journal International, 122(1), 108-124.
  https://doi.org/10.1111/j.1365-246X.1995.tb03540.x
- Krischer, L., Megies, T., Barsch, R., Beyreuther, M., Lecocq, T., Caudron, C., and Wassermann, J.
  (2015). ObsPy: a bridge for seismology into the scientific Python ecosystem. Computational Science
  and Discovery, 8(1), 014003. https://doi.org/10.1088/1749-4699/8/1/014003
- Lomax, A., Virieux, J., Volant, P., and Berge-Thierry, C. (2000). Probabilistic earthquake location
  in 3D and layered models. In Advances in Seismic Event Location, Kluwer, 101-134.
  https://doi.org/10.1007/978-94-015-9536-0_5
- Maeda, N. (1985). A method for reading and checking phase times in autoprocessing system of seismic
  wave data. Zisin, 38, 365-379. https://doi.org/10.4294/zisin1948.38.3_365
- Marzorati, S., Moretti, M., Margheriti, L., Pondrelli, S., et al. (2023). Seismic Data acquired by
  the SISMIKO Emergency Group, Central Italy 2016, T12 [data set, network 8P]. Istituto Nazionale di
  Geofisica e Vulcanologia (INGV). https://doi.org/10.13127/SD/2PNSQ5UATQ
- MedNet Project Partner Institutions (1990). Mediterranean Very Broadband Seismographic Network
  (MedNet) [data set, network MN]. Istituto Nazionale di Geofisica e Vulcanologia (INGV).
  https://doi.org/10.13127/SD/FBBBTDTD6Q
- Michele, M., Chiaraluce, L., Di Stefano, R., and Waldhauser, F. (2020). Fine-scale structure of the
  2016-2017 Central Italy seismic sequence from data recorded at the Italian National Network.
  Journal of Geophysical Research: Solid Earth, 125(4), e2019JB018440.
  https://doi.org/10.1029/2019JB018440
- Paige, C. C., and Saunders, M. A. (1982). LSQR: An algorithm for sparse linear equations and sparse
  least squares. ACM Transactions on Mathematical Software, 8(1), 43-71.
  https://doi.org/10.1145/355984.355989
- Podvin, P., and Lecomte, I. (1991). Finite difference computation of traveltimes in very
  contrasted velocity models: a massively parallel approach and its associated tools. Geophysical
  Journal International, 105(1), 271-284. https://doi.org/10.1111/j.1365-246X.1991.tb03461.x
- Stachnik, J. C., Abers, G. A., and Christensen, D. H. (2004). Seismic attenuation and mantle wedge
  temperatures in the Alaska subduction zone. Journal of Geophysical Research, 109, B10304.
  https://doi.org/10.1029/2004JB003018
- Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., et al. (2020). SciPy 1.0:
  fundamental algorithms for scientific computing in Python. Nature Methods, 17(3), 261-272.
  https://doi.org/10.1038/s41592-019-0686-2
- Waldhauser, F. (2001). hypoDD: A program to compute double-difference hypocenter locations. U.S.
  Geological Survey Open-File Report 01-113. https://doi.org/10.3133/ofr01113
- Waldhauser, F., and Ellsworth, W. L. (2000). A double-difference earthquake location algorithm:
  method and application to the northern Hayward fault, California. Bulletin of the Seismological
  Society of America, 90(6), 1353-1368. https://doi.org/10.1785/0120000006
- Waldhauser, F., Michele, M., Chiaraluce, L., Di Stefano, R., and Schaff, D. P. (2021). Fault planes,
  fault zone structure and detachment fragmentation resolved with high-precision aftershock
  locations of the 2016-2017 central Italy sequence. Geophysical Research Letters, 48(16),
  e2021GL092918. https://doi.org/10.1029/2021GL092918
- Wei, S. S., and Wiens, D. A. (2018). P-wave attenuation structure of the Lau back-arc basin and
  implications for mantle wedge processes. Earth and Planetary Science Letters, 502, 187-199.
  https://doi.org/10.1016/j.epsl.2018.09.005
- White, M. C. A., Fang, H., Nakata, N., and Ben-Zion, Y. (2020). PyKonal: A Python package for
  solving the eikonal equation in spherical and Cartesian coordinates using the fast marching
  method. Seismological Research Letters, 91(4), 2378-2389. https://doi.org/10.1785/0220190318
- Wölbern, I., Rümpker, G., and Leva, C. (2020). FOSA [data set, network 5M]. GFZ Data Services.
  https://doi.org/10.14470/0Z7560909466
