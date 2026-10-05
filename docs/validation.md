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
| Case | TomoSTAR method | Reference | Metric | Value | Limit | Result |
|---|---|---|---|---|---|---|
| eikonal | fast marching, ak135 first arrivals (P and S) | TauP (ObsPy) | max relative error | 0.98 % | 5.00 % | pass |
| eikonal | fast marching, real 3-D Vp and Vs model (central Italy 2016-2017), station to hypocentre, beyond 5 km | PyKonal (grid refined 4 times) | max relative error | 2.29 % | 5.00 % | pass |
| rays | ray paths traced back along the travel-time gradient, real 3-D model | PyKonal | max deviation / ray length | 3.15 % | 5.00 % | pass |
| rays | travel time integrated along the ray (row of G), real 3-D model | eikonal travel time | max relative error | 2.13 % | 5.00 % | pass |
| lsqr | damped LSQR on the real tomography problem (central Italy 2016-2017) | SciPy lsqr | max \|\|x - x_ref\|\| / \|\|x_ref\|\| | 2.5e-09 | 5.00 % | pass |
| qtomography | Qp tomography forward operator and model (t* = sum L s q + station term), real rays and t* | independent path integration of the written model along the stored rays | 95th percentile of relative t* difference | 1.29 % | 5.00 % | pass |
| location | absolute location (grid search + Levenberg-Marquardt), INGV picks, Carannante et al. 2013 model | NonLinLoc (Lomax et al. 2000), same picks and model | 95th percentile of \|hypocentre difference\| / median hypocentral distance | 1.02 % | 5.00 % | pass |
| location | absolute location, INGV picks, Carannante et al. 2013 model, no station corrections | published catalogue CAT1 of Chiaraluce et al. 2022 (NonLinLoc with station corrections) | 95th percentile of \|hypocentre difference\| / median hypocentral distance | 12.11 % |  | information |
| double_difference | double-difference relocation, catalogue differential times of the INGV picks | HypoDD 2.1b (Waldhauser 2001), same picks, pairs, weights and model | 95th percentile of \|hypocentre difference\| / median hypocentral distance (cluster offsets removed) | 4.07 % | 5.00 % | pass |
| double_difference | double-difference relocation, catalogue differential times of the INGV picks | published catalogue CAT2 (Michele et al. 2020; HypoDD with cross-correlation delays) | 95th percentile of \|hypocentre difference\| / median hypocentral distance (mean offset removed) | 6.76 % |  | information |
| tstar | joint multitaper t* inversion of an event (Brune source, shared ln Omega0, alpha 0.27, t* >= 0) | AttenTIon inversion() (Stachnik et al. 2004; Wei and Wiens 2018), same spectra and corner frequency | max over events of \|\|t*_TomoSTAR - t*_AttenTIon\|\| / \|\|t*_AttenTIon\|\| | 1.3e-08 | 5.00 % | pass |
| tstar | joint multitaper t* inversion with its own corner-frequency search | AttenTIon bestfc() + inversion(), same spectra, search records and corner range | median over events of \|\|t*_TomoSTAR - t*_AttenTIon\|\| / \|\|t*_AttenTIon\|\| | 6.69 % |  | information |
| tstar | joint multitaper t* inversion with its corner-frequency search (least squares) | AttenTIon d, G and nnls on its corner grid with the least-squares criterion, same spectra, search records and range | median over events of \|\|t*_TomoSTAR - t*_AttenTIon\|\| / \|\|t*_AttenTIon\|\| | 0.41 % | 5.00 % | pass |
| picker | automatic P picks (STA/LTA trigger, AIC onset) on real waveforms | INGV analyst picks of the same records, most precise class (0.1 s) | 95th percentile of \|t_auto - t_analyst\| / travel time | 4.11 % | 5.00 % | pass |
| picker | automatic P picks (STA/LTA trigger, AIC onset) on real waveforms | INGV analyst picks of the same records, all classes (0.1 to 1 s) | 95th percentile of \|t_auto - t_analyst\| / travel time | 6.97 % |  | information |
| picker | automatic S picks (STA/LTA trigger, AIC onset) on real waveforms | INGV analyst picks of the same records, most precise class (0.1 s) | 95th percentile of \|t_auto - t_analyst\| / travel time | 5.57 % | 5.00 % | FAIL |
| picker | automatic S picks (STA/LTA trigger, AIC onset) on real waveforms | INGV analyst picks of the same records, all classes (0.1 to 1 s) | 95th percentile of \|t_auto - t_analyst\| / travel time | 8.38 % |  | information |
| signal | instrument response removal to velocity | ObsPy remove_response | max relative L2 difference | 0.14 % | 5.00 % | pass |
| signal | Butterworth band-pass 1-15 Hz, orders 2 and 4, causal and zero phase | SciPy butter + sosfilt | max relative L2 difference | 6.0e-08 | 5.00 % | pass |
| signal | Slepian (DPSS) tapers | SciPy dpss | max L2 difference of the unit-norm tapers | 5.2e-13 | 5.00 % | pass |
| signal | recursive STA/LTA | ObsPy recursive_sta_lta | max relative difference after 8 LTA | 9.5e-08 | 5.00 % | pass |
| signal | AIC function of the onset picker (Maeda 1985) | ObsPy aic_simple (aligned by one sample) | max relative difference | 4.3e-16 | 5.00 % | pass |
| formats | miniSEED reader (samples, start times) | ObsPy read | max relative sample difference | 0 | 5.00 % | pass |
| formats | QuakeML reader (origins) | ObsPy read_events | max difference (degrees, s, and depth / 100 km) | 0 | 5.00 % | pass |
| formats | QuakeML reader (associated P and S picks) | ObsPy read_events | max time difference (s) | 0 | 5.00 % | pass |
| formats | StationXML reader (station coordinates) | ObsPy read_inventory | max difference (degrees, km) | 0 | 5.00 % | pass |
| formats | instrument response amplitude (stage gains, poles and zeros, digital filters) | ObsPy evalresp (all stages) | max relative error in the pass band | 0.83 % | 5.00 % | pass |
| published_model | Vp of the central Italy 2016-2017 example (TomoSTAR) | published 3-D Vp of Carannante et al. 2013 (SIMULPS14, other data) | median \|Vp - Vp_published\| / Vp_published | 2.80 % |  | information |
<!-- results end -->

## Notes on the comparisons

<!-- notes begin -->
### Defects found and fixed

The comparisons found five defects in TomoSTAR, all fixed before the results above were produced:

1. **Double difference.** A step was judged by the misfit of the absolute times alone and halved
   when that misfit rose, although the step had been computed to lower the differential misfit; the
   events moved two thirds as much as with HypoDD. The step is now judged by the objective of the
   linear system (absolute and differential rows with their weights).
2. **Picker window.** The search window grew with the time from the start of the trace instead of
   the origin time, and the first STA/LTA trigger in it was taken, often the coda or the arrival of
   another earthquake during the sequence (a third of the P picks of classes 1 to 3 more than a
   second early). Of the arrivals at least half as strong as the strongest, the one nearest the
   predicted time is now taken.
3. **Stage gain frequency.** For channels whose stage gain is declared at another frequency than the
   normalisation of the poles and zeros (IV.SAP2 and IV.ATPI in the INGV metadata), the response was
   scaled to the sensitivity at a frequency where it did not hold (25 times too high for IV.SAP2).
   The response is now built from the stage gains and normalisation factors as evalresp builds it.
4. **Digital filters.** The decimation filters were assumed flat in their passband; the binomial
   filters of the IV digitisers of 2008-2009 halve the amplitude at 20 Hz on a 100 Hz channel. FIR
   and coefficient stages are now part of the response.
5. **Uncorrected filter delay.** An asymmetric FIR stage without a time correction (IV.GIGS in 2016)
   shifts the record by three samples; its phase is now kept (symmetric filters stay zero phase, as
   in evalresp).

Defects 3 to 5 change the t\* of the example by 0.1 ms (median) and 2 ms (95th percentile); the
t\* and Qp results of the README were recomputed with the corrected responses.

### Absolute location

TomoSTAR and NonLinLoc locate the same 2991 events from the same 159,736 INGV picks with the same
pick uncertainties (0.1, 0.3 and 0.6 s by the bulletin's quality classes), in the same 1-D gradient
model, without station corrections, both by weighted least squares (NonLinLoc: oct-tree search with
the GAU_ANALYTIC likelihood; TomoSTAR: grid search and Levenberg-Marquardt). NonLinLoc's travel
times come from its finite-difference solver (Podvin and Lecomte, 1991) on a 0.1 km grid in a
spherical projection; TomoSTAR's from fast marching on a 0.01 degree by 0.5 km grid. The median
hypocentre difference is 150 m (32 m horizontally), with TomoSTAR 120 m shallower on average, which
is within the 600 m of the published vertical uncertainty. The two programs report different RMS
values for the same solutions because NonLinLoc's is weighted by the pick uncertainties and
TomoSTAR's is not.

Both programs differ from the published CAT1 hypocentres by about 1 km (median), by the same amount:
CAT1 used station corrections and the picks of 24 temporary stations that are not public.

### Double difference

HypoDD and TomoSTAR relocate the same events from the same starting hypocentres with the catalogue
differential times of the same picks, the same pair selection (ph2dt: 10 km, 10 neighbours, 8 to 50
links) and the same three iteration sets (no cutoff; 6 MADs and 10 km; 5 MADs and 6 km). HypoDD uses
1 km layers of the gradient model, as Michele et al. (2020) did; TomoSTAR uses the gradient model on
its grid and keeps the absolute times at a token weight (0.01) where HypoDD fixes the cluster
centroid. HypoDD moves the events by 0.88 km (median) and TomoSTAR by 0.84 km; the relocations
differ by 373 m (median). Against the published CAT2, which adds 4.4 million cross-correlation
delays, both programs differ by the same amount (6.8 and 6.6 % at the 95th percentile).

### t\*

AttenTIon's own functions (`inversion`, `bestfc`, `buildd`, `buildG`) are given the multitaper
spectra TomoSTAR measures on the real waveforms (243 events with at least five records, 2,786
records). At the same corner
frequency the two inversions agree to rounding. With each program's own corner search the t\* differ
by 6.7 % (median): AttenTIon minimises the residual norm divided by the sum of the data, which favours
low corners (TomoSTAR's corners are 12 % higher); with AttenTIon's corner grid and the least-squares
criterion the corners coincide and t\* agree to 0.4 %.

### Picker

The automatic picks are made without the analyst picks, around the times the example's relocated
3-D model predicts (window of 1 s plus 3 % of the travel time). Against the analysts' most precise
picks (0.1 s), 95 % of the P picks are within 4.1 % of the travel time (median error 30 ms); the S
picks are within 5.6 %, just above the limit, on 130 picks. S picking is the weakest step of the
pipeline and its picks enter the inversions with larger uncertainties (1.5 times those of P). Part
of the disagreement is on the analysts' side: some analyst picks were made on accelerometers, where
an emergent onset visible on the broadband sensor is below the noise (IV.NRCA, event 10740261: a
weak onset 0.83 s before the analyst pick).

### Forward problem, LSQR and Q

The eikonal solver is compared with TauP in ak135 and with PyKonal in the real 3-D model of the
example, from the eight stations with most picks to every relocated earthquake. LSQR is compared
with SciPy on the matrix of the last iteration of the example's inversion with its real residuals.
The Q tomography's forward operator is recomputed by an independent integration of the written Qp
and Vp models along the stored rays. The inverted Vp model of the example is compared, for
information, with the published model of Carannante et al. (2013), obtained from other data: the
median difference is 2.8 % over the 42 nodes resolved by both, larger near the surface (6.4 %) than
at 4 and 8 km (2.0 and 2.8 %).

### Signal processing and formats

On the example's data (StationXML of 209 stations with 1002 channel epochs, 93 miniSEED traces of
which 90 have a response, QuakeML of 42 events with 5537 picks), the instrument response agrees with
evalresp to 0.8 % and the response removal with ObsPy to 0.14 %
(with the same demean, detrend and 5 % cosine taper), the readers agree exactly with ObsPy, and the
filters, Slepian tapers, STA/LTA and AIC agree with SciPy and ObsPy to rounding (the recursive
STA/LTA after its start-up of 8 LTA, the AIC with ObsPy's split one sample later).
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
