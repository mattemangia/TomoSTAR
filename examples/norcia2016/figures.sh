#!/usr/bin/env bash
# Figures of the Norcia 2016 example (the ones in the README), from the results of norcia.tomo.
# Needs Python 3 with numpy and matplotlib.
#
#     bash figures.sh [RESULTS_FOLDER] [FIGURE_FOLDER]      (defaults: ./out and ./figures)
set -euo pipefail
OUT=${1:-out}
FIG=${2:-figures}
PLOT="python3 $(dirname "$0")/../../tools/plot_tomostar.py"
mkdir -p "$FIG"
EV="$OUT/vel/events_relocated.csv"
ST="$OUT/vel/stations.csv"
DEPTHS=2,6,10

# Vp, Vp/Vs and their changes from the minimum 1-D model, masked where DWS < 50 km of ray.
$PLOT slices "$OUT/vel/volumes/vel_dVp.qvol" --dws "$OUT/vel/volumes/vel_DWS_P.qvol" --min-dws 50 \
  --depths $DEPTHS --events "$EV" --stations "$ST" --limits=-8,8 --cmap RdBu --label "dVp/Vp (%)" \
  --title "Vp change from the minimum 1-D model, Amatrice-Visso-Norcia 2016-2017" --out "$FIG/norcia_dvp.png"
$PLOT slices "$OUT/vel/volumes/vel_VpVs.qvol" --dws "$OUT/vel/volumes/vel_DWS_S.qvol" --min-dws 50 \
  --depths $DEPTHS --events "$EV" --stations "$ST" --limits 1.65,2.05 --cmap RdYlBu_r --label "Vp/Vs" \
  --title "Vp/Vs" --out "$FIG/norcia_vpvs.png"

# A section across the fault system (SW to NE through Norcia), events within 3 km.
$PLOT section "$OUT/vel/volumes/vel_Vp.qvol" --dws "$OUT/vel/volumes/vel_DWS_P.qvol" --min-dws 50 \
  --from 12.87,42.735 --to 13.43,42.925 --width-km 4 --max-depth 16 --events "$EV" --cmap viridis_r --label "Vp (km/s)" \
  --title "Vp, WSW to ENE section through Norcia, across the fault system" --out "$FIG/norcia_section_vp.png"
$PLOT section "$OUT/vel/volumes/vel_VpVs.qvol" --dws "$OUT/vel/volumes/vel_DWS_S.qvol" --min-dws 50 \
  --from 12.87,42.735 --to 13.43,42.925 --width-km 4 --max-depth 16 --events "$EV" --limits 1.65,2.05 --cmap RdYlBu_r --label "Vp/Vs" \
  --title "Vp/Vs, WSW to ENE section through Norcia, across the fault system" --out "$FIG/norcia_section_vpvs.png"

# Adaptive grid: the model and the size of its cells.
$PLOT slices "$OUT/vel_adaptive/volumes/vel_adaptive_dVp.qvol" --dws "$OUT/vel_adaptive/volumes/vel_adaptive_DWS_P.qvol" --min-dws 50 \
  --depths $DEPTHS --events "$EV" --stations "$ST" --limits=-8,8 --cmap RdBu --label "dVp/Vp (%)" \
  --title "Vp change on the adaptive (octree) grid" --out "$FIG/norcia_adaptive_dvp.png"
$PLOT slices "$OUT/vel_adaptive/volumes/vel_adaptive_CellSize.qvol" \
  --depths $DEPTHS --stations "$ST" --cmap cividis --label "cell size (km)" \
  --title "Size of the adaptive cells (fine where the rays are dense)" --out "$FIG/norcia_adaptive_cells.png"

# Checkerboard test on the real geometry.
$PLOT pair "$OUT/checker_vp/volumes/checker_vp_dVp_true.qvol" "$OUT/checker_vp/volumes/checker_vp_dVp_recovered.qvol" \
  --dws "$OUT/checker_vp/dws.qvol" --min-dws 50 --depths $DEPTHS --stations "$ST" --label "dVp/Vp (%)" \
  --title "Checkerboard test (12 x 12 x 6 km cells, 5 %) on the real source-receiver geometry" --out "$FIG/norcia_checkerboard.png"

# Qp from t*, when the Q tomography was run.
if [ -f "$OUT/q/volumes/q_Qp.qvol" ]; then
  $PLOT slices "$OUT/q/volumes/q_Qp.qvol" --dws "$OUT/q/volumes/q_DWS_Q.qvol" --min-dws 50 \
    --depths $DEPTHS --events "$EV" --stations "$ST" --cmap magma --label "Qp" \
    --title "Qp from t* of the M >= 3.5 events" --out "$FIG/norcia_qp.png"
fi

# The L-curve figure written by TomoSTAR itself.
cp "$OUT/lcurve/lcurve.svg" "$FIG/norcia_lcurve.svg"
echo "Figures in $FIG."
