#!/usr/bin/env bash
# Downloads the real data of the Norcia 2016 example from the INGV FDSN web services
# (https://webservices.ingv.it): the M >= 2.5 earthquakes of the 2016-2017 Amatrice-Visso-Norcia
# sequence with the analyst picks of the INGV bulletin (QuakeML, one file per event, since the
# INGV event service returns the picks of one event at a time), the stations of the area
# (StationXML with responses) and, for the M >= 3.5 events, the vertical and horizontal
# waveforms used to measure t* (miniSEED, one folder per event).
#
#     bash download.sh [DATA_FOLDER]        (default: ./data)
#
# The data are distributed by INGV under CC BY 4.0; cite the network (IV, doi:10.13127/SD/X0FXNH7QFY)
# and the bulletin (ISIDe working group, doi:10.13127/ISIDE) when you use them.
set -euo pipefail
OUT=${1:-data}
BASE=https://webservices.ingv.it/fdsnws
BOX="minlat=42.45&maxlat=43.1&minlon=12.9&maxlon=13.55"
mkdir -p "$OUT/quakeml" "$OUT/waveforms"

# Events (text list), then the QuakeML with arrivals of each one, eight at a time.
curl -sS --compressed -o "$OUT/events.txt" \
  "$BASE/event/1/query?starttime=2016-08-24&endtime=2017-03-01&$BOX&maxdepth=20&minmag=2.5&format=text"
tail -n +2 "$OUT/events.txt" | cut -d'|' -f1 | \
  xargs -P 8 -I{} sh -c "[ -s '$OUT/quakeml/{}.xml' ] || curl -sS --compressed --retry 4 -o '$OUT/quakeml/{}.xml' '$BASE/event/1/query?eventid={}&includearrivals=true&format=xml'"

# Stations of every network with data in the area during the sequence, with instrument responses.
curl -sS --compressed --retry 4 -o "$OUT/stations.xml" \
  "https://federator.orfeus-eu.org/fdsnws/station/1/query?starttime=2016-08-24&endtime=2017-03-01&minlat=41.85&maxlat=43.7&minlon=12.1&maxlon=14.35&channel=HH?,EH?,HN?&level=response"

# Waveforms of the M >= 3.5 events, 20 s before to 60 s after the origin, broadband and short-period
# channels of the IV stations inside the study area (a POST request to the dataselect service).
curl -sS --retry 4 -o "$OUT/iv_stations.txt" \
  "$BASE/station/1/query?network=IV&starttime=2016-08-24&endtime=2017-03-01&minlat=42.3&maxlat=43.25&minlon=12.7&maxlon=13.75&format=text"
# One POST per event, six events at a time (the service answers each request in about a minute).
fetch_event() {
  id=$1; time=$2
  [ -s "$OUT/waveforms/$id/$id.mseed" ] && return 0
  mkdir -p "$OUT/waveforms/$id"
  start=$(date -u -d "${time}Z - 20 seconds" +%Y-%m-%dT%H:%M:%S)
  end=$(date -u -d "${time}Z + 60 seconds" +%Y-%m-%dT%H:%M:%S)
  tail -n +2 "$OUT/iv_stations.txt" | cut -d'|' -f2 | sort -u | \
    awk -v s="$start" -v e="$end" '{ print "IV " $1 " * HH?,EH? " s " " e }' > "$OUT/waveforms/$id/request.txt"
  curl -sS --retry 4 -o "$OUT/waveforms/$id/$id.mseed.part" --data-binary @"$OUT/waveforms/$id/request.txt" \
    "$BASE/dataselect/1/query" && mv "$OUT/waveforms/$id/$id.mseed.part" "$OUT/waveforms/$id/$id.mseed"
  rm -f "$OUT/waveforms/$id/request.txt" "$OUT/waveforms/$id/$id.mseed.part"
}
export -f fetch_event
export OUT BASE
awk -F'|' 'NR > 1 && $11 >= 3.5 { print $1, $2 }' "$OUT/events.txt" | xargs -P 6 -n 2 bash -c 'fetch_event "$0" "$1"'
touch "$OUT/waveforms/.complete"
echo "Data in $OUT: $(ls "$OUT/quakeml" | wc -l) events, stations.xml, waveforms of $(ls "$OUT/waveforms" | wc -l) events."
